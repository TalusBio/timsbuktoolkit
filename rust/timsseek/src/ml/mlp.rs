//! Small in-process MLP for rescoring.
//!
//! Tensors are flat row-major `f32` buffers. Training and scoring consume
//! regenerated batches, so a full candidate or fold matrix is never retained.

// Index-based loops are deliberate throughout. These are numeric kernels where
// one index walks several parallel arrays at once (`w`/`gw`, `mean`/`var`/
// `inv_std`) and where the inner AXPY loops are shaped for vectorization.
// Rewriting them as iterator zips needs 3- and 4-way nests and sinks the
// codegen — the same tradeoff `ml/lda.rs` documents at its own hot loops.
#![allow(clippy::needless_range_loop)]

use matrixmultiply::sgemm;
use rand::rngs::StdRng;
use rand::{
    Rng as _,
    SeedableRng,
};

/// Checked, safe boundary around [`sgemm`] for the layouts used by this MLP.
/// Inputs may be row-major or transposed views via positive strides; output is
/// always contiguous row-major. Separate shared/mutable borrows rule out output
/// aliasing either input at every safe call site.
#[allow(clippy::too_many_arguments)]
fn gemm_f32(
    m: usize,
    k: usize,
    n: usize,
    a: &[f32],
    a_row_stride: usize,
    a_col_stride: usize,
    b: &[f32],
    b_row_stride: usize,
    b_col_stride: usize,
    beta: f32,
    c: &mut [f32],
) {
    fn required_len(rows: usize, cols: usize, row_stride: usize, col_stride: usize) -> usize {
        if rows == 0 || cols == 0 {
            return 0;
        }
        (rows - 1)
            .checked_mul(row_stride)
            .and_then(|last_row| {
                (cols - 1)
                    .checked_mul(col_stride)
                    .and_then(|last_col| last_row.checked_add(last_col))
            })
            .and_then(|last| last.checked_add(1))
            .expect("GEMM matrix extent overflows usize")
    }

    assert!(
        a.len() >= required_len(m, k, a_row_stride, a_col_stride),
        "GEMM left input is smaller than its declared strided extent"
    );
    assert!(
        b.len() >= required_len(k, n, b_row_stride, b_col_stride),
        "GEMM right input is smaller than its declared strided extent"
    );
    assert!(
        c.len()
            >= m.checked_mul(n)
                .expect("GEMM output extent overflows usize"),
        "GEMM output is smaller than its declared row-major extent"
    );
    let ars = isize::try_from(a_row_stride).expect("GEMM left row stride exceeds isize");
    let acs = isize::try_from(a_col_stride).expect("GEMM left column stride exceeds isize");
    let brs = isize::try_from(b_row_stride).expect("GEMM right row stride exceeds isize");
    let bcs = isize::try_from(b_col_stride).expect("GEMM right column stride exceeds isize");
    let crs = isize::try_from(n).expect("GEMM output row stride exceeds isize");

    // SAFETY: the extent checks prove every strided input access and contiguous
    // output access is in bounds. Safe borrowing prevents C from aliasing A or
    // B, and C's fixed (n, 1) strides cannot alias distinct output elements.
    unsafe {
        sgemm(
            m,
            k,
            n,
            1.0,
            a.as_ptr(),
            ars,
            acs,
            b.as_ptr(),
            brs,
            bcs,
            beta,
            c.as_mut_ptr(),
            crs,
            1,
        );
    }
}

// ---------------------------------------------------------------- tensor

/// Flat row-major `[rows][cols]` buffer.
#[derive(Clone, Default)]
pub struct Tensor {
    data: Vec<f32>,
    rows: usize,
    cols: usize,
}

impl Tensor {
    pub fn new(rows: usize, cols: usize) -> Self {
        Tensor {
            data: vec![0.0; rows * cols],
            rows,
            cols,
        }
    }

    /// Grow to fit `rows x cols`. Capacity is never given back, so after the
    /// first few steps this stops touching the allocator entirely.
    pub fn reshape(&mut self, rows: usize, cols: usize) {
        if self.data.len() < rows * cols {
            self.data.resize(rows * cols, 0.0);
        }
        self.rows = rows;
        self.cols = cols;
    }

    #[cfg(test)]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f32] {
        &self.data[i * self.cols..i * self.cols + self.cols]
    }

    #[inline]
    pub fn row_mut(&mut self, i: usize) -> &mut [f32] {
        let c = self.cols;
        &mut self.data[i * c..i * c + c]
    }

    #[inline]
    fn live(&self) -> &[f32] {
        &self.data[..self.rows * self.cols]
    }

    #[inline]
    fn live_mut(&mut self) -> &mut [f32] {
        let n = self.rows * self.cols;
        &mut self.data[..n]
    }
}

// ---------------------------------------------------------------- layer trait

/// `(param, grad, apply_weight_decay)`. Biases and norm parameters return
/// `false` — decaying them toward zero fights what those layers are for.
type ParamSlots<'a> = Vec<(&'a mut [f32], &'a mut [f32], bool)>;

pub trait Layer: Send {
    fn out_dim(&self) -> usize;
    fn forward(&mut self, x: &Tensor, out: &mut Tensor);
    /// `x` = layer input, `y` = layer output, `gy` = dL/dy, writes dL/dx into `gx`.
    fn backward(&mut self, x: &Tensor, y: &Tensor, gy: &Tensor, gx: &mut Tensor);
    fn params_and_grads(&mut self) -> ParamSlots<'_> {
        Vec::new()
    }
    /// `sum_o |w[j, o]|` per input, for the layers that have weights over
    /// inputs. Only `Linear` answers; everything else returns `None` and the
    /// model takes the first `Some`.
    ///
    /// This exists so the trait needs no `as_any` + downcast just to reach the
    /// first layer's weights for feature importance.
    fn abs_in_weights(&self) -> Option<Vec<f32>> {
        None
    }
}

// ---------------------------------------------------------------- linear

pub struct Linear {
    n_in: usize,
    n_out: usize,
    /// Flat, row-major, stride `n_out`: `w[i * n_out + o]`.
    w: Vec<f32>,
    b: Vec<f32>,
    gw: Vec<f32>,
    gb: Vec<f32>,
}

impl Linear {
    /// He-normal weights (`std = sqrt(2/n_in)`), Torch-style uniform biases
    /// (`U(-1/sqrt(n_in), 1/sqrt(n_in))`).
    ///
    /// The biases are NOT initialized to zero, which is the more common
    /// convention. Zero biases put every unit's ReLU boundary through the
    /// origin, and standardized inputs are centered on the origin, so at init
    /// the whole first layer hinges on the exact center of the data. A unit
    /// that lands with all its inputs on the negative side outputs zero, and
    /// a dead ReLU has EXACTLY zero gradient — an absorbing state that no
    /// learning rate, epoch count, or minibatch noise recovers from. A random
    /// bias moves those boundaries away from the shared input origin.
    pub fn new(n_in: usize, n_out: usize, rng: &mut StdRng) -> Self {
        let scale = (2.0 / n_in as f32).sqrt(); // He init
        let bound = 1.0 / (n_in as f32).sqrt();
        Linear {
            n_in,
            n_out,
            w: (0..n_in * n_out).map(|_| normal(rng) * scale).collect(),
            b: (0..n_out)
                .map(|_| rng.random_range(-bound..bound))
                .collect(),
            gw: vec![0.0; n_in * n_out],
            gb: vec![0.0; n_out],
        }
    }
}

impl Layer for Linear {
    fn out_dim(&self) -> usize {
        self.n_out
    }

    fn forward(&mut self, x: &Tensor, out: &mut Tensor) {
        let (n, ni, no) = (x.rows, self.n_in, self.n_out);
        out.reshape(n, no);
        for row in out.data[..n * no].chunks_exact_mut(no) {
            row.copy_from_slice(&self.b);
        }
        gemm_f32(
            n,
            ni,
            no,
            &x.data,
            ni,
            1,
            &self.w,
            no,
            1,
            1.0,
            &mut out.data,
        );
    }

    fn backward(&mut self, x: &Tensor, _y: &Tensor, gy: &Tensor, gx: &mut Tensor) {
        let (n, ni, no) = (gy.rows, self.n_in, self.n_out);
        gx.reshape(n, ni);
        for gr in gy.data[..n * no].chunks_exact(no) {
            for (gb, &value) in self.gb.iter_mut().zip(gr) {
                *gb += value;
            }
        }
        // dW = X^T dY; beta=1 honors gradient accumulation.
        gemm_f32(
            ni,
            n,
            no,
            &x.data,
            1,
            ni,
            &gy.data,
            no,
            1,
            1.0,
            &mut self.gw,
        );
        // dX = dY W^T; beta=0 initializes every live output element.
        gemm_f32(
            n,
            no,
            ni,
            &gy.data,
            no,
            1,
            &self.w,
            1,
            no,
            0.0,
            &mut gx.data,
        );
    }

    fn params_and_grads(&mut self) -> ParamSlots<'_> {
        vec![
            (self.w.as_mut_slice(), self.gw.as_mut_slice(), true),
            (self.b.as_mut_slice(), self.gb.as_mut_slice(), false),
        ]
    }

    /// With standardized inputs the first layer's weights are already on a
    /// comparable scale, which is what makes this the MLP's analogue of LDA's
    /// `|coef|`.
    fn abs_in_weights(&self) -> Option<Vec<f32>> {
        Some(
            (0..self.n_in)
                .map(|i| {
                    self.w[i * self.n_out..(i + 1) * self.n_out]
                        .iter()
                        .map(|v| v.abs())
                        .sum()
                })
                .collect(),
        )
    }
}

// ---------------------------------------------------------------- relu

/// Default negative slope. Small enough to be a rectifier in every way that
/// matters, nonzero so a saturated unit still receives gradient.
pub const LEAKY_SLOPE: f32 = 0.01;

/// Rectifier with an optional negative slope.
///
/// The default architecture uses [`Relu::leaky`], not the pure form, because
/// a pure ReLU's zero region is an ABSORBING state: its gradient is exactly
/// zero, so a unit that dies contributes nothing and receives nothing, forever.
/// That is not a slow-learning problem that a smaller learning rate or more
/// epochs or minibatch noise can fix — there is no gradient to be noisy about.
pub struct Relu {
    dim: usize,
    slope: f32,
}

impl Relu {
    pub fn leaky(dim: usize, slope: f32) -> Self {
        Relu { dim, slope }
    }
}

impl Layer for Relu {
    fn out_dim(&self) -> usize {
        self.dim
    }

    fn forward(&mut self, x: &Tensor, out: &mut Tensor) {
        out.reshape(x.rows, x.cols);
        let (src, dst) = (x.live(), out.live_mut());
        for k in 0..src.len() {
            dst[k] = if src[k] > 0.0 {
                src[k]
            } else {
                self.slope * src[k]
            };
        }
    }

    /// No mask needed: the output already encodes it. With a positive slope
    /// `y` keeps the sign of `x`, so `y > 0` still recovers the active set.
    fn backward(&mut self, _x: &Tensor, y: &Tensor, gy: &Tensor, gx: &mut Tensor) {
        gx.reshape(gy.rows, gy.cols);
        let (yv, gv, dst) = (y.live(), gy.live(), gx.live_mut());
        for k in 0..dst.len() {
            dst[k] = if yv[k] > 0.0 {
                gv[k]
            } else {
                self.slope * gv[k]
            };
        }
    }
}

// ---------------------------------------------------------------- loss

#[inline]
fn sigmoid(z: f32) -> f32 {
    if z >= 0.0 {
        1.0 / (1.0 + (-z).exp())
    } else {
        let e = z.exp();
        e / (1.0 + e)
    }
}

#[inline]
fn bce_from_logit(z: f32, y: f32) -> f32 {
    z.max(0.0) - z * y + (-z.abs()).exp().ln_1p()
}

/// Weighted BCE-with-logits, writing dL/dlogit into `grad`.
fn loss_and_grad(logits: &Tensor, y: &[f32], w: &[f32], grad: &mut Tensor) -> f32 {
    let n = logits.rows;
    debug_assert_eq!(logits.cols, 1);
    debug_assert_eq!(y.len(), n);
    debug_assert_eq!(w.len(), n);
    grad.reshape(n, 1);

    let wsum: f32 = w.iter().sum();
    let inv_w = if wsum > 0.0 { 1.0 / wsum } else { 0.0 };
    let mut total = 0.0;
    for b in 0..n {
        let z = logits.row(b)[0];
        total += w[b] * bce_from_logit(z, y[b]);
        grad.row_mut(b)[0] = w[b] * inv_w * (sigmoid(z) - y[b]);
    }
    total * inv_w
}

// ---------------------------------------------------------------- optimizer

/// Adam / AdamW. `decoupled_wd = true` gives AdamW.
pub struct Adam {
    lr: f32,
    beta1: f32,
    beta2: f32,
    eps: f32,
    weight_decay: f32,
    decoupled_wd: bool,
    t: i32,
    m: Vec<Vec<f32>>,
    v: Vec<Vec<f32>>,
}

impl Adam {
    pub fn new(lr: f32) -> Self {
        Adam {
            lr,
            beta1: 0.9,
            beta2: 0.999,
            eps: 1e-8,
            weight_decay: 0.0,
            decoupled_wd: true,
            t: 0,
            m: Vec::new(),
            v: Vec::new(),
        }
    }

    pub fn with_weight_decay(mut self, wd: f32) -> Self {
        self.weight_decay = wd;
        self
    }

    fn state_for(&mut self, slot: usize, len: usize) {
        if self.m.len() <= slot {
            self.m.resize(slot + 1, Vec::new());
            self.v.resize(slot + 1, Vec::new());
        }
        if self.m[slot].len() != len {
            self.m[slot] = vec![0.0; len];
            self.v[slot] = vec![0.0; len];
        }
    }

    pub fn step(&mut self, params: ParamSlots<'_>) {
        self.t += 1; // once per step, not per tensor
        let bc1 = 1.0 - self.beta1.powi(self.t);
        let bc2 = 1.0 - self.beta2.powi(self.t);

        for (slot, (p, g, decay)) in params.into_iter().enumerate() {
            self.state_for(slot, p.len());
            let m = &mut self.m[slot];
            let v = &mut self.v[slot];
            let wd = if decay { self.weight_decay } else { 0.0 };

            for i in 0..p.len() {
                let mut grad = g[i];
                if wd != 0.0 && !self.decoupled_wd {
                    grad += wd * p[i]; // classic Adam: L2 inside the moments
                }
                m[i] = self.beta1 * m[i] + (1.0 - self.beta1) * grad;
                v[i] = self.beta2 * v[i] + (1.0 - self.beta2) * grad * grad;
                let m_hat = m[i] / bc1;
                let v_hat = v[i] / bc2;
                if wd != 0.0 && self.decoupled_wd {
                    p[i] -= self.lr * wd * p[i]; // AdamW: decoupled
                }
                p[i] -= self.lr * m_hat / (v_hat.sqrt() + self.eps);
                g[i] = 0.0;
            }
        }
    }
}

// ---------------------------------------------------------------- model

pub struct Mlp {
    layers: Vec<Box<dyn Layer>>,
    acts: Vec<Tensor>,  // acts[k] = output of layer k
    grads: Vec<Tensor>, // grads[k] = dL/d(acts[k])
    scratch: Tensor,    // sink for dL/dx of layer 0
    /// Throwaway gradient sink for [`Self::eval_loss`], persistent for the same
    /// reason every other buffer here is: with early stopping on, `eval_loss`
    /// runs once per epoch, so a local `Tensor` would be one
    /// `val_rows * 1` alloc/free per epoch for a value that is never read. Kept
    /// SEPARATE from `grads` so an eval can never leave a gradient behind that a
    /// later backward pass might read.
    eval_sink: Tensor,
}

impl Mlp {
    pub fn new(layers: Vec<Box<dyn Layer>>) -> Self {
        let acts = layers.iter().map(|l| Tensor::new(0, l.out_dim())).collect();
        let grads = layers.iter().map(|l| Tensor::new(0, l.out_dim())).collect();
        Mlp {
            layers,
            acts,
            grads,
            scratch: Tensor::default(),
            eval_sink: Tensor::default(),
        }
    }

    /// `n_in -> hidden[0] -> ... -> 1`, with leaky ReLU hidden layers.
    pub fn feedforward(n_in: usize, hidden: &[usize], rng: &mut StdRng) -> Self {
        let mut layers: Vec<Box<dyn Layer>> = Vec::new();
        let mut prev = n_in;
        for &h in hidden {
            layers.push(Box::new(Linear::new(prev, h, rng)));
            layers.push(Box::new(Relu::leaky(h, LEAKY_SLOPE)));
            prev = h;
        }
        layers.push(Box::new(Linear::new(prev, 1, rng)));
        Mlp::new(layers)
    }

    pub fn forward(&mut self, x: &Tensor) -> &Tensor {
        for k in 0..self.layers.len() {
            if k == 0 {
                self.layers[0].forward(x, &mut self.acts[0]);
            } else {
                let (prev, cur) = self.acts.split_at_mut(k);
                self.layers[k].forward(&prev[k - 1], &mut cur[0]);
            }
        }
        self.acts.last().unwrap()
    }

    /// Expects dL/d(output) already written into `grads.last()`.
    pub fn backward(&mut self, x: &Tensor) {
        for k in (0..self.layers.len()).rev() {
            let xk: &Tensor = if k == 0 { x } else { &self.acts[k - 1] };
            let yk = &self.acts[k];
            if k == 0 {
                self.layers[0].backward(xk, yk, &self.grads[0], &mut self.scratch);
            } else {
                let (lo, hi) = self.grads.split_at_mut(k);
                self.layers[k].backward(xk, yk, &hi[0], &mut lo[k - 1]);
            }
        }
    }

    pub fn params_and_grads(&mut self) -> ParamSlots<'_> {
        self.layers
            .iter_mut()
            .flat_map(|l| l.params_and_grads())
            .collect()
    }

    /// One forward/backward/update over a batch. Returns the batch loss.
    pub fn train_step(&mut self, x: &Tensor, y: &[f32], w: &[f32], opt: &mut Adam) -> f32 {
        self.forward(x);
        let last = self.acts.len() - 1;
        let l = loss_and_grad(&self.acts[last], y, w, &mut self.grads[last]);
        self.backward(x);
        opt.step(self.params_and_grads());
        l
    }

    /// Mean loss over a set WITHOUT touching parameters, gradients, or the RNG.
    ///
    /// Reuses [`loss_and_grad`] rather than re-deriving the loss, so the
    /// reported number cannot drift from the one being optimized. The gradient
    /// it writes goes to a throwaway buffer; the caller's `grads` are untouched,
    /// and `train_step` overwrites `grads.last()` before every backward pass
    /// regardless.
    ///
    /// The sink is the persistent [`Self::eval_sink`] rather than a local, so an
    /// every-epoch early-stopping evaluation allocates nothing in steady state —
    /// the claim the module header makes. Reusing it cannot change the returned
    /// number: `loss_and_grad` reshapes the sink and writes every row of it
    /// before returning, and the loss it returns never reads the gradient.
    pub fn eval_loss(&mut self, x: &Tensor, y: &[f32], w: &[f32]) -> f32 {
        self.forward(x);
        let last = self.acts.len() - 1;
        loss_and_grad(&self.acts[last], y, w, &mut self.eval_sink)
    }

    /// Copy every parameter into `out`, in [`Self::params_and_grads`] order.
    ///
    /// `out` is reused across epochs, so a whole training run snapshots into one
    /// allocation. The order is a pure function of the layer stack, which never
    /// changes after construction, so [`Self::load_params_from`] can walk it
    /// back with nothing but the length of each slot.
    fn snapshot_params_into(&mut self, out: &mut Vec<f32>) {
        out.clear();
        for (p, _, _) in self.params_and_grads() {
            out.extend_from_slice(p);
        }
    }

    /// Inverse of [`Self::snapshot_params_into`].
    fn load_params_from(&mut self, src: &[f32]) {
        let mut at = 0;
        for (p, _, _) in self.params_and_grads() {
            let n = p.len();
            p.copy_from_slice(&src[at..at + n]);
            at += n;
        }
        debug_assert_eq!(
            at,
            src.len(),
            "snapshot length disagrees with the layer stack"
        );
    }

    /// Production training path over regenerated, batch-sized chunks. Candidate
    /// order was randomized before fold assignment, so each chunk already has
    /// random membership; epochs shuffle chunk IDs only. `visit_train` supplies
    /// each transformed chunk to the consumer and retains none of them.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(crate) fn train_reporting_streaming(
        &mut self,
        cfg: &MlpConfig,
        n_train_batches: usize,
        visit_train: &mut dyn FnMut(&[usize], &mut dyn FnMut(&MlpBatch)),
        mut visit_val: Option<&mut dyn FnMut(&mut dyn FnMut(&MlpBatch))>,
        opt: &mut Adam,
        rng: &mut StdRng,
        tag: &str,
        epoch_finished: &(dyn Fn() + Sync),
    ) -> TrainOutcome {
        let patience = match (cfg.early_stopping_patience, visit_val.is_some()) {
            (Some(p), true) => Some(p),
            _ => None,
        };
        let mut order: Vec<usize> = (0..n_train_batches).collect();
        let mut last_epoch = 0.0f32;
        let mut best_val = f32::INFINITY;
        let mut best_epoch = None;
        let mut best_train = 0.0f32;
        let mut best_params = Vec::new();
        let mut stale = 0usize;
        let mut epochs_run = 0usize;

        for epoch in 0..cfg.epochs {
            for i in (1..order.len()).rev() {
                let j = rng.random_range(0..=i);
                order.swap(i, j);
            }

            let mut running = 0.0f32;
            let mut seen = 0usize;
            {
                let mut consume_train = |b: &MlpBatch| {
                    running += self.train_step(&b.x, &b.y, &b.w, opt);
                    seen += 1;
                };
                visit_train(&order, &mut consume_train);
            }
            last_epoch = if seen == 0 {
                0.0
            } else {
                running / seen as f32
            };
            epochs_run = epoch + 1;

            let traced = tracing::enabled!(tracing::Level::DEBUG);
            let val_loss = match (visit_val.as_mut(), patience.is_some() || traced) {
                (Some(visit), true) => {
                    let mut weighted = 0.0f64;
                    let mut total_weight = 0.0f64;
                    {
                        let mut consume_val = |b: &MlpBatch| {
                            let weight: f32 = b.w.iter().sum();
                            let loss = self.eval_loss(&b.x, &b.y, &b.w);
                            weighted += loss as f64 * weight as f64;
                            total_weight += weight as f64;
                        };
                        (**visit)(&mut consume_val);
                    }
                    Some(if total_weight > 0.0 {
                        (weighted / total_weight) as f32
                    } else {
                        0.0
                    })
                }
                _ => None,
            };

            if traced {
                match val_loss {
                    Some(vl) => tracing::debug!(
                        "{tag}epoch {}/{}: train_loss={:.6} held_out_loss={:.6}",
                        epoch + 1,
                        cfg.epochs,
                        last_epoch,
                        vl,
                    ),
                    None => tracing::debug!(
                        "{tag}epoch {}/{}: train_loss={:.6} (no held-out set)",
                        epoch + 1,
                        cfg.epochs,
                        last_epoch,
                    ),
                }
            }

            epoch_finished();
            let Some(p) = patience else { continue };
            let vl = val_loss.expect("early stopping requires a held-out measurement");
            if vl < best_val {
                best_val = vl;
                best_epoch = Some(epoch);
                best_train = last_epoch;
                self.snapshot_params_into(&mut best_params);
                stale = 0;
            } else {
                stale += 1;
                if stale >= p {
                    break;
                }
            }
        }

        let restored = match best_epoch {
            Some(b) if patience.is_some() && b + 1 != epochs_run => {
                self.load_params_from(&best_params);
                true
            }
            _ => false,
        };
        if restored {
            tracing::debug!(
                "{tag}early stop: ran {} of {} epochs, restored epoch {} (held_out_loss={:.6})",
                epochs_run,
                cfg.epochs,
                best_epoch.unwrap() + 1,
                best_val,
            );
        }

        TrainOutcome {
            final_train_loss: if restored { best_train } else { last_epoch },
            epochs_run,
            best_epoch,
            best_val_loss: best_epoch.map(|_| best_val),
            restored,
        }
    }

    /// `sum_o |W1[j, o]|` over the first weight-bearing layer, per input
    /// column. The closest analogue to LDA's `|coef|` that a net admits —
    /// there is no gain equivalent.
    ///
    /// Indexed by TRANSFORMED input, so callers map back through
    /// [`ColumnTransform::lane_of_input`] to get a lane-indexed vector.
    pub fn input_importance(&self) -> Vec<f32> {
        self.layers
            .iter()
            .find_map(|l| l.abs_in_weights())
            .unwrap_or_default()
    }
}

/// One reusable transformed batch buffer. Production fitting owns exactly two:
/// the producer fills one from candidate rows while forward/backward consumes
/// the other, then the buffers trade places.
pub struct MlpBatch {
    pub x: Tensor,
    pub y: Vec<f32>,
    pub w: Vec<f32>,
}

impl MlpBatch {
    /// Reusable producer/consumer buffer with capacity for one full batch and
    /// initially zero live rows.
    pub fn buffer(row_capacity: usize, cols: usize) -> Self {
        let mut x = Tensor::new(row_capacity, cols);
        x.reshape(0, cols);
        Self {
            x,
            y: Vec::with_capacity(row_capacity),
            w: Vec::with_capacity(row_capacity),
        }
    }
}

/// Result of one streaming fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrainOutcome {
    /// Mean training loss of the epoch whose weights the model now holds — the
    /// last epoch run, or the RESTORED epoch when early stopping rolled back.
    /// Reporting the last epoch's loss after a rollback would describe weights
    /// that were thrown away.
    pub final_train_loss: f32,
    /// Epochs actually executed, `<= cfg.epochs`.
    pub epochs_run: usize,
    /// Zero-based epoch with the best held-out loss, or `None` when there was
    /// no held-out set (or every measurement was `NaN`).
    pub best_epoch: Option<usize>,
    /// The held-out loss at `best_epoch`.
    pub best_val_loss: Option<f32>,
    /// Whether the parameters were rolled back to `best_epoch`.
    pub restored: bool,
}

// ---------------------------------------------------------------- input transform

/// What the column stats pass decided about one lane column.
#[derive(Debug, Clone, Copy, PartialEq)]
enum ColKind {
    /// Never non-finite on the training rows: standardize, no companion.
    Clean,
    /// Sometimes non-finite: standardize, impute to the train-fold mean (i.e.
    /// 0 post-standardization), and emit an `_isna` companion.
    Missable,
}

struct ColSpec {
    lane: usize,
    mean: f64,
    inv_std: f64,
    kind: ColKind,
}

/// Per-column cull / standardize / impute, fitted on training rows only.
///
/// Three jobs, one pass, because they all read the same per-column moments:
///
///  1. **Cull** columns with no finite value at all, and columns whose measured
///     std is `<= MIN_STD`. What that GUARANTEES is narrower than "constant
///     columns are dropped" — see below.
///  2. **Standardize** the survivors, which is what lets the default
///     architecture skip BatchNorm entirely.
///  3. **Impute** non-finite values to the column mean and emit an `_isna`
///     companion, so missingness survives as signal instead of being laundered
///     into the mean.
///
/// Every statistic here is fitted on the rows it is handed and then applied
/// unchanged to held-out rows. Fitting across all rows would use held-out
/// feature values.
///
/// # What the cull does NOT guarantee
/// The variance is the textbook-unstable `sumsq/n - mean^2` form, so a
/// bit-for-bit CONSTANT column of realistic magnitude keeps enough
/// floating-point residue to clear `MIN_STD` and SURVIVES. Measured on the
/// suite's all-constant fixture: 8 and 16 train rows cull all 101 columns, while
/// 30 leaves 10 of them alive. So the cull reliably drops all-non-finite
/// columns, and drops constant ones only sometimes — as a function of the row
/// count, not of the data.
///
/// This is numerically benign, which is why it is documented rather than fixed.
/// The `.max(0.0)` closes the only real hazard (a negative variance would give a
/// `NAN` std, hence `NAN` inputs), and a surviving constant column standardizes
/// to values bounded by roughly `sqrt(f64::EPSILON)`, about `1.5e-8` — noise the
/// net cannot learn from. What it costs is honesty in the report: such a column
/// silently occupies an MLP input and a row in the importance sidecar, reading as
/// "uninformative" rather than as "constant". Switching to Welford or a two-pass
/// variance would make the cull dependable; that is a numerics change to a
/// shipping path, not a doc fix, and has not been made.
pub struct ColumnTransform {
    ncols_lane: usize,
    cols: Vec<ColSpec>,
    /// Lane indices dropped by the cull, in lane order.
    culled: Vec<usize>,
    /// Number of `_isna` companions, appended after the standardized block.
    n_isna: usize,
}

/// Below this, a column counts as constant and gets culled.
const MIN_STD: f64 = 1e-12;

impl ColumnTransform {
    /// Fit over `rows` of the row-major lane matrix `feat` (`feat[i*ncols + j]`).
    #[cfg(test)]
    pub fn fit(feat: &[f64], ncols_lane: usize, rows: &[usize]) -> Self {
        Self::fit_streaming(ncols_lane, rows.iter().copied(), |i, out| {
            out.copy_from_slice(&feat[i * ncols_lane..(i + 1) * ncols_lane]);
        })
    }

    /// Fit from rows supplied into one reusable scratch buffer. This is
    /// numerically the same one-pass moment calculation as [`Self::fit`], but it
    /// does not require a retained raw lane matrix.
    pub fn fit_streaming(
        ncols_lane: usize,
        rows: impl IntoIterator<Item = usize>,
        mut write_row: impl FnMut(usize, &mut [f64]),
    ) -> Self {
        let mut sum = vec![0.0f64; ncols_lane];
        let mut sumsq = vec![0.0f64; ncols_lane];
        let mut finite = vec![0u64; ncols_lane];
        let mut nonfinite = vec![0u64; ncols_lane];
        let mut row = vec![0.0f64; ncols_lane];

        for i in rows {
            write_row(i, &mut row);
            for (j, &v) in row.iter().enumerate() {
                if v.is_finite() {
                    sum[j] += v;
                    sumsq[j] += v * v;
                    finite[j] += 1;
                } else {
                    nonfinite[j] += 1;
                }
            }
        }

        let mut cols = Vec::with_capacity(ncols_lane);
        let mut culled = Vec::new();
        let mut n_isna = 0;
        for j in 0..ncols_lane {
            if finite[j] == 0 {
                culled.push(j);
                continue;
            }
            let n = finite[j] as f64;
            let mean = sum[j] / n;
            let var = (sumsq[j] / n - mean * mean).max(0.0);
            let std = var.sqrt();
            if std <= MIN_STD {
                culled.push(j);
                continue;
            }
            let kind = if nonfinite[j] > 0 {
                n_isna += 1;
                ColKind::Missable
            } else {
                ColKind::Clean
            };
            cols.push(ColSpec {
                lane: j,
                mean,
                inv_std: 1.0 / std,
                kind,
            });
        }

        ColumnTransform {
            ncols_lane,
            cols,
            culled,
            n_isna,
        }
    }

    /// Width of the transformed row: surviving columns, then `_isna` companions.
    pub fn width(&self) -> usize {
        self.cols.len() + self.n_isna
    }

    /// Lane indices dropped by the cull. Callers report these — a culled
    /// column's `0.0` importance is otherwise indistinguishable from an
    /// uninformative one's.
    pub fn culled(&self) -> &[usize] {
        &self.culled
    }

    /// Lane column count this transform was fitted against.
    pub fn ncols_lane(&self) -> usize {
        self.ncols_lane
    }

    /// Non-finite values seen on a column the fit saw as clean. Always zero on
    /// the fitting rows by construction; nonzero here means a scored row
    /// carried a NaN into a column that never had one during training.
    ///
    /// Deliberately a returned count rather than a `debug_assert!`: assertions
    /// compile out in release, which is exactly where a production run would
    /// hit this.
    #[cfg(test)]
    pub fn check_clean(&self, feat: &[f64], rows: &[usize]) -> Vec<usize> {
        let mut bad = Vec::new();
        for spec in self.cols.iter().filter(|c| c.kind == ColKind::Clean) {
            let hit = rows
                .iter()
                .any(|&i| !feat[i * self.ncols_lane + spec.lane].is_finite());
            if hit {
                bad.push(spec.lane);
            }
        }
        bad
    }

    /// Mark clean-at-fit lane columns that are non-finite in this one row.
    /// Streaming inference accumulates these marks while transforming rows,
    /// avoiding the old gather + diagnostic pass over a raw matrix.
    pub fn mark_dirty_clean(&self, row: &[f64], dirty: &mut [bool]) {
        debug_assert_eq!(row.len(), self.ncols_lane);
        debug_assert_eq!(dirty.len(), self.ncols_lane);
        for spec in self.cols.iter().filter(|c| c.kind == ColKind::Clean) {
            if !row[spec.lane].is_finite() {
                dirty[spec.lane] = true;
            }
        }
    }

    /// Transform one lane row into `out` (`out.len() == self.width()`).
    pub fn apply(&self, row: &[f64], out: &mut [f32]) {
        debug_assert_eq!(row.len(), self.ncols_lane);
        debug_assert_eq!(out.len(), self.width());
        let mut isna_at = self.cols.len();
        for (k, spec) in self.cols.iter().enumerate() {
            let v = row[spec.lane];
            if v.is_finite() {
                out[k] = ((v - spec.mean) * spec.inv_std) as f32;
                if spec.kind == ColKind::Missable {
                    out[isna_at] = 0.0;
                    isna_at += 1;
                }
            } else {
                // Impute to the mean, i.e. exactly 0 in standardized space.
                out[k] = 0.0;
                if spec.kind == ColKind::Missable {
                    out[isna_at] = 1.0;
                    isna_at += 1;
                }
            }
        }
    }

    /// Map a transformed-input index back to its lane column, if it is a
    /// standardized feature rather than an `_isna` companion.
    pub fn lane_of_input(&self, k: usize) -> Option<usize> {
        self.cols.get(k).map(|c| c.lane)
    }

    /// The other half of [`Self::lane_of_input`]: the lane column an `_isna`
    /// COMPANION belongs to, or `None` when `k` is a standardized feature (or
    /// past the end).
    ///
    /// Companions have no lane of their own — they are appended after the
    /// standardized block, in the lane order of the missable columns they
    /// flag — so a caller that wants a lane-indexed vector has no way to place
    /// them without this. It is a lookup rather than a public field so the
    /// companion LAYOUT stays an internal detail of [`Self::apply`]; the two
    /// must agree, and this is the only other place that encodes the order.
    pub fn isna_lane_of_input(&self, k: usize) -> Option<usize> {
        let nth = k.checked_sub(self.cols.len())?;
        self.cols
            .iter()
            .filter(|c| c.kind == ColKind::Missable)
            .nth(nth)
            .map(|c| c.lane)
    }
}

// ---------------------------------------------------------------- config

/// Fixed rescoring hyperparameters. The input width is derived at fit time.
#[derive(Debug, Clone, PartialEq)]
pub struct MlpConfig {
    pub hidden: Vec<usize>,
    pub lr: f32,
    pub weight_decay: f32,
    pub epochs: usize,
    pub batch_size: usize,
    pub seed: u64,
    /// Stop after this many non-improving validation epochs and restore the best
    /// weights. `None` runs the full epoch budget.
    pub early_stopping_patience: Option<usize>,
}

impl Default for MlpConfig {
    fn default() -> Self {
        Self::compiled_default()
    }
}

impl MlpConfig {
    pub fn compiled_default() -> Self {
        MlpConfig {
            hidden: vec![32, 16],
            lr: 1.2e-3,
            weight_decay: 1e-4,
            epochs: 60,
            batch_size: 1024,
            seed: 0x2545_F491_4F6C_DD1D,
            early_stopping_patience: Some(3),
        }
    }

    /// Per-fold RNG: the configured seed mixed with the fold index, so folds
    /// differ but a rerun of the same fold does not.
    pub fn rng_for_fold(&self, fold: usize) -> StdRng {
        StdRng::seed_from_u64(self.seed ^ (fold as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }
}

// ---------------------------------------------------------------- rng helper

/// Box-Muller, kept for He init.
fn normal(rng: &mut StdRng) -> f32 {
    let u1: f32 = rng.random::<f32>().max(1e-7);
    let u2: f32 = rng.random::<f32>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f32::consts::PI * u2).cos()
}

#[cfg(test)]
mod test {
    use super::*;

    fn tensor(rows: usize, cols: usize, values: &[f32]) -> Tensor {
        let mut out = Tensor::new(rows, cols);
        out.live_mut().copy_from_slice(values);
        out
    }

    #[test]
    fn gemm_handles_forward_and_transposed_layouts() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let b = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let mut c = vec![0.0; 4];
        gemm_f32(2, 3, 2, &a, 3, 1, &b, 2, 1, 0.0, &mut c);
        assert_eq!(c, vec![58.0, 64.0, 139.0, 154.0]);

        let mut ata = vec![0.0; 9];
        gemm_f32(3, 2, 3, &a, 1, 3, &a, 3, 1, 0.0, &mut ata);
        assert_eq!(
            ata,
            vec![17.0, 22.0, 27.0, 22.0, 29.0, 36.0, 27.0, 36.0, 45.0]
        );
    }

    #[test]
    fn linear_input_gradient_matches_finite_differences() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut layer = Linear::new(3, 2, &mut rng);
        let x = tensor(2, 3, &[0.2, -0.7, 1.1, -0.4, 0.3, 0.9]);
        let gy = tensor(2, 2, &[0.5, -0.2, 0.7, 0.1]);
        let mut y = Tensor::default();
        layer.forward(&x, &mut y);
        let mut gx = Tensor::default();
        layer.backward(&x, &y, &gy, &mut gx);

        let eps = 1e-3;
        for i in 0..x.live().len() {
            let mut plus = x.clone();
            let mut minus = x.clone();
            plus.live_mut()[i] += eps;
            minus.live_mut()[i] -= eps;
            let mut yp = Tensor::default();
            let mut ym = Tensor::default();
            layer.forward(&plus, &mut yp);
            layer.forward(&minus, &mut ym);
            let lp: f32 = yp.live().iter().zip(gy.live()).map(|(a, b)| a * b).sum();
            let lm: f32 = ym.live().iter().zip(gy.live()).map(|(a, b)| a * b).sum();
            let numerical = (lp - lm) / (2.0 * eps);
            assert!((gx.live()[i] - numerical).abs() < 2e-3);
        }
    }

    #[test]
    fn leaky_relu_preserves_a_negative_gradient() {
        let x = tensor(1, 2, &[-2.0, 3.0]);
        let gy = tensor(1, 2, &[4.0, 5.0]);
        let mut relu = Relu::leaky(2, LEAKY_SLOPE);
        let mut y = Tensor::default();
        let mut gx = Tensor::default();
        relu.forward(&x, &mut y);
        relu.backward(&x, &y, &gy, &mut gx);
        assert_eq!(y.live(), &[-2.0 * LEAKY_SLOPE, 3.0]);
        assert_eq!(gx.live(), &[4.0 * LEAKY_SLOPE, 5.0]);
    }

    #[test]
    fn bce_gradient_matches_finite_differences_and_stays_finite() {
        let logits = tensor(4, 1, &[-20.0, -0.7, 0.9, 20.0]);
        let y = [0.0, 1.0, 0.0, 1.0];
        let w = [1.0, 0.5, 1.0, 0.5];
        let mut grad = Tensor::default();
        let loss = loss_and_grad(&logits, &y, &w, &mut grad);
        assert!(loss.is_finite());
        assert!(grad.live().iter().all(|v| v.is_finite()));

        let eps = 1e-3;
        for i in 0..4 {
            let mut plus = logits.clone();
            let mut minus = logits.clone();
            plus.live_mut()[i] += eps;
            minus.live_mut()[i] -= eps;
            let mut sink = Tensor::default();
            let lp = loss_and_grad(&plus, &y, &w, &mut sink);
            let lm = loss_and_grad(&minus, &y, &w, &mut sink);
            assert!((grad.live()[i] - (lp - lm) / (2.0 * eps)).abs() < 2e-3);
        }
    }

    #[test]
    fn adam_updates_parameters_and_clears_gradients() {
        let mut p = vec![1.0, -1.0];
        let mut g = vec![0.5, -0.25];
        let mut adam = Adam::new(1e-2).with_weight_decay(1e-4);
        adam.step(vec![(p.as_mut_slice(), g.as_mut_slice(), true)]);
        assert_ne!(p, vec![1.0, -1.0]);
        assert_eq!(g, vec![0.0, 0.0]);
    }

    #[test]
    fn transform_culls_dead_columns_and_tracks_missingness() {
        let feat = vec![1.0, 5.0, f64::NAN, 2.0, 5.0, 8.0, 3.0, 5.0, 10.0];
        let tx = ColumnTransform::fit(&feat, 3, &[0, 1, 2]);
        assert_eq!(tx.culled(), &[1]);
        assert_eq!(tx.width(), 3);
        assert_eq!(tx.lane_of_input(0), Some(0));
        assert_eq!(tx.lane_of_input(1), Some(2));
        assert_eq!(tx.isna_lane_of_input(2), Some(2));

        let mut out = vec![0.0; tx.width()];
        tx.apply(&feat[..3], &mut out);
        assert_eq!(out[1], 0.0);
        assert_eq!(out[2], 1.0);
    }

    #[test]
    fn streaming_transform_matches_materialized_fit() {
        let feat = vec![1.0, 2.0, 3.0, 2.0, f64::NAN, 4.0, 4.0, 8.0, 9.0];
        let rows = [0, 1, 2];
        let expected = ColumnTransform::fit(&feat, 3, &rows);
        let actual = ColumnTransform::fit_streaming(3, rows, |i, out| {
            out.copy_from_slice(&feat[i * 3..i * 3 + 3]);
        });
        assert_eq!(actual.width(), expected.width());
        assert_eq!(actual.culled(), expected.culled());
        for row in feat.chunks_exact(3) {
            let mut a = vec![0.0; actual.width()];
            let mut b = vec![0.0; expected.width()];
            actual.apply(row, &mut a);
            expected.apply(row, &mut b);
            assert_eq!(a, b);
        }
    }

    #[test]
    fn transform_statistics_are_train_only() {
        let feat = vec![0.0, 2.0, 1000.0];
        let tx = ColumnTransform::fit(&feat, 1, &[0, 1]);
        let mut out = [0.0];
        tx.apply(&feat[2..], &mut out);
        assert!(out[0] > 100.0);
    }

    #[test]
    fn feedforward_is_seed_reproducible() {
        let x = tensor(2, 3, &[1.0, 2.0, 3.0, -1.0, 0.5, 2.0]);
        let mut r1 = StdRng::seed_from_u64(9);
        let mut r2 = StdRng::seed_from_u64(9);
        let mut a = Mlp::feedforward(3, &[4, 2], &mut r1);
        let mut b = Mlp::feedforward(3, &[4, 2], &mut r2);
        assert_eq!(a.forward(&x).live(), b.forward(&x).live());
    }

    #[test]
    fn streaming_training_learns_a_separable_batch() {
        let batch = MlpBatch {
            x: tensor(4, 1, &[-2.0, -1.0, 1.0, 2.0]),
            y: vec![0.0, 0.0, 1.0, 1.0],
            w: vec![1.0; 4],
        };
        let cfg = MlpConfig {
            hidden: vec![4],
            epochs: 80,
            batch_size: 4,
            early_stopping_patience: None,
            ..MlpConfig::default()
        };
        let mut rng = cfg.rng_for_fold(0);
        let mut net = Mlp::feedforward(1, &cfg.hidden, &mut rng);
        let before = net.eval_loss(&batch.x, &batch.y, &batch.w);
        let mut opt = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
        let mut visit = |order: &[usize], consume: &mut dyn FnMut(&MlpBatch)| {
            for _ in order {
                consume(&batch);
            }
        };
        net.train_reporting_streaming(&cfg, 1, &mut visit, None, &mut opt, &mut rng, "", &|| {});
        let after = net.eval_loss(&batch.x, &batch.y, &batch.w);
        assert!(after < before, "{before} -> {after}");
    }

    #[test]
    fn fold_rng_is_stable_and_fold_specific() {
        let cfg = MlpConfig::default();
        let mut a = cfg.rng_for_fold(2);
        let mut b = cfg.rng_for_fold(2);
        let mut c = cfg.rng_for_fold(3);
        assert_eq!(a.random::<u64>(), b.random::<u64>());
        assert_ne!(a.random::<u64>(), c.random::<u64>());
    }
}
