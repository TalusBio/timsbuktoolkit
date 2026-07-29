//! Dependency-free MLP for rescoring, ported from `dev/mlp.rs`.
//!
//! Differences from the reference implementation, all deliberate:
//!
//!  * **Single-logit head with BCE-with-logits**, not a 2-class softmax. The
//!    score handed to `assign_qval` only needs a monotone ranking, so the
//!    probability transform is dead weight.
//!  * **Per-row sample weights**, so the MLP inherits the same 0.5/1.0
//!    target/decoy class balance the GBM uses (`cv.rs`).
//!  * **No BatchNorm in the default architecture.** Inputs are standardized once
//!    by [`ColumnTransform`], which is where the leading BatchNorm's job went.
//!    `BatchNorm` is still ported and gradient-tested — it is the first thing to
//!    try if training turns out unstable — but keeping it off the default path
//!    means a non-finite value cannot poison a batch mean and corrupt every
//!    other row in the batch.
//!  * **Seeded `StdRng`** rather than the vendored xorshift, so the crate has
//!    one RNG story.
//!  * **No BatchNorm folding / `FoldedMlp` / `Sgd` / softmax / argmax**.
//!    Inference here is batched and in-process, so folding buys nothing.
//!
//! All tensors are flat `Vec<f32>`, row-major, `[batch][features]`. The `Mlp`
//! owns every activation and gradient buffer; buffers grow but never shrink, so
//! steady-state training does no allocation. The public boundary is `f64`,
//! matching the lane matrices the rest of `ml::` passes around.

// Index-based loops are deliberate throughout. These are numeric kernels where
// one index walks several parallel arrays at once (`w`/`gw`, `mean`/`var`/
// `inv_std`) and where the inner AXPY loops are shaped for vectorization.
// Rewriting them as iterator zips needs 3- and 4-way nests and sinks the
// codegen — the same tradeoff `ml/lda.rs` documents at its own hot loops.
#![allow(clippy::needless_range_loop)]

use rand::rngs::StdRng;
use rand::{
    Rng as _,
    SeedableRng,
};

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

pub trait Layer {
    fn out_dim(&self) -> usize;
    fn forward(&mut self, x: &Tensor, out: &mut Tensor, training: bool);
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
    /// learning rate, epoch count, or minibatch noise recovers from.
    ///
    /// Observed, not theorized: with zero biases a 2-16-8-1 net on XOR settles
    /// at exactly 0.75 accuracy and `ln(2)/2` loss for some seeds and never
    /// leaves. `dev/mlp.rs` never hit it because it puts a `BatchNorm` at layer
    /// 0, whose learnable `beta` supplies the shift zero biases do not; the
    /// default architecture here has no BatchNorm, so the bias has to.
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

    fn forward(&mut self, x: &Tensor, out: &mut Tensor, _training: bool) {
        let (n, ni, no) = (x.rows, self.n_in, self.n_out);
        out.reshape(n, no);
        for b in 0..n {
            let xr = x.row(b);
            let or = out.row_mut(b);
            or.copy_from_slice(&self.b);
            // AXPY order: inner loop is contiguous and vectorizes.
            for i in 0..ni {
                let xi = xr[i];
                let wr = &self.w[i * no..i * no + no];
                for (v, w) in or.iter_mut().zip(wr) {
                    *v += xi * w;
                }
            }
        }
    }

    fn backward(&mut self, x: &Tensor, _y: &Tensor, gy: &Tensor, gx: &mut Tensor) {
        let (n, ni, no) = (gy.rows, self.n_in, self.n_out);
        gx.reshape(n, ni);
        for b in 0..n {
            let gr = gy.row(b);
            let xr = x.row(b);
            let gxr = gx.row_mut(b);
            for o in 0..no {
                self.gb[o] += gr[o];
            }
            for i in 0..ni {
                // Two separate loops: the AXPY vectorizes, the reduction
                // doesn't, and fusing them would sink both.
                let xi = xr[i];
                let gwr = &mut self.gw[i * no..i * no + no];
                for (g, d) in gwr.iter_mut().zip(gr) {
                    *g += xi * d;
                }
                let wr = &self.w[i * no..i * no + no];
                let mut acc = 0.0;
                for (d, w) in gr.iter().zip(wr) {
                    acc += d * w;
                }
                gxr[i] = acc;
            }
        }
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
///
/// Measured on the XOR test in this module: with pure ReLU, seeds 7 and 13 both
/// pinned at exactly 0.75 accuracy / `ln(2)/2` loss (half the samples emitting
/// logit 0) and never moved across 400 epochs. Non-zero bias init made it
/// rarer — seed 7 recovered — but did not remove it, because it changes where
/// units start, not what happens once one dies.
pub struct Relu {
    dim: usize,
    slope: f32,
}

impl Relu {
    /// Pure ReLU, zero negative slope. Retained for gradient tests and for
    /// anyone who wants the exact classical behavior; see the type docs for
    /// why it is not the default.
    pub fn new(dim: usize) -> Self {
        Relu { dim, slope: 0.0 }
    }

    pub fn leaky(dim: usize, slope: f32) -> Self {
        Relu { dim, slope }
    }
}

impl Layer for Relu {
    fn out_dim(&self) -> usize {
        self.dim
    }

    fn forward(&mut self, x: &Tensor, out: &mut Tensor, _training: bool) {
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

// ---------------------------------------------------------------- square

/// `f(x) = x^2`, as an alternative to the rectifier.
///
/// Reported to work well by the project owner in cases where DEAD UNITS were the
/// problem, which is the case for it: `x^2` has zero gradient at exactly one
/// point rather than over a whole half-line, so there is no absorbing dead state
/// of the kind [`Relu`] documents. Nonlinear on its own — a single squaring layer
/// already separates XOR — and even, so it responds to the MAGNITUDE of a
/// standardized feature rather than its sign, which is a different inductive bias
/// from a rectifier, not merely a smoother one.
///
/// # It is unbounded and even, and that is a real hazard
/// A rectifier is 1-Lipschitz and cannot amplify its input; squaring can, and
/// does. Inputs here are standardized, so |x| = 3 is an ordinary value, and one
/// layer turns that into 9. Two stacked squaring layers compound it (9 -> 81
/// before the second layer's weights are even applied), so the forward pass can
/// reach a non-finite logit from a weight update that a rectifier would merely
/// have made large. That is a live failure mode for this activation in a way it
/// is not for [`Relu`].
///
/// The failure is DIAGNOSED rather than silent: `MlpFoldError::NonFiniteScore`
/// checks every scored logit and aborts naming the row, so a diverged net does
/// not reach `assign_qval`. But diagnosed is not prevented — expect to need a
/// SMALLER `lr`, or a larger `weight_decay`, than the rectifier defaults in
/// [`MlpConfig::compiled_default`], which were swept with a rectifier in place
/// and carry no guarantee here.
///
/// # What was actually measured, and what was not
/// The 2-input XOR fixture in this module's tests does NOT diverge: it reaches
/// the same accuracy as the rectifier at the shipped `lr` of 1.2e-3 on all eight
/// seeds, at a lower training loss on every one, and it still does not go
/// non-finite at 400x that step size. Read that as "the layer and its gradient
/// are right", not as reassurance — 2 inputs of magnitude ~1 is the case where
/// squaring cannot amplify much.
///
/// The 24-column overfitting fixture is the closer analogue and it shows the
/// amplification plainly. At the same `lr`, `hidden` and epoch budget, the
/// squaring net reaches a LOWER train loss and a held-out loss of 11-21 against
/// the rectifier's 1.1-1.6 (three seeds), i.e. it emits far larger logits on rows
/// it has not seen; at `lr = 1e-1` it drives train loss to exactly 0 with a
/// held-out loss of 4e6-4.5e7. Nothing there is non-finite, so the guard rail
/// never trips — the realistic hazard is not a crash but a net whose held-out
/// logits are enormous, and the 101- and 128-column rescore lanes have not been
/// measured at all.
pub struct Square {
    dim: usize,
}

impl Square {
    pub fn new(dim: usize) -> Self {
        Square { dim }
    }
}

impl Layer for Square {
    fn out_dim(&self) -> usize {
        self.dim
    }

    fn forward(&mut self, x: &Tensor, out: &mut Tensor, _training: bool) {
        out.reshape(x.rows, x.cols);
        let (src, dst) = (x.live(), out.live_mut());
        for k in 0..src.len() {
            dst[k] = src[k] * src[k];
        }
    }

    /// `dL/dx = gy * 2x`, and it READS `x` — unlike [`Relu::backward`], which
    /// deliberately ignores `x` because `y > 0` recovers its mask. Squaring
    /// destroys the sign, so `y` alone cannot say whether the local slope is
    /// `+2|x|` or `-2|x|`; `sqrt(y)` would recover the magnitude and lose the
    /// sign, which is exactly the wrong half. The trait hands both in, so this
    /// costs nothing but is worth naming: a future refactor that drops the `x`
    /// argument because "`Relu` does not use it" would silently break this layer.
    fn backward(&mut self, x: &Tensor, _y: &Tensor, gy: &Tensor, gx: &mut Tensor) {
        gx.reshape(gy.rows, gy.cols);
        let (xv, gv, dst) = (x.live(), gy.live(), gx.live_mut());
        for k in 0..dst.len() {
            dst[k] = gv[k] * 2.0 * xv[k];
        }
    }
}

/// Which activation [`Mlp::feedforward`] puts between hidden layers.
///
/// Default is [`Activation::LeakyRelu`] at [`LEAKY_SLOPE`], which is what every
/// measurement in [`MlpConfig::compiled_default`] was taken with.
/// [`Activation::Square`] is EXPERIMENTAL — see [`Square`] for the divergence
/// hazard and the learning-rate caveat that comes with it.
///
/// Neither variant draws from the RNG, so switching between them does not move
/// the initialization stream: a `LeakyRelu` fit is bit-identical to what it was
/// before this choice existed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Activation {
    #[default]
    LeakyRelu,
    Square,
}

impl Activation {
    fn layer(self, dim: usize) -> Box<dyn Layer> {
        match self {
            Activation::LeakyRelu => Box::new(Relu::leaky(dim, LEAKY_SLOPE)),
            Activation::Square => Box::new(Square::new(dim)),
        }
    }
}

// ---------------------------------------------------------------- batch norm

/// Ported and gradient-tested, but NOT part of the default architecture — see
/// the module docs. Kept as the first escalation if standardized inputs turn
/// out not to be enough to train stably.
pub struct BatchNorm {
    dim: usize,
    gamma: Vec<f32>,
    beta: Vec<f32>,
    running_mean: Vec<f32>,
    running_var: Vec<f32>,
    momentum: f32,
    eps: f32,
    // caches / scratch, allocated once
    xhat: Tensor,
    inv_std: Vec<f32>,
    mean: Vec<f32>,
    var: Vec<f32>,
    s1: Vec<f32>,
    s2: Vec<f32>,
    g_gamma: Vec<f32>,
    g_beta: Vec<f32>,
}

impl BatchNorm {
    pub fn new(dim: usize) -> Self {
        BatchNorm {
            dim,
            gamma: vec![1.0; dim],
            beta: vec![0.0; dim],
            running_mean: vec![0.0; dim],
            running_var: vec![1.0; dim],
            momentum: 0.1,
            eps: 1e-5,
            xhat: Tensor::new(0, dim),
            inv_std: vec![0.0; dim],
            mean: vec![0.0; dim],
            var: vec![0.0; dim],
            s1: vec![0.0; dim],
            s2: vec![0.0; dim],
            g_gamma: vec![0.0; dim],
            g_beta: vec![0.0; dim],
        }
    }
}

impl Layer for BatchNorm {
    fn out_dim(&self) -> usize {
        self.dim
    }

    fn forward(&mut self, x: &Tensor, out: &mut Tensor, training: bool) {
        let (n, d) = (x.rows, self.dim);
        out.reshape(n, d);

        if !training {
            for b in 0..n {
                let (xr, or) = (x.row(b), out.row_mut(b));
                for j in 0..d {
                    let inv = 1.0 / (self.running_var[j] + self.eps).sqrt();
                    or[j] = self.gamma[j] * (xr[j] - self.running_mean[j]) * inv + self.beta[j];
                }
            }
            return;
        }

        let nf = n as f32;
        // Inner loops run over features so every access is contiguous.
        for v in self.mean.iter_mut() {
            *v = 0.0;
        }
        for b in 0..n {
            let xr = x.row(b);
            for j in 0..d {
                self.mean[j] += xr[j];
            }
        }
        for v in self.mean.iter_mut() {
            *v /= nf;
        }

        for v in self.var.iter_mut() {
            *v = 0.0;
        }
        for b in 0..n {
            let xr = x.row(b);
            for j in 0..d {
                let dv = xr[j] - self.mean[j];
                self.var[j] += dv * dv;
            }
        }
        for v in self.var.iter_mut() {
            *v /= nf;
        }

        for j in 0..d {
            self.inv_std[j] = 1.0 / (self.var[j] + self.eps).sqrt();
            self.running_mean[j] =
                (1.0 - self.momentum) * self.running_mean[j] + self.momentum * self.mean[j];
            self.running_var[j] =
                (1.0 - self.momentum) * self.running_var[j] + self.momentum * self.var[j];
        }

        self.xhat.reshape(n, d);
        for b in 0..n {
            let xr = x.row(b);
            let hr = self.xhat.row_mut(b);
            for j in 0..d {
                hr[j] = (xr[j] - self.mean[j]) * self.inv_std[j];
            }
            let hr = self.xhat.row(b);
            let or = out.row_mut(b);
            for j in 0..d {
                or[j] = self.gamma[j] * hr[j] + self.beta[j];
            }
        }
    }

    fn backward(&mut self, _x: &Tensor, _y: &Tensor, gy: &Tensor, gx: &mut Tensor) {
        let (n, d) = (gy.rows, self.dim);
        let nf = n as f32;
        gx.reshape(n, d);

        // Pass 1: per-feature reductions, feature-contiguous inner loop.
        for j in 0..d {
            self.s1[j] = 0.0;
            self.s2[j] = 0.0;
        }
        for b in 0..n {
            let (gr, hr) = (gy.row(b), self.xhat.row(b));
            for j in 0..d {
                let dxhat = gr[j] * self.gamma[j];
                self.s1[j] += dxhat;
                self.s2[j] += dxhat * hr[j];
                self.g_gamma[j] += gr[j] * hr[j];
                self.g_beta[j] += gr[j];
            }
        }

        // Pass 2: dx = inv_std/N * (N*dxhat - sum(dxhat) - xhat*sum(dxhat*xhat))
        for b in 0..n {
            let (gr, hr) = (gy.row(b), self.xhat.row(b));
            let gxr = gx.row_mut(b);
            for j in 0..d {
                let dxhat = gr[j] * self.gamma[j];
                gxr[j] = self.inv_std[j] / nf * (nf * dxhat - self.s1[j] - hr[j] * self.s2[j]);
            }
        }
    }

    fn params_and_grads(&mut self) -> ParamSlots<'_> {
        vec![
            (
                self.gamma.as_mut_slice(),
                self.g_gamma.as_mut_slice(),
                false,
            ),
            (self.beta.as_mut_slice(), self.g_beta.as_mut_slice(), false),
        ]
    }
}

// ---------------------------------------------------------------- loss

/// Loss shape for the single-logit head.
///
/// [`MlpLoss::Focal`] at `gamma = 0, alpha = 0.5` is exactly [`MlpLoss::Bce`],
/// so the focal form is a strict generalization and the default stays the plain
/// loss.
///
/// # A removed variant, and why
/// There used to be an `AsymFocal` here that weighted every row by
/// `p^gamma_class` — by the predicted TARGET probability regardless of label,
/// with a per-class exponent — on the argument that target labels are noisy
/// positives, so a confidently-wrong target should be DISCOUNTED rather than
/// emphasized. That argument is wrong for this pipeline: the rows it discounts
/// are the marginal ones, and the 1% FDR cut lives exactly among them. Removing
/// gradient signal from the decision boundary removes it from the only region
/// that determines the identification count. The replacement is the STANDARD
/// alpha-balanced focal loss, which up-weights hard examples on BOTH sides — the
/// opposite of what the removed variant did on the target side.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum MlpLoss {
    #[default]
    Bce,
    /// Standard alpha-balanced focal loss (Lin et al., 2017):
    ///
    /// ```text
    ///   FL(p_t) = -alpha_t * (1 - p_t)^gamma * log(p_t)
    ///   p_t     = p       if y == 1 else 1 - p
    ///   alpha_t = alpha   if y == 1 else 1 - alpha
    /// ```
    ///
    /// `1 - p_t` is the probability mass the model put on the WRONG class, so
    /// `(1 - p_t)^gamma` is near 1 for a row the model gets wrong and near 0 for
    /// one it gets easily right: `gamma` trades away the easy bulk to concentrate
    /// the gradient on the hard rows, in both classes. `gamma = 0` disables the
    /// modulation entirely.
    ///
    /// # `alpha` is normalized so `gamma = 0, alpha = 0.5` IS `Bce`
    /// `alpha_t` is applied as `2 * alpha` (targets) / `2 * (1 - alpha)`
    /// (decoys), not as the bare `alpha`. The textbook form would leave a
    /// constant factor of 0.5 at the neutral setting; scaling by 2 makes the
    /// neutral setting bit-identical to [`MlpLoss::Bce`] instead of merely
    /// proportional to it. That matters because the loss is compared ACROSS
    /// configurations here — the early-stopping rule and the per-fold summary
    /// both read absolute loss values — so a silent 2x on the reported number
    /// would make an `alpha`-sweep's rows incomparable. `alpha` still has its
    /// usual meaning: `> 0.5` favors targets, `< 0.5` favors decoys.
    ///
    /// # `alpha` and the per-row sample weights MULTIPLY — tune them together
    /// The MLP already receives per-row weights (`0.5` targets / `1.0` decoys,
    /// inherited from the GBM's convention via `cv::fold_weights`), and those
    /// weights multiply this `alpha_t`. They are two dials on the same quantity:
    /// setting `alpha = 0.5` does NOT mean the classes are balanced, it means
    /// they are balanced 0.5/1.0 by the sample weights alone. An `alpha` chosen
    /// as if the weights were uniform will double-count the imbalance. Tune the
    /// pair, never one in isolation.
    ///
    /// This is well-trodden in this domain rather than a novelty: Percolator
    /// adapts a hinge loss with different target and decoy weights for the same
    /// reason, so a class-asymmetric loss on top of class-asymmetric row weights
    /// is the normal shape of a target/decoy discriminant, not an exotic one.
    Focal { gamma: f32, alpha: f32 },
}

/// Numerically stable `sigma(z)`.
#[inline]
fn sigmoid(z: f32) -> f32 {
    if z >= 0.0 {
        1.0 / (1.0 + (-z).exp())
    } else {
        let e = z.exp();
        e / (1.0 + e)
    }
}

/// Numerically stable `-[y log p + (1-y) log(1-p)]` from the logit.
#[inline]
fn bce_from_logit(z: f32, y: f32) -> f32 {
    z.max(0.0) - z * y + (-z.abs()).exp().ln_1p()
}

impl MlpLoss {
    /// Weighted mean loss over the batch, writing `dL/dlogit` into `grad`.
    ///
    /// `logits` and `grad` are `[n][1]`. `y` is 0.0 (decoy) / 1.0 (target).
    /// `w` is the per-row sample weight; the mean is taken over `sum(w)` so the
    /// loss scale does not move when the class balance does.
    ///
    /// # Gradient of the focal form
    /// Write `p = sigma(z)`, `p_t = p` for a target and `1 - p` for a decoy,
    /// `q = 1 - p_t` (the wrong-class mass), `s = +1` for a target and `-1` for a
    /// decoy, `bce = -log(p_t)`, and `m = alpha_t * q^gamma`, so `L_i = m * bce`.
    ///
    /// Two standard identities do the work. First `dp_t/dz = s * p_t * q`, since
    /// `dp/dz = p(1-p) = p_t * q` and the decoy case flips the sign. Hence
    ///
    /// ```text
    ///   d(bce)/dz = -(1/p_t) * dp_t/dz = -s * q = p - y
    /// ```
    ///
    /// — the plain BCE-with-logits gradient. Second, the MODULATING FACTOR is
    /// itself a function of `z`, so it contributes its own term:
    ///
    /// ```text
    ///   dm/dz = alpha_t * gamma * q^(gamma-1) * (-dp_t/dz)
    ///         = -gamma * s * p_t * m
    /// ```
    ///
    /// Product rule:
    ///
    /// ```text
    ///   dL/dz = m * [ (p - y) - gamma * s * p_t * bce ]
    /// ```
    ///
    /// The `- gamma * s * p_t * bce` term is what makes this NOT a reweighted BCE
    /// gradient: dropping it (the easy mistake) leaves a gradient that is wrong
    /// everywhere `gamma != 0`, yet still trains, because it still points
    /// downhill. Finite differences are the only cheap check — see the tests.
    ///
    /// At `gamma = 0` the bracket collapses to `p - y` and `q^0 = 1`, so with
    /// `alpha = 0.5` (hence `alpha_t = 1`) both the loss and the gradient are the
    /// `Bce` arm's, to the bit.
    pub fn loss_and_grad(&self, logits: &Tensor, y: &[f32], w: &[f32], grad: &mut Tensor) -> f32 {
        let n = logits.rows;
        debug_assert_eq!(logits.cols, 1, "single-logit head");
        debug_assert_eq!(y.len(), n);
        debug_assert_eq!(w.len(), n);
        grad.reshape(n, 1);

        let wsum: f32 = w.iter().sum();
        let inv_w = if wsum > 0.0 { 1.0 / wsum } else { 0.0 };

        let mut total = 0.0f32;
        for b in 0..n {
            let z = logits.row(b)[0];
            let yi = y[b];
            let bce = bce_from_logit(z, yi);

            let (m, dm_term) = match *self {
                MlpLoss::Bce => (1.0, 0.0),
                MlpLoss::Focal { gamma, alpha } => {
                    let target = yi > 0.5;
                    // `p_t` and `q` are each taken from `sigmoid` DIRECTLY
                    // rather than one from `1.0 - other`: at |z| ~ 20 one of
                    // them rounds to 1.0 and the complement would come out as
                    // exactly 0 instead of the ~1e-9 it is. `sigmoid(-z)`
                    // keeps the small side's magnitude.
                    let (p_t, q) = if target {
                        (sigmoid(z), sigmoid(-z))
                    } else {
                        (sigmoid(-z), sigmoid(z))
                    };
                    // 2x, so `alpha = 0.5` is the neutral setting — see the
                    // variant docs.
                    let a_t = if target {
                        2.0 * alpha
                    } else {
                        2.0 * (1.0 - alpha)
                    };
                    let s = if target { 1.0 } else { -1.0 };
                    // `q.powf(0.0)` is exactly 1.0 for every `q`, including 0,
                    // and `0.0 * anything_finite` is 0.0, so `gamma = 0` leaves
                    // both of these at the `Bce` arm's values.
                    (a_t * q.powf(gamma), -gamma * s * p_t * bce)
                }
            };

            let p = sigmoid(z);
            total += w[b] * m * bce;
            grad.row_mut(b)[0] = w[b] * inv_w * (m * (p - yi) + m * dm_term);
        }
        total * inv_w
    }
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

    /// `n_in -> hidden[0] -> ... -> 1`, `act` between hidden layers, single
    /// logit out. No BatchNorm: inputs arrive standardized.
    ///
    /// `act` comes from [`MlpConfig::activation`] on every production path.
    /// It is a positional argument rather than a builder default so that adding
    /// [`Activation::Square`] could not leave a call site silently on the
    /// rectifier while its config said otherwise.
    pub fn feedforward(n_in: usize, hidden: &[usize], act: Activation, rng: &mut StdRng) -> Self {
        let mut layers: Vec<Box<dyn Layer>> = Vec::new();
        let mut prev = n_in;
        for &h in hidden {
            layers.push(Box::new(Linear::new(prev, h, rng)));
            layers.push(act.layer(h));
            prev = h;
        }
        layers.push(Box::new(Linear::new(prev, 1, rng)));
        Mlp::new(layers)
    }

    pub fn forward(&mut self, x: &Tensor, training: bool) -> &Tensor {
        for k in 0..self.layers.len() {
            if k == 0 {
                self.layers[0].forward(x, &mut self.acts[0], training);
            } else {
                let (prev, cur) = self.acts.split_at_mut(k);
                self.layers[k].forward(&prev[k - 1], &mut cur[0], training);
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
    pub fn train_step(
        &mut self,
        x: &Tensor,
        y: &[f32],
        w: &[f32],
        loss: &MlpLoss,
        opt: &mut Adam,
    ) -> f32 {
        self.forward(x, true);
        let last = self.acts.len() - 1;
        let l = loss.loss_and_grad(&self.acts[last], y, w, &mut self.grads[last]);
        self.backward(x);
        opt.step(self.params_and_grads());
        l
    }

    /// Mean loss over a set WITHOUT touching parameters, gradients, or the RNG.
    ///
    /// Runs the forward pass in eval mode and reuses
    /// [`MlpLoss::loss_and_grad`] rather than re-deriving the loss, so the
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
    pub fn eval_loss(&mut self, x: &Tensor, y: &[f32], w: &[f32], loss: &MlpLoss) -> f32 {
        self.forward(x, false);
        let last = self.acts.len() - 1;
        loss.loss_and_grad(&self.acts[last], y, w, &mut self.eval_sink)
    }

    /// Copy every parameter into `out`, in [`Self::params_and_grads`] order.
    ///
    /// `out` is reused across epochs, so a whole training run snapshots into one
    /// allocation. The order is a pure function of the layer stack, which never
    /// changes after construction, so [`Self::load_params_from`] can walk it
    /// back with nothing but the length of each slot.
    ///
    /// **Parameters only.** [`BatchNorm`]'s running mean/variance are NOT
    /// parameters and are therefore NOT snapshotted, so restoring a snapshot
    /// into a net containing one would pair epoch-`k` weights with epoch-`n`
    /// running stats. No architecture this crate builds contains a `BatchNorm`
    /// ([`Self::feedforward`] emits `Linear` plus a parameterless activation —
    /// [`Relu`] or [`Square`] — only) and the type is kept for gradient tests, so
    /// this is documented rather than solved.
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

    /// Shuffled-minibatch training. Returns the mean loss of the final epoch.
    ///
    /// The per-epoch shuffle is seeded from the caller's RNG (which
    /// [`MlpConfig::rng_for_fold`] derives from the config seed and the fold
    /// index), so folds differ but a rerun of the same fold does not.
    ///
    /// Minibatching is not only a throughput choice: full-batch descent on a
    /// ReLU net can settle into a dead-unit configuration and stay there,
    /// because every step sees the identical gradient. The per-batch noise is
    /// what walks it back out.
    ///
    /// No held-out set, therefore NO EARLY STOPPING regardless of
    /// [`MlpConfig::early_stopping_patience`]: the stopping rule has nothing to
    /// measure. Callers that want it go through [`Self::train_reporting`].
    pub fn train(
        &mut self,
        cfg: &MlpConfig,
        x: &Tensor,
        y: &[f32],
        w: &[f32],
        opt: &mut Adam,
        rng: &mut StdRng,
    ) -> f32 {
        self.train_reporting(cfg, x, y, w, opt, rng, None, "")
            .final_train_loss
    }

    /// [`Self::train`] plus a per-epoch loss trace at `debug` level and,
    /// when [`MlpConfig::early_stopping_patience`] is set and `val` is present,
    /// EARLY STOPPING with best-weight restore.
    ///
    /// `val` is a HELD-OUT set. It reaches [`Self::eval_loss`] and nothing
    /// else — never the optimizer, never the input transform, never the RNG —
    /// so the only ways it can influence the result are the log and the choice
    /// of WHICH epoch's weights are kept. Both are what early stopping is.
    ///
    /// Enable the trace with:
    /// ```text
    /// RUST_LOG=timsseek::ml::mlp=debug
    /// ```
    ///
    /// The point of tracing BOTH curves is that they answer different
    /// questions. Train loss still falling means the fixed epoch budget is the
    /// binding constraint; held-out loss flattening or turning up while train
    /// loss keeps falling means the budget is already too large and the epochs
    /// after the turn are spent overfitting. Only the second reading justifies
    /// early stopping, and it cannot be read off the train curve alone.
    ///
    /// # The stopping rule
    /// After every epoch the held-out loss is measured. An epoch IMPROVES iff
    /// `loss < best` — a STRICT comparison, so a tie is not an improvement.
    /// Two consequences, both wanted: the EARLIEST epoch attaining the minimum
    /// is the one kept (fewer epochs, and no dependence on how a later epoch's
    /// rounding happens to land), and a dead-flat plateau still runs the
    /// patience counter down instead of training forever. A `NaN` held-out loss
    /// compares false and so counts as no improvement; if EVERY epoch is `NaN`
    /// there is no snapshot to restore and the last weights stand.
    ///
    /// When `patience` epochs pass with no improvement, the best-seen weights
    /// are restored and training stops. Restoring is not optional: stopping at
    /// `best + patience` and keeping those weights leaves exactly the
    /// overfitted parameters the rule exists to avoid. The restore also happens
    /// when the epoch budget runs out before patience does, so with early
    /// stopping on you always get the best-held-out-loss weights that were seen.
    #[allow(clippy::too_many_arguments)]
    pub fn train_reporting(
        &mut self,
        cfg: &MlpConfig,
        x: &Tensor,
        y: &[f32],
        w: &[f32],
        opt: &mut Adam,
        rng: &mut StdRng,
        val: Option<ValSet<'_>>,
        tag: &str,
    ) -> TrainOutcome {
        let n = x.rows;
        debug_assert_eq!(y.len(), n);
        debug_assert_eq!(w.len(), n);

        // Early stopping needs something to stop ON. With no held-out set the
        // knob is inert, which is what makes `Mlp::train` and the `val = &[]`
        // paths behave as they always did.
        let patience = match (cfg.early_stopping_patience, val.is_some()) {
            (Some(p), true) => Some(p),
            _ => None,
        };

        let mut order: Vec<usize> = (0..n).collect();
        let batch = cfg.batch_size.max(1);
        let (mut xb, mut yb, mut wb) = (Tensor::default(), Vec::new(), Vec::new());
        let mut last_epoch = 0.0;

        let mut best_val = f32::INFINITY;
        let mut best_epoch: Option<usize> = None;
        let mut best_train = 0.0f32;
        let mut best_params: Vec<f32> = Vec::new();
        let mut stale = 0usize;
        let mut epochs_run = 0usize;

        for epoch in 0..cfg.epochs {
            // Fisher-Yates, so the shuffle consumes a fixed number of draws
            // per epoch and the sequence stays reproducible.
            for i in (1..order.len()).rev() {
                let j = rng.random_range(0..=i);
                order.swap(i, j);
            }

            let (mut running, mut nb) = (0.0, 0);
            for chunk in order.chunks(batch) {
                xb.reshape(chunk.len(), x.cols);
                yb.clear();
                wb.clear();
                for (bi, &i) in chunk.iter().enumerate() {
                    xb.row_mut(bi).copy_from_slice(x.row(i));
                    yb.push(y[i]);
                    wb.push(w[i]);
                }
                running += self.train_step(&xb, &yb, &wb, &cfg.loss, opt);
                nb += 1;
            }
            last_epoch = if nb > 0 { running / nb as f32 } else { 0.0 };
            epochs_run = epoch + 1;

            // The held-out forward pass runs when the STOPPING RULE needs it
            // (every epoch, because it is the decision) or when the trace is on
            // (a diagnostic nobody asked to see should not cost a forward pass).
            // `tracing::enabled!` is still what gates the LOG, exactly as before.
            let traced = tracing::enabled!(tracing::Level::DEBUG);
            let val_loss = match (&val, patience.is_some() || traced) {
                (Some(v), true) => Some(self.eval_loss(v.x, v.y, v.w, &cfg.loss)),
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

            let Some(p) = patience else { continue };
            // `val_loss` is `Some` whenever `patience` is: both arms of the
            // match above require `val` to be present, and `patience` is only
            // `Some` when it is.
            let vl = val_loss.expect("early stopping requires a held-out measurement");
            // STRICT: a tie is not an improvement. See the rule in the docs.
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

        // Restore whenever the kept epoch is not the one we happen to be
        // standing on — i.e. after a patience stop, and also after a budget
        // exhaustion whose last epochs were not improvements.
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

/// A held-out set for the per-epoch loss trace AND the early-stopping
/// decision, borrowed for the call.
///
/// Deliberately a distinct type from the training triple so a caller cannot
/// swap the two by argument order: handing the training set here would turn
/// early stopping into "stop when the train loss stops falling", which on this
/// net is never, and handing the held-out set to `train` would be an actual
/// leak. The rows here reach nothing but [`Mlp::eval_loss`] — no optimizer
/// step, no input transform, no RNG draw — so the only thing they can move is
/// which epoch's weights are kept.
pub struct ValSet<'a> {
    pub x: &'a Tensor,
    pub y: &'a [f32],
    pub w: &'a [f32],
}

/// What one [`Mlp::train_reporting`] call did.
///
/// A struct rather than the bare loss because with early stopping the loss no
/// longer describes the run: the caller's log wants to say WHERE it stopped and
/// which epoch it kept, and the tests need those two numbers to distinguish
/// "stopped early" from "ran the budget and got lucky".
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
    pub fn fit(feat: &[f64], ncols_lane: usize, rows: &[usize]) -> Self {
        let mut sum = vec![0.0f64; ncols_lane];
        let mut sumsq = vec![0.0f64; ncols_lane];
        let mut finite = vec![0u64; ncols_lane];
        let mut nonfinite = vec![0u64; ncols_lane];

        for &i in rows {
            let row = &feat[i * ncols_lane..(i + 1) * ncols_lane];
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

/// Hyperparameters. Not TOML-exposed, matching how `GBMConfig` is handled —
/// only the input dimension is dynamic, and it is derived from the lane matrix
/// at fit time rather than configured.
///
/// # Dev-only env overrides
/// [`Self::default`] applies the `TIMSSEEK_MLP_*` overrides described on
/// [`Self::from_env`]. That is an EXPERIMENT ESCAPE HATCH for sweeps, in the
/// spirit of `TIMSSEEK_LDA_DUMP`, not a config surface: it is undocumented in
/// `--help`, unsupported, and with nothing set the config is bit-identical to
/// [`Self::compiled_default`].
#[derive(Debug, Clone, PartialEq)]
pub struct MlpConfig {
    pub hidden: Vec<usize>,
    /// Adam step size. COUPLED TO [`Self::batch_size`] — see
    /// [`Self::compiled_default`], where the pair was measured together and
    /// moving one alone is documented as the worst arm of the sweep.
    ///
    /// # A superseded warning, kept because it is still true of its own regime
    /// This field used to carry "1e-2/1e-3 both walked a 2-input ReLU net into a
    /// dead-unit plateau in testing", which is why the default was 3e-4 for so
    /// long. That observation stands, but it is about the XOR fixture in this
    /// module's tests — 2 inputs, pure ReLU, few steps — and NOT about the 101-
    /// or 128-column rescore lanes, where 1.2e-3 was measured better on two real
    /// searches. Do not read it as a ceiling on this field. See [`Relu`] for the
    /// plateau's actual mechanism, which a learning rate cannot fix anyway.
    pub lr: f32,
    pub weight_decay: f32,
    /// UPPER BOUND on epochs, not a target: with
    /// [`Self::early_stopping_patience`] set, training normally stops before
    /// this.
    pub epochs: usize,
    /// Minibatch size. COUPLED TO [`Self::lr`]: it sets how many optimizer steps
    /// an epoch is worth, so changing it without rescaling the learning rate
    /// changes the amount of training, not just its granularity. See
    /// [`Self::compiled_default`] for the measurement.
    pub batch_size: usize,
    pub loss: MlpLoss,
    /// Activation between hidden layers. Defaults to
    /// [`Activation::LeakyRelu`], which is what every number quoted on
    /// [`Self::compiled_default`] was measured with; [`Activation::Square`] is an
    /// experiment, and [`Square`] documents both what recommends it and the
    /// divergence hazard it brings.
    ///
    /// NOT INDEPENDENT OF [`Self::lr`] AND [`Self::weight_decay`]. `x^2` can
    /// amplify its input where a rectifier cannot, so the step size that is right
    /// for one is not automatically right for the other. Changing this field alone
    /// and reading the result as "the activation is worse" is the same mistake the
    /// `lr`/`batch_size` coupling above warns about.
    pub activation: Activation,
    pub seed: u64,
    /// Stop after this many consecutive epochs with no improvement in held-out
    /// loss, and restore the best-seen weights. `None` disables the whole
    /// mechanism: the full `epochs` budget runs, nothing is snapshotted, and
    /// the fit is bit-identical to what it was before early stopping existed.
    ///
    /// # Why it defaults to on
    /// Measured on a 114k-candidate search, all three rescore folds bottom out
    /// well inside the then-30-epoch budget and then get worse, while the TRAIN
    /// loss falls monotonically to the end:
    ///
    /// | fold | best held-out | at epoch | held-out @30 |
    /// |------|---------------|----------|--------------|
    /// | 0    | 0.465810      | 16       | 0.469477     |
    /// | 1    | 0.459646      | 23       | 0.461799     |
    /// | 2    | 0.458787      | 17       | 0.462773     |
    ///
    /// THIS TABLE IS PRE-RETUNE: it was taken at `hidden: [64, 32]`, `lr: 3e-4`,
    /// `batch_size: 256`, `epochs: 30`, none of which are the current defaults
    /// (see [`Self::compiled_default`]). Where the turn now falls has NOT been
    /// re-measured. What it establishes is the MECHANISM, which is not
    /// hyperparameter-specific: held-out loss turns while train loss keeps
    /// falling, so 40-47% of that budget was spent memorizing.
    ///
    /// Lowering `epochs` instead would NOT do the same job: the best epoch is 16,
    /// 23 and 17 on the three folds of one single run, so any fixed budget is too
    /// long for one fold and too short for another. The turn is a per-fold
    /// property and only a per-fold rule can find it.
    ///
    /// # Why 5
    /// The three curves above are within ~0.2% of their minimum for several
    /// epochs around it, so a patience of 1 or 2 would stop on measurement noise
    /// well before the turn. Five is long enough to walk through a flat stretch
    /// and short enough to cut most of the memorizing tail. It has NOT been swept
    /// on real data; it is a default, not a tuned value.
    pub early_stopping_patience: Option<usize>,
}

impl Default for MlpConfig {
    /// [`Self::compiled_default`], plus the dev-only `TIMSSEEK_MLP_*` overrides
    /// (see [`Self::from_env`]) in non-test builds.
    ///
    /// # Why `cfg(test)` skips the environment
    /// This crate's tests run in-process and share one environment, and a large
    /// number of them build their fixture config as `MlpConfig { .. ,
    /// ..MlpConfig::default() }` — including the seed-sweeping learning tests,
    /// whose whole job is to pin training behaviour at KNOWN hyperparameters. A
    /// stray `TIMSSEEK_MLP_EPOCHS` in the developer's shell would silently
    /// retune every one of them at once, so under `cfg(test)` `default()` is
    /// exactly the compiled default and the environment is never read.
    /// [`Self::apply_overrides`] is what the parsing tests exercise instead: it
    /// takes the lookup as an argument, so the whole override path is testable
    /// without a process-wide (and, in edition 2024, `unsafe`) `set_var`.
    fn default() -> Self {
        #[cfg(test)]
        {
            Self::compiled_default()
        }
        #[cfg(not(test))]
        {
            Self::from_env()
        }
    }
}

/// Gate for the once-per-process effective-config log in [`MlpConfig::from_env`].
/// `default()` is called once per rescorer entry point, but the hybrids and the
/// tests can reach it more than once in a process and the sweep wants ONE line
/// to grep, not one per call.
static ENV_OVERRIDE_LOGGED: std::sync::Once = std::sync::Once::new();

impl MlpConfig {
    /// The hyperparameters AS COMPILED. No environment, no logging.
    ///
    /// This is the authoritative default — the values here are what a run with
    /// no `TIMSSEEK_MLP_*` set trains with, and every override starts from this
    /// struct and replaces individual fields.
    /// # Provenance of these numbers
    /// `hidden`, `lr`, `batch_size` and `epochs` are the winning arm of a sweep
    /// on two REAL searches (`mlp_all`, Phase 5 wall time / targets at 1% FDR):
    /// a 114138-candidate run and a 2174837-candidate run of the same Astral
    /// mzML. They replaced conventional, unmeasured values. Against the GBM the
    /// winner is 1.31x faster with +211 targets at 114k, and 2.04x faster with
    /// -55 targets (-0.04%) at 2.17M; against the previous MLP defaults it is
    /// 1.73x-2.26x faster while finding MORE targets at both scales.
    ///
    /// Sweep-wide caveats, which apply to every figure quoted on the fields
    /// below: both runs are the same FAIMS acquisition, so there is no ion
    /// mobility and 25 of the 101 linear / 27 of the 128 all-lane columns were
    /// culled as dead; two datasets, one machine; and Phase 5 is a fraction of a
    /// whole search, so a 2x here is not a 2x end to end. See
    /// [`crate::ml::RescoreModel`] for the cross-model table.
    pub fn compiled_default() -> Self {
        MlpConfig {
            // MEASURED, and counter-intuitively so: `[64, 32]` was
            // OVERPARAMETERIZED for this data. Halving it to `[32, 16]` improved
            // speed AND accuracy at both scales. The two arms compared, both
            // named in the sweep ledger -- baseline (`hidden: [64, 32]`,
            // `epochs: 30`, batch 256, `lr` 3e-4) against arm "D"
            // (`hidden: [32, 16]`, `epochs: 60`, batch and `lr` left at those old
            // values): at 114k candidates 16.614s/25425 -> 9.165s/25476, and at
            // 2.17M 487.042s/134345 -> 395.990s/134302 (-43 targets, -0.03%, for
            // a 1.23x speedup at that scale).
            //
            // THAT PAIR IS CONFOUNDED, and no arm of the sweep isolated `hidden`
            // on its own: `epochs` moved from 30 to 60 alongside it. The
            // conclusion survives because the epoch raise is INERT -- the three
            // folds instrumented (the 114k run) stopped at 21, 28 and 22 epochs,
            // under either ceiling, so a raise from 30 to 60 cannot be what moved
            // the time (see `epochs` below). Read the pair as "`hidden` plus a
            // knob that did nothing", not as a clean single-variable measurement,
            // and note that only that one run's stop epochs were recorded: at
            // 2.17M nothing counted them, so inertness there is an inference.
            //
            // Either way, this is NOT a chosen point on a speed/accuracy
            // trade-off curve; the smaller net is simply better here, so "it is
            // undersized, give it capacity back" is a change that has already
            // been tried and lost.
            hidden: vec![32, 16],
            // COUPLED TO `batch_size`. NEVER MOVE ONE WITHOUT THE OTHER: this is
            // the old 3e-4 (which was tuned at batch 256) scaled linearly with
            // the 4x larger batch. Raising the batch alone -- 1024 with `lr`
            // left at 3e-4 -- was the WORST arm of the entire sweep: 31.612s on
            // the 114k run against 16.614s for changing nothing at all, because
            // 4x fewer optimizer steps per epoch never converges and spends the
            // whole epoch budget not converging. Lowering `lr` instead is worse
            // still: 1e-4 over 100 epochs cost 35.512s for fewer targets.
            lr: 1.2e-3,
            weight_decay: 1e-4,
            // A CEILING THAT EARLY STOPPING IS EXPECTED TO CUT SHORT, not a
            // target -- see [`Self::early_stopping_patience`]. Raising it does
            // not mean 60 epochs are wanted or expected to run.
            //
            // THE BUDGET WAS NEVER THE BINDING CONSTRAINT, at 30 or at 60. In the
            // one run whose stop epochs were recorded (114k candidates, `epochs:
            // 30`, and the pre-retune `hidden`/`lr`/`batch_size` -- nothing has
            // counted epochs since) the folds ran 21, 28 and 22, which is exactly
            // `best + patience` for the measured best epochs 16, 23 and 17: all
            // three terminated ON THE CRITERION and 30 truncated nothing. (An
            // earlier version of this comment claimed the opposite. It was wrong,
            // and the speed win came from `hidden`/`batch_size`/`lr`, not from the
            // ceiling.)
            //
            // What the 28-of-30 fold does leave is HEADROOM, which is the only
            // argument for 60: two spare epochs is not enough slack to tell
            // whether that fold's held-out curve had finished turning, and a fold
            // that turns later would be cut off rather than measured. 60 removes
            // that doubt.
            //
            // NO DIRECT A/B SUPPORTS 60 at the `hidden`/`lr`/`batch_size` above;
            // it rests on the headroom argument alone. It is also free whenever
            // early stopping fires, and it fired on every fold ever instrumented,
            // so being generous here has cost nothing observed. Treat 60 as
            // UNVALIDATED -- neither known-good nor known-wrong.
            epochs: 60,
            // COUPLED TO `lr` -- see there for the measurement. 1024 pays ONLY
            // with the learning rate scaled to match (256 -> 1024 is 4x, so
            // 3e-4 -> 1.2e-3). Swept as an independent axis it looks strictly
            // harmful, and anyone who tunes it that way will conclude, wrongly,
            // that bigger batches are slower on this workload.
            batch_size: 1024,
            loss: MlpLoss::Bce,
            // The rectifier the whole sweep above ran on. See `activation`.
            activation: Activation::LeakyRelu,
            seed: 0x2545_F491_4F6C_DD1D,
            early_stopping_patience: Some(5),
        }
    }

    /// [`Self::compiled_default`] with the dev-only `TIMSSEEK_MLP_*` overrides
    /// applied. **THE ONLY PLACE IN THE CRATE THAT READS THE ENVIRONMENT FOR
    /// HYPERPARAMETERS** — one reader is what makes "which knobs exist" a
    /// question with a single answer, in [`Self::apply_overrides`].
    ///
    /// This exists because `MlpConfig` is deliberately NOT TOML-exposed (see the
    /// type docs), which makes every hyperparameter experiment a recompile —
    /// unworkable for sweeping a 2M-candidate search. It is a dev affordance in
    /// the spirit of `TIMSSEEK_LDA_DUMP`, NOT a supported interface, and it is
    /// not a step towards putting these in the config file.
    ///
    /// | Variable | Field | Format |
    /// |---|---|---|
    /// | `TIMSSEEK_MLP_EPOCHS` | `epochs` | positive integer, e.g. `60` |
    /// | `TIMSSEEK_MLP_LR` | `lr` | positive finite float, e.g. `1e-3` |
    /// | `TIMSSEEK_MLP_BATCH_SIZE` | `batch_size` | positive integer, e.g. `512` |
    /// | `TIMSSEEK_MLP_HIDDEN` | `hidden` | comma-separated positive integers (`64,32`), or `none` for no hidden layer |
    /// | `TIMSSEEK_MLP_WEIGHT_DECAY` | `weight_decay` | non-negative finite float, e.g. `0` or `1e-3` |
    /// | `TIMSSEEK_MLP_PATIENCE` | `early_stopping_patience` | positive integer, or `none`/`off` to disable early stopping entirely |
    /// | `TIMSSEEK_MLP_ACTIVATION` | `activation` | `leaky_relu` or `square`, case-insensitive (`-` and `_` interchangeable) |
    /// | `TIMSSEEK_MLP_SEED` | `seed` | `u64`, decimal or `0x`-prefixed hex |
    /// | `TIMSSEEK_MLP_LOSS` | `loss` | `bce`, or `focal:GAMMA:ALPHA` (`focal:2:0.25`) |
    ///
    /// `TIMSSEEK_MLP_LOSS` spells BOTH focal parameters or none of them: they are
    /// coupled to each other and to the per-row sample weights (see
    /// [`MlpLoss::Focal`]), so there is no defensible default for one given the
    /// other, and a bare `focal` that silently picked a `gamma` would be the
    /// mislabelled sweep row this whole mechanism aborts to avoid. This is still
    /// a sweep hatch and not the real config surface `MlpLoss` eventually wants.
    ///
    /// # Unset changes nothing
    /// A variable that is not set is not read into anything — see
    /// [`Self::apply_overrides`], which touches a field only on `Some`. With
    /// none set this returns [`Self::compiled_default`] and logs nothing.
    ///
    /// # A malformed value ABORTS
    /// `TIMSSEEK_MLP_LR=banana` panics with the variable name, the offending
    /// value and the expected format. It does NOT warn-and-continue, because the
    /// failure this override exists inside of is a SWEEP: a warned-past variable
    /// produces a run that completes normally, lands in the results table under
    /// the label the operator believes they set, and is indistinguishable from
    /// the compiled-default row next to it. A `warn!` line is 1 of ~10^4 in a
    /// full search log and would be found after the conclusion was drawn.
    /// Aborting costs one re-run and cannot be missed.
    ///
    /// `timsseek_cli` calls this at startup when the selected rescore model uses
    /// the MLP, so the abort lands in the first second rather than at Phase 5.
    ///
    /// # Determinism
    /// Untouched: this only chooses the numbers. For a FIXED environment the fit
    /// is the same pure function of (config, fold, rows) it was before — see
    /// [`Self::rng_for_fold`].
    pub fn from_env() -> Self {
        let mut cfg = Self::compiled_default();
        let applied = cfg.apply_overrides(|k| std::env::var(k).ok());
        if !applied.is_empty() {
            // Once per process, at `info`: the default filter is `info`, so a
            // sweep's run log records WHAT ACTUALLY RAN with no extra flags,
            // which is the whole point. Not `warn` — deliberately setting an
            // override is not a problem, and the malformed case aborts instead.
            ENV_OVERRIDE_LOGGED.call_once(|| {
                tracing::info!(
                    "MLP config: DEV OVERRIDE ACTIVE ({}); effective config: hidden={:?} \
                     lr={:e} weight_decay={:e} epochs={} batch_size={} patience={:?} \
                     activation={:?} seed={:#018x} loss={:?}",
                    applied.join(" "),
                    cfg.hidden,
                    cfg.lr,
                    cfg.weight_decay,
                    cfg.epochs,
                    cfg.batch_size,
                    cfg.early_stopping_patience,
                    cfg.activation,
                    cfg.seed,
                    cfg.loss,
                );
            });
        }
        cfg
    }

    /// The pure core of [`Self::from_env`]: `get` stands in for
    /// `std::env::var`, so the parsing, the validation and the "unset changes
    /// nothing" property are all testable without mutating the process
    /// environment (`unsafe` in edition 2024, and racy across parallel tests
    /// regardless).
    ///
    /// Returns one `NAME=value` string per variable that was actually read, in
    /// declaration order — the audit trail [`Self::from_env`] logs. Empty means
    /// nothing was touched.
    ///
    /// # Panics
    /// On any malformed or out-of-range value. See [`Self::from_env`].
    fn apply_overrides(&mut self, get: impl Fn(&str) -> Option<String>) -> Vec<String> {
        let mut applied = Vec::new();
        // One closure per variable so a variable that is UNSET never reaches an
        // assignment: `read` yields `None` and the `if let` body is skipped.
        let mut read = |name: &str| -> Option<String> {
            let raw = get(name)?;
            // Quoted, so a value that contains a space or is empty is still
            // unambiguous in the one-line audit trail.
            applied.push(format!("{name}={raw:?}"));
            Some(raw)
        };

        if let Some(raw) = read("TIMSSEEK_MLP_EPOCHS") {
            self.epochs = parse_positive_usize("TIMSSEEK_MLP_EPOCHS", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_LR") {
            self.lr = parse_float("TIMSSEEK_MLP_LR", &raw, false);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_BATCH_SIZE") {
            self.batch_size = parse_positive_usize("TIMSSEEK_MLP_BATCH_SIZE", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_HIDDEN") {
            self.hidden = parse_hidden("TIMSSEEK_MLP_HIDDEN", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_WEIGHT_DECAY") {
            self.weight_decay = parse_float("TIMSSEEK_MLP_WEIGHT_DECAY", &raw, true);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_PATIENCE") {
            self.early_stopping_patience = parse_patience("TIMSSEEK_MLP_PATIENCE", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_ACTIVATION") {
            self.activation = parse_activation("TIMSSEEK_MLP_ACTIVATION", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_SEED") {
            self.seed = parse_seed("TIMSSEEK_MLP_SEED", &raw);
        }
        if let Some(raw) = read("TIMSSEEK_MLP_LOSS") {
            self.loss = parse_loss("TIMSSEEK_MLP_LOSS", &raw);
        }

        applied
    }

    /// Per-fold RNG: the configured seed mixed with the fold index, so folds
    /// differ but a rerun of the same fold does not.
    pub fn rng_for_fold(&self, fold: usize) -> StdRng {
        StdRng::seed_from_u64(self.seed ^ (fold as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }
}

/// The one abort message shape, so every malformed variable explains itself the
/// same way and says why it is not continuing.
fn bad_env(name: &str, raw: &str, expected: &str) -> ! {
    panic!(
        "{name}={raw:?} is not usable: expected {expected}.\n\
         This is a DEV-ONLY MLP hyperparameter override and a malformed value ABORTS rather \
         than falling back to the compiled default: a run that quietly used the default would \
         land in your sweep table labelled with a value it never trained on. Fix the variable \
         or unset it."
    );
}

fn parse_positive_usize(name: &str, raw: &str) -> usize {
    match raw.trim().parse::<usize>() {
        Ok(v) if v > 0 => v,
        _ => bad_env(name, raw, "a positive integer, e.g. `60`"),
    }
}

/// `allow_zero` distinguishes `weight_decay` (0 disables it, a legitimate sweep
/// arm) from `lr` (0 trains nothing and is never what anyone meant).
fn parse_float(name: &str, raw: &str, allow_zero: bool) -> f32 {
    let lo = if allow_zero {
        "non-negative"
    } else {
        "positive"
    };
    match raw.trim().parse::<f32>() {
        Ok(v) if v.is_finite() && (v > 0.0 || (allow_zero && v == 0.0)) => v,
        _ => bad_env(name, raw, &format!("a finite {lo} float, e.g. `1e-3`")),
    }
}

fn parse_hidden(name: &str, raw: &str) -> Vec<usize> {
    let t = raw.trim();
    // No hidden layer is a real sweep arm (`feedforward` then builds a single
    // `n_in -> 1` linear layer), but it must be SPELLED: an empty string is
    // indistinguishable from a shell that dropped the value.
    if t.eq_ignore_ascii_case("none") {
        return Vec::new();
    }
    let mut out = Vec::new();
    for part in t.split(',') {
        match part.trim().parse::<usize>() {
            Ok(v) if v > 0 => out.push(v),
            _ => bad_env(
                name,
                raw,
                "a comma-separated list of positive integers, e.g. `64,32`, or `none` for no \
                 hidden layer",
            ),
        }
    }
    if out.is_empty() {
        bad_env(
            name,
            raw,
            "a comma-separated list of positive integers, e.g. `64,32`, or `none` for no hidden \
             layer",
        );
    }
    out
}

fn parse_patience(name: &str, raw: &str) -> Option<usize> {
    let t = raw.trim();
    if t.eq_ignore_ascii_case("none") || t.eq_ignore_ascii_case("off") {
        return None;
    }
    match t.parse::<usize>() {
        // 0 is rejected rather than mapped to either meaning: as a number it
        // would stop on the first non-improving epoch, and as a synonym for
        // "off" it would silently disagree with the `none` spelling.
        Ok(v) if v > 0 => Some(v),
        _ => bad_env(
            name,
            raw,
            "a positive integer, or `none`/`off` to disable early stopping",
        ),
    }
}

/// `leaky_relu` / `square`. Case is ignored and `-` normalizes to `_`, so the
/// three spellings a shell script is likely to produce (`SQUARE`, `leaky-relu`,
/// `leaky_relu`) all land; anything else ABORTS like every other malformed
/// override, because "the activation I asked for was ignored" is invisible in a
/// sweep table.
fn parse_activation(name: &str, raw: &str) -> Activation {
    let t = raw.trim().to_ascii_lowercase().replace('-', "_");
    match t.as_str() {
        "leaky_relu" => Activation::LeakyRelu,
        "square" => Activation::Square,
        _ => bad_env(name, raw, "`leaky_relu` or `square` (case-insensitive)"),
    }
}

fn parse_seed(name: &str, raw: &str) -> u64 {
    let t = raw.trim();
    // Hex accepted because the compiled default is written that way, so
    // `TIMSSEEK_MLP_SEED=0x2545F4914F6CDD1D` reproduces it verbatim.
    let parsed = match t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        Some(hex) => u64::from_str_radix(hex, 16),
        None => t.parse::<u64>(),
    };
    parsed.unwrap_or_else(|_| {
        bad_env(
            name,
            raw,
            "a u64, decimal (`12345`) or hex (`0x2545F4914F6CDD1D`)",
        )
    })
}

/// `bce` / `focal:GAMMA:ALPHA`. Case is ignored and whitespace around each field
/// is trimmed, so `FOCAL: 2 : 0.25` lands like `focal:2:0.25`.
///
/// The focal arm takes EXACTLY two numbers, both spelled out: `focal` on its own
/// aborts rather than filling in a `gamma`, because the two parameters are
/// coupled to each other and to the per-row sample weights (see
/// [`MlpLoss::Focal`]) and an invented value would land in the sweep table under
/// the operator's label.
///
/// # The ranges
/// `gamma >= 0` (finite): 0 disables the hard-example modulation, negative would
/// invert it and up-weight the easy bulk. `alpha` strictly inside `(0, 1)`
/// (finite): it is applied normalized, as `2 * alpha` for targets and
/// `2 * (1 - alpha)` for decoys, so 0 or 1 would multiply one entire class by
/// zero — training on one label, which is never what an `alpha` sweep means.
/// `focal:0:0.5` is the neutral point and is bit-identical to `bce`.
fn parse_loss(name: &str, raw: &str) -> MlpLoss {
    const EXPECTED: &str = "`bce`, or `focal:GAMMA:ALPHA` with a finite GAMMA >= 0 and a finite \
                            ALPHA strictly inside (0, 1) — e.g. `focal:2:0.25` (case-insensitive)";
    let t = raw.trim().to_ascii_lowercase();
    if t == "bce" {
        return MlpLoss::Bce;
    }
    // The trailing `None` is what makes the field count exact: `focal:2` and
    // `focal:2:0.25:9` both fall through to the abort.
    let mut fields = t.split(':').map(str::trim);
    let (gamma, alpha) = match (fields.next(), fields.next(), fields.next(), fields.next()) {
        (Some("focal"), Some(g), Some(a), None) => (g, a),
        _ => bad_env(name, raw, EXPECTED),
    };
    let gamma = match gamma.parse::<f32>() {
        Ok(v) if v.is_finite() && v >= 0.0 => v,
        _ => bad_env(name, raw, EXPECTED),
    };
    let alpha = match alpha.parse::<f32>() {
        Ok(v) if v.is_finite() && v > 0.0 && v < 1.0 => v,
        _ => bad_env(name, raw, EXPECTED),
    };
    MlpLoss::Focal { gamma, alpha }
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

    /// Central-difference step. `f32` params with a `1e-2` step is the sweet
    /// spot: smaller and catastrophic cancellation dominates, larger and the
    /// second-order term does.
    const H: f32 = 1e-2;
    const TOL: f32 = 2e-2;

    fn seeded() -> StdRng {
        StdRng::seed_from_u64(42)
    }

    fn filled(rows: usize, cols: usize, rng: &mut StdRng) -> Tensor {
        let mut t = Tensor::new(rows, cols);
        for v in t.data.iter_mut() {
            *v = normal(rng);
        }
        t
    }

    /// Scalar objective `L = sum(y .* gy)` for a fixed upstream gradient `gy`,
    /// which makes `dL/dx` exactly what `backward` should write into `gx`.
    fn scalar_loss(layer: &mut dyn Layer, x: &Tensor, gy: &Tensor) -> f32 {
        let mut out = Tensor::default();
        layer.forward(x, &mut out, true);
        out.live().iter().zip(gy.live()).map(|(a, b)| a * b).sum()
    }

    fn check_input_grad(layer: &mut dyn Layer, mut x: Tensor, gy: &Tensor, what: &str) {
        let mut out = Tensor::default();
        layer.forward(&x, &mut out, true);
        let mut gx = Tensor::default();
        layer.backward(&x, &out, gy, &mut gx);
        let analytic = gx.live().to_vec();

        for k in 0..x.rows * x.cols {
            let orig = x.data[k];
            x.data[k] = orig + H;
            let up = scalar_loss(layer, &x, gy);
            x.data[k] = orig - H;
            let dn = scalar_loss(layer, &x, gy);
            x.data[k] = orig;
            let numeric = (up - dn) / (2.0 * H);
            let scale = analytic[k].abs().max(numeric.abs()).max(1.0);
            assert!(
                (analytic[k] - numeric).abs() / scale < TOL,
                "{what}: dL/dx[{k}] analytic {} vs numeric {}",
                analytic[k],
                numeric
            );
        }
    }

    #[test]
    fn linear_input_gradient_matches_finite_differences() {
        let mut rng = seeded();
        let mut lin = Linear::new(5, 3, &mut rng);
        let x = filled(4, 5, &mut rng);
        let gy = filled(4, 3, &mut rng);
        check_input_grad(&mut lin, x, &gy, "Linear");
    }

    #[test]
    fn linear_weight_gradient_matches_finite_differences() {
        let mut rng = seeded();
        let mut lin = Linear::new(4, 3, &mut rng);
        let x = filled(6, 4, &mut rng);
        let gy = filled(6, 3, &mut rng);

        let mut out = Tensor::default();
        lin.forward(&x, &mut out, true);
        let mut gx = Tensor::default();
        lin.backward(&x, &out, &gy, &mut gx);
        let analytic = lin.gw.clone();

        for k in 0..lin.w.len() {
            let orig = lin.w[k];
            lin.w[k] = orig + H;
            let up = scalar_loss(&mut lin, &x, &gy);
            lin.w[k] = orig - H;
            let dn = scalar_loss(&mut lin, &x, &gy);
            lin.w[k] = orig;
            let numeric = (up - dn) / (2.0 * H);
            let scale = analytic[k].abs().max(numeric.abs()).max(1.0);
            assert!(
                (analytic[k] - numeric).abs() / scale < TOL,
                "Linear: dL/dw[{k}] analytic {} vs numeric {}",
                analytic[k],
                numeric
            );
        }
    }

    #[test]
    fn relu_gradient_matches_finite_differences() {
        for slope in [0.0f32, LEAKY_SLOPE, 0.1] {
            let mut rng = seeded();
            let mut relu = Relu::leaky(4, slope);
            // Offset away from 0 so no element sits on the kink, where the
            // one-sided derivatives legitimately disagree.
            let mut x = filled(5, 4, &mut rng);
            for v in x.data.iter_mut() {
                *v += if *v >= 0.0 { 1.0 } else { -1.0 };
            }
            let gy = filled(5, 4, &mut rng);
            check_input_grad(&mut relu, x, &gy, &format!("Relu(slope={slope})"));
        }
    }

    /// The property that makes the leaky form worth defaulting to: a saturated
    /// unit still passes gradient, so it can come back. A pure ReLU cannot,
    /// which is what makes its dead state absorbing rather than merely slow.
    #[test]
    fn leaky_relu_passes_gradient_where_pure_relu_does_not() {
        let gy = {
            let mut t = Tensor::new(1, 1);
            t.row_mut(0)[0] = 1.0;
            t
        };
        let mut x = Tensor::new(1, 1);
        x.row_mut(0)[0] = -5.0; // deep in the dead region

        let grad_at = |slope: f32| {
            let mut act = Relu::leaky(1, slope);
            let mut out = Tensor::default();
            act.forward(&x, &mut out, true);
            let mut gx = Tensor::default();
            act.backward(&x, &out, &gy, &mut gx);
            gx.row(0)[0]
        };

        assert_eq!(grad_at(0.0), 0.0, "pure ReLU: no way back");
        assert_eq!(grad_at(LEAKY_SLOPE), LEAKY_SLOPE);
    }

    /// `dL/dx = gy * 2x`. Unlike [`Relu`], this needs the layer INPUT, so the
    /// check is also what would catch a `backward` that tried to re-derive the
    /// slope from `y = x^2` (which has lost the sign).
    #[test]
    fn square_gradient_matches_finite_differences() {
        let mut rng = seeded();
        let mut sq = Square::new(4);
        // Both signs are present in the draw, which is the point: a
        // `y`-derived implementation agrees on the positive half and is wrong
        // on the negative one.
        let x = filled(5, 4, &mut rng);
        assert!(
            x.live().iter().any(|v| *v < 0.0) && x.live().iter().any(|v| *v > 0.0),
            "the fixture must straddle zero or it cannot see a sign error"
        );
        let gy = filled(5, 4, &mut rng);
        check_input_grad(&mut sq, x, &gy, "Square");
    }

    /// The sign asymmetry, pinned directly: at `+a` and `-a` the OUTPUT is the
    /// same and the gradient is opposite. This is the whole reason `Square`
    /// reads `x` where `Relu` reads `y`.
    #[test]
    fn square_gradient_sign_follows_the_input_not_the_output() {
        let mut gy = Tensor::new(1, 1);
        gy.row_mut(0)[0] = 1.0;
        let at = |xv: f32| {
            let mut x = Tensor::new(1, 1);
            x.row_mut(0)[0] = xv;
            let mut sq = Square::new(1);
            let mut out = Tensor::default();
            sq.forward(&x, &mut out, true);
            let mut gx = Tensor::default();
            sq.backward(&x, &out, &gy, &mut gx);
            (out.row(0)[0], gx.row(0)[0])
        };
        assert_eq!(at(3.0), (9.0, 6.0));
        assert_eq!(at(-3.0), (9.0, -6.0));
        // The one point with no gradient — a single point, not the half-line
        // that makes a pure ReLU's dead state absorbing.
        assert_eq!(at(0.0), (0.0, 0.0));
    }

    #[test]
    fn batchnorm_gradient_matches_finite_differences() {
        let mut rng = seeded();
        let mut bn = BatchNorm::new(3);
        let x = filled(8, 3, &mut rng);
        let gy = filled(8, 3, &mut rng);
        // BatchNorm's backward accumulates into g_gamma/g_beta and its forward
        // updates running stats; neither affects dL/dx in training mode, which
        // is what this checks.
        check_input_grad(&mut bn, x, &gy, "BatchNorm");
    }

    // ------------------------------------------------------------ loss

    fn loss_only(loss: &MlpLoss, z: &[f32], y: &[f32], w: &[f32]) -> f32 {
        let mut logits = Tensor::new(z.len(), 1);
        for (b, &v) in z.iter().enumerate() {
            logits.row_mut(b)[0] = v;
        }
        let mut grad = Tensor::default();
        loss.loss_and_grad(&logits, y, w, &mut grad)
    }

    fn loss_grad(loss: &MlpLoss, z: &[f32], y: &[f32], w: &[f32]) -> Vec<f32> {
        let mut logits = Tensor::new(z.len(), 1);
        for (b, &v) in z.iter().enumerate() {
            logits.row_mut(b)[0] = v;
        }
        let mut grad = Tensor::default();
        loss.loss_and_grad(&logits, y, w, &mut grad);
        grad.live().to_vec()
    }

    fn check_loss_grad(loss: &MlpLoss, z: &[f32], y: &[f32], w: &[f32], what: &str) {
        let analytic = loss_grad(loss, z, y, w);
        let h = 1e-3;
        for k in 0..z.len() {
            let mut zp = z.to_vec();
            let mut zm = z.to_vec();
            zp[k] += h;
            zm[k] -= h;
            let numeric = (loss_only(loss, &zp, y, w) - loss_only(loss, &zm, y, w)) / (2.0 * h);
            let scale = analytic[k].abs().max(numeric.abs()).max(1e-3);
            assert!(
                (analytic[k] - numeric).abs() / scale < 5e-2,
                "{what}: dL/dz[{k}] analytic {} vs numeric {}",
                analytic[k],
                numeric
            );
        }
    }

    #[test]
    fn bce_gradient_matches_finite_differences() {
        let z = [-2.0f32, -0.3, 0.0, 0.7, 3.1];
        let y = [0.0f32, 1.0, 0.0, 1.0, 1.0];
        let w = [1.0f32, 0.5, 1.0, 0.5, 0.5];
        check_loss_grad(&MlpLoss::Bce, &z, &y, &w, "Bce");
    }

    /// The modulating factor depends on `p`, hence on the logit, so it
    /// contributes its own term to `dL/dz` and the BCE-with-logits gradient does
    /// NOT carry over. Finite differences are what catch a missing or
    /// sign-flipped term — an analytically wrong gradient here still points
    /// downhill and still trains, so no learning test would notice.
    ///
    /// `(0, 0.5)` is included so the degenerate arm is checked by the same
    /// machinery as the rest, not only against `Bce`.
    #[test]
    fn focal_gradient_matches_finite_differences() {
        let z = [-2.0f32, -0.3, 0.0, 0.7, 3.1];
        let y = [0.0f32, 1.0, 0.0, 1.0, 1.0];
        let w = [1.0f32, 0.5, 1.0, 0.5, 0.5];
        for (gamma, alpha) in [
            (0.0f32, 0.5f32),
            (0.0, 0.25),
            (1.0, 0.5),
            (2.0, 0.5),
            (2.0, 0.75),
            (0.5, 0.25),
            (5.0, 0.5),
        ] {
            let loss = MlpLoss::Focal { gamma, alpha };
            check_loss_grad(&loss, &z, &y, &w, &format!("Focal(g={gamma},a={alpha})"));
        }
    }

    /// The cheapest guard against a sign error in the modulating factor's
    /// derivative: at `gamma = 0` with the neutral `alpha` the general case must
    /// reproduce plain BCE, loss AND gradient. Asserted on the BITS, which the
    /// 2x `alpha` normalization is what buys — see the [`MlpLoss::Focal`] docs.
    #[test]
    fn focal_at_zero_gamma_is_exactly_bce() {
        let z = [-2.0f32, -0.3, 0.0, 0.7, 3.1];
        let y = [0.0f32, 1.0, 0.0, 1.0, 1.0];
        let w = [1.0f32, 0.5, 1.0, 0.5, 0.5];
        let focal = MlpLoss::Focal {
            gamma: 0.0,
            alpha: 0.5,
        };

        let (lb, lf) = (
            loss_only(&MlpLoss::Bce, &z, &y, &w),
            loss_only(&focal, &z, &y, &w),
        );
        assert_eq!(lb.to_bits(), lf.to_bits(), "loss {lb} vs {lf}");

        for (a, b) in loss_grad(&MlpLoss::Bce, &z, &y, &w)
            .iter()
            .zip(loss_grad(&focal, &z, &y, &w).iter())
        {
            assert!((a - b).abs() < 1e-7, "grad {a} vs {b}");
        }
    }

    /// The property that replaced the removed `AsymFocal`: `gamma` up-weights
    /// hard rows RELATIVE to easy ones on BOTH sides. Stated as a ratio of
    /// ratios, because `(1 - p_t)^gamma <= 1` shrinks every row in absolute
    /// terms — what focal changes is the balance between them.
    ///
    /// This is also where the sign convention is pinned: an implementation that
    /// modulated by `p_t^gamma` instead of `(1 - p_t)^gamma` (i.e. the discount
    /// the removed variant applied to targets) fails both halves.
    #[test]
    fn focal_up_weights_hard_examples_in_both_classes() {
        let w = [1.0f32];
        let at = |gamma: f32, z: f32, y: f32| {
            loss_only(&MlpLoss::Focal { gamma, alpha: 0.5 }, &[z], &[y], &w)
        };

        // |z| = 2.1972246 is p = 0.9 / 0.1.
        for (y, hard_z, easy_z) in [
            (1.0f32, -2.1972246f32, 2.1972246f32),
            (0.0, 2.1972246, -2.1972246),
        ] {
            let hard_ratio = at(2.0, hard_z, y) / at(0.0, hard_z, y);
            let easy_ratio = at(2.0, easy_z, y) / at(0.0, easy_z, y);
            let class = if y > 0.5 { "target" } else { "decoy" };
            assert!(
                hard_ratio > easy_ratio * 10.0,
                "{class}: the confidently-wrong row must keep far more of its weight \
                 than the easy one: hard {hard_ratio} vs easy {easy_ratio}"
            );
        }
    }

    /// `alpha` tilts the two classes against each other, and the neutral 0.5
    /// really is neutral. Same logit magnitude on both sides, so the only
    /// difference between the two numbers is `alpha_t`.
    #[test]
    fn alpha_tilts_the_class_balance() {
        let w = [1.0f32];
        let at = |alpha: f32, y: f32| {
            loss_only(
                &MlpLoss::Focal { gamma: 2.0, alpha },
                &[if y > 0.5 { -1.0 } else { 1.0 }],
                &[y],
                &w,
            )
        };

        assert_eq!(
            at(0.5, 1.0).to_bits(),
            at(0.5, 0.0).to_bits(),
            "alpha 0.5 must treat the two classes identically"
        );
        assert!(
            at(0.9, 1.0) > at(0.5, 1.0) && at(0.9, 0.0) < at(0.5, 0.0),
            "alpha > 0.5 must favor targets"
        );
        assert!(
            at(0.1, 1.0) < at(0.5, 1.0) && at(0.1, 0.0) > at(0.5, 0.0),
            "alpha < 0.5 must favor decoys"
        );
    }

    /// Extreme logits are where a naive `log(p)` or a `1 - p` complement would
    /// give `inf`, `NaN` or a hard zero. `q` and `p_t` both come straight from
    /// `sigmoid`, so the modulation stays finite and the loss stays ~linear in
    /// `|z|` the way `bce_from_logit` is.
    #[test]
    fn focal_is_finite_at_extreme_logits() {
        let w = [1.0f32];
        for gamma in [0.0f32, 1.0, 2.0, 5.0] {
            for z in [-90.0f32, -40.0, -20.0, 0.0, 20.0, 40.0, 90.0] {
                for y in [0.0f32, 1.0] {
                    let loss = MlpLoss::Focal { gamma, alpha: 0.6 };
                    let l = loss_only(&loss, &[z], &[y], &w);
                    let g = loss_grad(&loss, &[z], &[y], &w)[0];
                    assert!(
                        l.is_finite() && l >= 0.0 && g.is_finite(),
                        "gamma {gamma} z {z} y {y}: loss {l} grad {g}"
                    );
                }
            }
        }
    }

    // ------------------------------------------------------------ transform

    #[test]
    fn transform_culls_dead_and_constant_columns() {
        // col0 informative, col1 all-NaN, col2 constant, col3 informative.
        let ncols = 4;
        let feat: Vec<f64> = vec![
            1.0,
            f64::NAN,
            7.0,
            -1.0, //
            2.0,
            f64::NAN,
            7.0,
            0.5, //
            3.0,
            f64::NAN,
            7.0,
            2.0, //
            4.0,
            f64::NAN,
            7.0,
            9.0, //
        ];
        let rows: Vec<usize> = (0..4).collect();
        let t = ColumnTransform::fit(&feat, ncols, &rows);

        assert_eq!(t.culled(), &[1, 2]);
        assert_eq!(t.width(), 2, "two survivors, neither missable");
        assert_eq!(t.lane_of_input(0), Some(0));
        assert_eq!(t.lane_of_input(1), Some(3));
    }

    #[test]
    fn transform_imputes_missable_and_flags_companion() {
        let ncols = 2;
        // col0 clean, col1 missable (one NaN).
        let feat: Vec<f64> = vec![
            1.0,
            10.0, //
            2.0,
            f64::NAN, //
            3.0,
            30.0, //
            4.0,
            40.0, //
        ];
        let rows: Vec<usize> = (0..4).collect();
        let t = ColumnTransform::fit(&feat, ncols, &rows);

        assert!(t.culled().is_empty());
        // 2 standardized + 1 companion for the missable column.
        assert_eq!(t.width(), 3);

        let mut out = vec![0.0f32; t.width()];
        // Row 1 carries the NaN: imputed to the mean, i.e. 0 standardized, and
        // the companion trips.
        t.apply(&feat[2..4], &mut out);
        assert_eq!(out[1], 0.0, "NaN imputes to the column mean");
        assert_eq!(out[2], 1.0, "companion flags the imputation");

        // Row 0 is finite: companion stays clear.
        t.apply(&feat[0..2], &mut out);
        assert_eq!(out[2], 0.0);
    }

    /// `lane_of_input` / `isna_lane_of_input` must agree with the layout
    /// `apply` actually writes, or a lane-indexed importance vector puts a
    /// companion's weight on the wrong feature. Asserted against `apply`
    /// itself rather than against a hand-copied layout: lane 1 is the only
    /// missable column here, so the slot that trips when lane 1 is NaN IS its
    /// companion by definition.
    #[test]
    fn isna_companions_map_back_to_the_column_they_flag() {
        let ncols = 3;
        // col0 clean, col1 missable, col2 clean.
        let feat: Vec<f64> = vec![
            1.0,
            10.0,
            -1.0, //
            2.0,
            f64::NAN,
            5.0, //
            3.0,
            30.0,
            2.0, //
            4.0,
            40.0,
            9.0, //
        ];
        let rows: Vec<usize> = (0..4).collect();
        let t = ColumnTransform::fit(&feat, ncols, &rows);
        assert_eq!(t.width(), 4, "three survivors + one companion");

        // Standardized block: lanes, no companions.
        assert_eq!(
            (0..3).map(|k| t.lane_of_input(k)).collect::<Vec<_>>(),
            vec![Some(0), Some(1), Some(2)]
        );
        assert_eq!(
            (0..3).map(|k| t.isna_lane_of_input(k)).collect::<Vec<_>>(),
            vec![None; 3]
        );

        // Companion block: the reverse.
        assert_eq!(t.lane_of_input(3), None);
        assert_eq!(t.isna_lane_of_input(3), Some(1));
        // Past the end on both.
        assert_eq!(t.lane_of_input(4), None);
        assert_eq!(t.isna_lane_of_input(4), None);

        // And slot 3 is the one `apply` actually trips for a NaN in lane 1.
        let mut out = vec![0.0f32; t.width()];
        t.apply(&feat[3..6], &mut out);
        assert_eq!(out[3], 1.0);
    }

    /// The property the whole input transform rests on: statistics come from
    /// the rows it was fitted with, never from the rows it will score.
    #[test]
    fn transform_statistics_ignore_held_out_rows() {
        let ncols = 1;
        // Train rows 0..4 are all ~1.0; held-out rows 4..8 are wildly larger.
        let feat: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 1000.0, 2000.0, 3000.0, 4000.0];
        let train: Vec<usize> = (0..4).collect();

        let fitted = ColumnTransform::fit(&feat, ncols, &train);
        let train_only: Vec<f64> = feat[..4].to_vec();
        let reference = ColumnTransform::fit(&train_only, ncols, &train);

        let mut a = vec![0.0f32; fitted.width()];
        let mut b = vec![0.0f32; reference.width()];
        fitted.apply(&feat[0..1], &mut a);
        reference.apply(&feat[0..1], &mut b);
        assert_eq!(a, b, "held-out rows must not move the standardization");
    }

    /// Compiles out in release if written as `debug_assert!`, which is exactly
    /// where it would matter — so it is a returned count, and this test runs
    /// under any profile.
    #[test]
    fn check_clean_catches_nan_in_a_column_fitted_as_clean() {
        let ncols = 1;
        let feat: Vec<f64> = vec![1.0, 2.0, 3.0, f64::NAN];
        let train: Vec<usize> = (0..3).collect();
        let t = ColumnTransform::fit(&feat, ncols, &train);

        assert!(t.check_clean(&feat, &train).is_empty());
        assert_eq!(
            t.check_clean(&feat, &[3]),
            vec![0],
            "a NaN in a column fitted as clean must be reported"
        );
    }

    // ------------------------------------------------------------ learning

    /// Build the noisy-XOR fixture, train `cfg` on it from init seed `seed`, and
    /// return `(accuracy, final_train_loss)`.
    ///
    /// ONE `StdRng`, consumed in a fixed order — 512 rows of data, then the
    /// init, then the per-epoch shuffles. That order is load-bearing: it is what
    /// makes a given `(cfg, seed)` a fixed experiment, so every caller here is
    /// comparing the same trajectory.
    fn xor_accuracy(cfg: &MlpConfig, seed: u64) -> (f32, f32) {
        let mut rng = StdRng::seed_from_u64(seed);
        let n = 512;
        let mut x = Tensor::new(n, 2);
        let mut y = vec![0.0f32; n];
        for i in 0..n {
            let a = if i % 2 == 0 { 1.0 } else { -1.0 };
            let b = if (i / 2) % 2 == 0 { 1.0 } else { -1.0 };
            let r = x.row_mut(i);
            r[0] = a + 0.1 * normal(&mut rng);
            r[1] = b + 0.1 * normal(&mut rng);
            y[i] = if a * b > 0.0 { 1.0 } else { 0.0 };
        }
        let w = vec![1.0f32; n];

        let mut model = Mlp::feedforward(2, &cfg.hidden, cfg.activation, &mut rng);
        let mut opt = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
        let loss = model.train(cfg, &x, &y, &w, &mut opt, &mut rng);

        let out = model.forward(&x, false);
        let correct = (0..n)
            .filter(|&i| (out.row(i)[0] > 0.0) == (y[i] > 0.5))
            .count();
        (correct as f32 / n as f32, loss)
    }

    /// The XOR config the learning tests share: everything but the activation is
    /// held fixed, `lr` and `weight_decay` INCLUDED FROM THE DEFAULT — see
    /// [`mlp_learns_xor_which_no_linear_model_can`] for why the inheritance is
    /// deliberate.
    fn xor_cfg(activation: Activation) -> MlpConfig {
        MlpConfig {
            hidden: vec![16, 8],
            epochs: 400,
            batch_size: 64,
            activation,
            ..MlpConfig::default()
        }
    }

    /// Without this, every other test still passes if the MLP silently
    /// collapses to a linear map — XOR is the cheapest thing a linear
    /// discriminant provably cannot separate.
    ///
    /// Run over several seeds on purpose. A single seed passing proves very
    /// little here: full-batch descent on this net was observed settling into a
    /// dead-ReLU plateau at exactly 0.75 accuracy / `ln(2)/2` loss for one
    /// specific init, and a one-seed test would have called that either a pass
    /// or a bug depending on luck.
    ///
    /// `lr` is INHERITED from [`MlpConfig::default`] rather than pinned, so the
    /// test trains at whatever step size the crate ships (3e-4 until the 2026-07
    /// retune, 1.2e-3 since, green across all eight seeds in both regimes).
    /// Keep the inheritance: pinning `lr` would leave nothing exercising the
    /// shipped value on a nonlinear problem at all.
    ///
    /// IT IS NOT A CANARY ON THE STEP SIZE, and an earlier version of this note
    /// claiming otherwise was measured false. All eight seeds pass at `lr` =
    /// 3e-4, 8e-4, 1.2e-3, 2e-3, 3e-3, 6e-3 and 1.2e-2 — insensitive across a 40x
    /// range, including the 1e-2 an older dead-plateau warning here named as the
    /// danger. A test that passes everywhere cannot detect a cliff, so a green run
    /// says nothing about whether the default is well chosen; that question is
    /// settled by the sweep in [`MlpConfig::compiled_default`], on real data.
    /// What this test does still catch is the collapse it was written for: a net
    /// degenerated to a linear map caps near 0.75 accuracy on XOR, well under the
    /// 0.95 asserted below.
    #[test]
    fn mlp_learns_xor_which_no_linear_model_can() {
        // Seed 7 is here on purpose: it is the init that trapped at exactly
        // 0.75 / `ln(2)/2` back when biases initialized to zero.
        // Verified green across 16 seeds; trimmed to 8 to keep the test under
        // ~10s. Seeds 7 and 13 are pinned members, not arbitrary.
        for seed in [7u64, 13, 42, 99, 1234, 20260728, 3, 31337] {
            let (acc, loss) = xor_accuracy(&xor_cfg(Activation::LeakyRelu), seed);
            assert!(
                acc > 0.95,
                "seed {seed}: XOR accuracy {acc} (loss {loss}) — collapsed to linear?"
            );
        }
    }

    /// [`Activation::Square`] can fit something a linear model cannot, at the
    /// SHIPPED `lr` — the same fixture, seed set and threshold as
    /// [`mlp_learns_xor_which_no_linear_model_can`], with only the activation
    /// changed, so the two are directly comparable.
    ///
    /// Swept over seeds for the reason every learning test here is: this module
    /// has two documented init-dependent training traps, and `x^2` adds a third
    /// hazard of its own — it is unbounded, so a bad init can diverge outright
    /// rather than merely converge slowly (see [`Square`]). A single seed would
    /// not distinguish "works" from "worked once".
    ///
    /// MEASURED, and the point of recording it: all eight seeds reach the same
    /// threshold as the rectifier at the shipped `lr` of 1.2e-3, and the loss it
    /// gets there is LOWER on every seed. That says the activation is usable on a
    /// small nonlinear problem; it says nothing about the 101-column rescore
    /// lanes, where the divergence hazard scales with the input dimension and
    /// nothing has been measured at all.
    #[test]
    fn square_activation_learns_xor() {
        for seed in [7u64, 13, 42, 99, 1234, 20260728, 3, 31337] {
            let cfg = xor_cfg(Activation::Square);
            let (acc, loss) = xor_accuracy(&cfg, seed);
            assert!(
                loss.is_finite(),
                "seed {seed}: the squaring net diverged (loss {loss}) at lr {}",
                cfg.lr
            );
            assert!(
                acc > 0.95,
                "seed {seed}: XOR accuracy {acc} (loss {loss}) at lr {} — \
                 report the lr it needs rather than retuning the default",
                cfg.lr
            );
        }
    }

    /// The only test that the MLP's per-row weights do anything — this is what
    /// carries the GBM's 0.5/1.0 target/decoy balance into the MLP.
    ///
    /// Swept over INIT seeds, like every other learning test here. It is a PAIRED
    /// comparison with matched init (both arms start from the same
    /// `StdRng::seed_from_u64(seed)`, so the two nets are bit-identical before the
    /// first step and the only difference is the weight vector), which makes it far
    /// more robust than an absolute-threshold test — but "robust" was an argument,
    /// not a measurement, and this module has two documented init-dependent
    /// training traps. The sweep is what turns it into a measurement.
    ///
    /// The DATA seed stays fixed at 11: varying the draws too would change what is
    /// being asserted from "weights move the boundary" to "weights move the
    /// boundary on any sample", which needs a tolerance rather than a strict
    /// inequality.
    /// THE `patience: None` contract: with early stopping off, handing
    /// `train_reporting` a held-out set changes the LOG and nothing else. That
    /// is the behavior the whole module had before early stopping existed, and
    /// `None` has to reproduce it exactly.
    ///
    /// This is a leak assertion, not a tidiness one. The set passed here is by
    /// construction rows the fitted model will later be asked to score, so if
    /// reporting could nudge a weight — through the optimizer, the RNG stream, or
    /// a stale gradient buffer left by `eval_loss` — the model would have been
    /// fitted on rows it is scoring, and every downstream leak test would still
    /// pass because they all check the PARTITION rather than the training loop.
    /// Asserted on raw bits over the whole score vector, at several seeds.
    ///
    /// `patience: None` is EXPLICIT here rather than inherited: the default is
    /// `Some(5)`, under which passing a held-out set is SUPPOSED to be able to
    /// change the fit. The second half of this test asserts that it can, so a
    /// `None` that silently early-stopped anyway — or a `Some` that silently did
    /// not — fails one arm or the other.
    #[test]
    fn patience_none_makes_a_held_out_set_reporting_only() {
        for seed in [7u64, 13, 42] {
            let mut dr = StdRng::seed_from_u64(seed);
            let n = 128;
            let mut x = Tensor::new(n, 3);
            let mut y = vec![0.0f32; n];
            for i in 0..n {
                let pos = i % 2 == 0;
                let r = x.row_mut(i);
                r[0] = if pos { 0.7 } else { -0.7 } + 0.5 * normal(&mut dr);
                r[1] = normal(&mut dr);
                r[2] = if pos { -0.3 } else { 0.3 } + 0.5 * normal(&mut dr);
                y[i] = if pos { 1.0 } else { 0.0 };
            }
            let w = vec![1.0f32; n];
            // A DIFFERENT set, standing in for the early-stop fold. Pure noise
            // against alternating labels, so its loss cannot fall for long and
            // the `Some` arm below is guaranteed something to stop on.
            let mut vx = Tensor::new(32, 3);
            for i in 0..32 {
                for j in 0..3 {
                    vx.row_mut(i)[j] = normal(&mut dr);
                }
            }
            let vy: Vec<f32> = (0..32)
                .map(|i| if i % 2 == 0 { 1.0 } else { 0.0 })
                .collect();
            let vw = vec![1.0f32; 32];

            let cfg = MlpConfig {
                hidden: vec![8],
                // Long enough that the noise validation set has certainly
                // turned by the end, which is what the control arm needs.
                epochs: 400,
                batch_size: 32,
                early_stopping_patience: None,
                ..MlpConfig::default()
            };

            let run = |cfg: &MlpConfig, val: Option<ValSet<'_>>| {
                let mut r = StdRng::seed_from_u64(99);
                let mut m = Mlp::feedforward(3, &cfg.hidden, cfg.activation, &mut r);
                let mut o = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
                let out = m.train_reporting(cfg, &x, &y, &w, &mut o, &mut r, val, "t");
                let p = m.forward(&x, false);
                let bits = (0..n).map(|i| p.row(i)[0].to_bits()).collect::<Vec<u32>>();
                (out, bits)
            };
            let vset = || {
                Some(ValSet {
                    x: &vx,
                    y: &vy,
                    w: &vw,
                })
            };

            let (out_without, without) = run(&cfg, None);
            let (out_with, with) = run(&cfg, vset());
            assert_eq!(
                without, with,
                "seed {seed}: with patience None a held-out set moved the fit — it must only log"
            );
            // Non-vacuity: the fit is not a constant that any two runs would match.
            assert!(
                without
                    .iter()
                    .collect::<std::collections::HashSet<_>>()
                    .len()
                    > 1,
                "seed {seed}: scores are all identical, so the comparison proves nothing"
            );
            // ...and `None` really means the whole budget, unsnapshotted.
            for o in [out_without, out_with] {
                assert_eq!(
                    o.epochs_run, cfg.epochs,
                    "seed {seed}: None must run the budget"
                );
                assert_eq!(
                    o.best_epoch, None,
                    "seed {seed}: None must not track a best"
                );
                assert!(!o.restored, "seed {seed}: None must not roll back");
            }

            // THE CONTROL: the same call with patience ON does change the fit,
            // so the equality above is a property of `None` and not of a
            // held-out set that could never matter.
            let es = MlpConfig {
                early_stopping_patience: Some(5),
                ..cfg.clone()
            };
            let (out_es, es_bits) = run(&es, vset());
            assert!(
                out_es.epochs_run < es.epochs,
                "seed {seed}: the fixture must overfit its noise validation set for the \
                 control to mean anything (ran {} of {} epochs)",
                out_es.epochs_run,
                es.epochs
            );
            assert!(
                es_bits != without,
                "seed {seed}: early stopping produced the SAME weights as the full budget, \
                 so `patience: None` reproducing them proves nothing"
            );
        }
    }

    // ------------------------------------------------------- early stopping

    /// `(train x, y, w, validation x, y, w)` — see [`overfitting_fixture`].
    type OverfitFixture = (Tensor, Vec<f32>, Vec<f32>, Tensor, Vec<f32>, Vec<f32>);

    /// A fixture that MUST overfit: 48 training rows, 24 columns of which
    /// exactly one carries a weak signal, against a 256-row validation set from
    /// the same distribution. There is far more capacity than signal, so the
    /// train loss walks to zero while the validation loss bottoms out early and
    /// then climbs — which is the situation early stopping exists for.
    fn overfitting_fixture(seed: u64) -> OverfitFixture {
        const NCOLS: usize = 24;
        let mut rng = StdRng::seed_from_u64(seed);
        let mut make = |n: usize| {
            let mut t = Tensor::new(n, NCOLS);
            let mut labels = vec![0.0f32; n];
            for i in 0..n {
                let pos = i % 2 == 0;
                let r = t.row_mut(i);
                // Column 0 is weakly informative; 1..24 are pure noise the net
                // can memorize the training rows with and nothing else.
                r[0] = if pos { 0.6 } else { -0.6 } + normal(&mut rng);
                for j in 1..NCOLS {
                    r[j] = normal(&mut rng);
                }
                labels[i] = if pos { 1.0 } else { 0.0 };
            }
            let w = vec![1.0f32; n];
            (t, labels, w)
        };
        let (x, y, w) = make(48);
        let (vx, vy, vw) = make(256);
        (x, y, w, vx, vy, vw)
    }

    fn overfitting_cfg(patience: Option<usize>) -> MlpConfig {
        MlpConfig {
            hidden: vec![32, 16],
            lr: 1e-2,
            weight_decay: 0.0,
            epochs: 200,
            batch_size: 16,
            early_stopping_patience: patience,
            ..MlpConfig::default()
        }
    }

    /// Train the overfitting fixture and hand back the outcome plus the
    /// validation loss OF THE WEIGHTS THE CALLER IS LEFT WITH.
    fn run_overfit(
        cfg: &MlpConfig,
        init_seed: u64,
        fx: &OverfitFixture,
        with_val: bool,
    ) -> (TrainOutcome, f32, Vec<u32>) {
        let (x, y, w, vx, vy, vw) = fx;
        let mut r = StdRng::seed_from_u64(init_seed);
        let mut m = Mlp::feedforward(x.cols, &cfg.hidden, cfg.activation, &mut r);
        let mut o = Adam::new(cfg.lr).with_weight_decay(cfg.weight_decay);
        let val = with_val.then(|| ValSet {
            x: vx,
            y: vy,
            w: vw,
        });
        let outcome = m.train_reporting(cfg, x, y, w, &mut o, &mut r, val, "es");
        let vl = m.eval_loss(vx, vy, vw, &cfg.loss);
        let bits = m
            .params_and_grads()
            .into_iter()
            .flat_map(|(p, _, _)| p.iter().map(|v| v.to_bits()).collect::<Vec<_>>())
            .collect();
        (outcome, vl, bits)
    }

    /// Early stopping stops BEFORE the budget on a fixture that overfits, and
    /// the weights it leaves behind are no worse on the held-out set than the
    /// ones the full budget produces.
    ///
    /// Swept over seeds, like every learning test in this module: this one is
    /// doubly init-dependent, since both the training trap and the epoch the
    /// validation curve turns at move with the initialization.
    #[test]
    fn early_stopping_stops_early_and_does_not_hurt_the_held_out_loss() {
        for seed in [7u64, 13, 42, 1234, 31337] {
            let fx = overfitting_fixture(seed ^ 0xF00D);
            let es = overfitting_cfg(Some(5));
            let full = overfitting_cfg(None);

            let (out, es_val, _) = run_overfit(&es, seed, &fx, true);
            let (full_out, full_val, _) = run_overfit(&full, seed, &fx, true);

            assert_eq!(full_out.epochs_run, full.epochs);
            assert!(
                out.epochs_run < es.epochs,
                "seed {seed}: ran the whole {} epoch budget, so nothing was stopped",
                es.epochs
            );
            assert!(
                out.restored,
                "seed {seed}: a run that stopped early must roll back"
            );
            // The stop is exactly `patience` epochs after the best one.
            assert_eq!(
                out.epochs_run,
                out.best_epoch.unwrap() + 1 + 5,
                "seed {seed}: stopped somewhere other than best + patience"
            );
            assert!(
                es_val <= full_val,
                "seed {seed}: early stopping made the held-out loss WORSE \
                 ({es_val} > {full_val}) — it is supposed to be the point"
            );
            // Non-vacuity: the full budget really does overfit here, so the
            // comparison above is not two identical numbers.
            assert!(
                full_val > es_val,
                "seed {seed}: the full budget did not overfit this fixture \
                 (held-out {full_val} vs {es_val}), so there was nothing to stop"
            );
        }
    }

    /// THE restore assertion: the weights left behind are the BEST-seen ones,
    /// not the ones training was standing on when patience ran out.
    ///
    /// The comparison arm is the same trajectory truncated at the same epoch
    /// WITHOUT a rollback — `patience: None` with `epochs` set to what the early
    /// stopped run actually executed. Nothing else differs: the RNG feeds only
    /// the init and the per-epoch shuffle, and `eval_loss` draws from it not at
    /// all, so the two runs are bit-identical up to the point one of them rolls
    /// back. An implementation that stopped without restoring would produce that
    /// arm's weights and fail both assertions below.
    #[test]
    fn early_stopping_restores_the_best_weights_not_the_last_ones() {
        for seed in [7u64, 13, 42, 1234, 31337] {
            let fx = overfitting_fixture(seed ^ 0xF00D);
            let es = overfitting_cfg(Some(5));
            let (out, es_val, es_bits) = run_overfit(&es, seed, &fx, true);

            // The held-out loss of the restored weights IS the best one seen,
            // to the bit — not merely close to it.
            assert_eq!(
                es_val.to_bits(),
                out.best_val_loss.unwrap().to_bits(),
                "seed {seed}: the restored weights do not reproduce the best held-out \
                 loss ({es_val} vs {})",
                out.best_val_loss.unwrap()
            );

            let stopped_here = MlpConfig {
                epochs: out.epochs_run,
                ..overfitting_cfg(None)
            };
            let (no_restore_out, no_restore_val, no_restore_bits) =
                run_overfit(&stopped_here, seed, &fx, true);
            assert_eq!(no_restore_out.epochs_run, out.epochs_run);
            assert!(
                es_bits != no_restore_bits,
                "seed {seed}: restoring changed no weight, so this test cannot tell a \
                 restoring implementation from a non-restoring one"
            );
            assert!(
                es_val < no_restore_val,
                "seed {seed}: stopping WITHOUT restoring would have left held-out loss \
                 {no_restore_val}; the restore must beat it, got {es_val}"
            );
        }
    }

    /// TIES ARE NOT IMPROVEMENTS, pinned on the one fixture where every epoch
    /// ties exactly: `lr = 0` freezes every parameter (both the Adam step and
    /// the decoupled decay scale by `lr`), so the held-out loss is bit-identical
    /// at every epoch.
    ///
    /// Under the strict `<` rule epoch 0 is the only improvement, the patience
    /// counter then runs uninterrupted, and training stops at epoch `patience`.
    /// Under a `<=` rule every epoch would improve, the counter would never
    /// advance, and the run would train the full budget and keep the LAST epoch
    /// — so this test distinguishes the two rules exactly.
    #[test]
    fn a_tie_is_not_an_improvement() {
        for patience in [1usize, 3, 5] {
            let fx = overfitting_fixture(11);
            let cfg = MlpConfig {
                lr: 0.0,
                weight_decay: 0.0,
                epochs: 100,
                early_stopping_patience: Some(patience),
                ..overfitting_cfg(Some(patience))
            };
            let (out, _, _) = run_overfit(&cfg, 5, &fx, true);
            assert_eq!(
                out.best_epoch,
                Some(0),
                "patience {patience}: with every epoch tied, the FIRST must be kept"
            );
            assert_eq!(
                out.epochs_run,
                patience + 1,
                "patience {patience}: a dead-flat curve must still run the counter down"
            );
        }
    }

    /// Determinism, including the stopping decision: two runs of the same build
    /// on the same input agree on the epoch they stopped at, the epoch they
    /// kept, and every parameter bit.
    #[test]
    fn early_stopping_is_bit_reproducible() {
        for seed in [7u64, 42, 31337] {
            let fx = overfitting_fixture(seed ^ 0xF00D);
            let cfg = overfitting_cfg(Some(5));
            let (a_out, a_val, a_bits) = run_overfit(&cfg, seed, &fx, true);
            let (b_out, b_val, b_bits) = run_overfit(&cfg, seed, &fx, true);
            assert_eq!(a_out, b_out, "seed {seed}: the stopping decision moved");
            assert_eq!(a_val.to_bits(), b_val.to_bits(), "seed {seed}");
            assert_eq!(a_bits, b_bits, "seed {seed}: the restored weights moved");
        }
    }

    /// With no held-out set the knob is INERT — there is nothing to measure a
    /// stopping decision against, so the full budget runs and `Mlp::train`
    /// behaves as it always did even on the default config.
    #[test]
    fn early_stopping_without_a_validation_set_is_a_no_op() {
        let fx = overfitting_fixture(11);
        let cfg = overfitting_cfg(Some(5));
        let (out, _, with_knob) = run_overfit(&cfg, 7, &fx, false);
        assert_eq!(out.epochs_run, cfg.epochs);
        assert_eq!(out.best_epoch, None);
        assert!(!out.restored);

        let (_, _, without_knob) = run_overfit(&overfitting_cfg(None), 7, &fx, false);
        assert_eq!(
            with_knob, without_knob,
            "with no validation set the patience knob must change nothing"
        );
    }

    #[test]
    fn eval_loss_leaves_the_parameters_untouched() {
        let mut rng = StdRng::seed_from_u64(5);
        let mut m = Mlp::feedforward(4, &[6], Activation::default(), &mut rng);
        let x = filled(16, 4, &mut rng);
        let y: Vec<f32> = (0..16).map(|i| (i % 2) as f32).collect();
        let w = vec![1.0f32; 16];

        let snapshot = |m: &mut Mlp| -> Vec<u32> {
            m.params_and_grads()
                .into_iter()
                .flat_map(|(p, _, _)| p.iter().map(|v| v.to_bits()).collect::<Vec<_>>())
                .collect()
        };

        let before = snapshot(&mut m);
        let l = m.eval_loss(&x, &y, &w, &MlpLoss::Bce);
        let after = snapshot(&mut m);
        assert_eq!(before, after, "eval_loss mutated a parameter");
        assert!(l.is_finite() && l > 0.0, "eval_loss returned {l}");
    }

    #[test]
    fn sample_weights_shift_the_decision_boundary() {
        // Two overlapping classes; up-weighting one must move the fitted
        // boundary toward the other.
        let mut rng = StdRng::seed_from_u64(11);
        let n = 400;
        let mut x = Tensor::new(n, 1);
        let mut y = vec![0.0f32; n];
        for i in 0..n {
            let pos = i % 2 == 0;
            x.row_mut(i)[0] = if pos { 0.5 } else { -0.5 } + 0.8 * normal(&mut rng);
            y[i] = if pos { 1.0 } else { 0.0 };
        }

        let mean_logit = |w: &[f32], rng: &mut StdRng| {
            let mut m = Mlp::feedforward(1, &[8], Activation::default(), rng);
            let mut o = Adam::new(1e-2);
            for _ in 0..200 {
                m.train_step(&x, &y, w, &MlpLoss::Bce, &mut o);
            }
            let out = m.forward(&x, false);
            (0..n).map(|i| out.row(i)[0]).sum::<f32>() / n as f32
        };

        let balanced = vec![1.0f32; n];
        let decoy_heavy: Vec<f32> = (0..n).map(|i| if i % 2 == 0 { 0.5 } else { 1.0 }).collect();

        for seed in [3u64, 7, 42] {
            // MATCHED INIT: same seed on both sides, so the two runs differ in
            // nothing but `w`.
            let mut r1 = StdRng::seed_from_u64(seed);
            let mut r2 = StdRng::seed_from_u64(seed);
            let m_bal = mean_logit(&balanced, &mut r1);
            let m_dec = mean_logit(&decoy_heavy, &mut r2);
            assert!(
                m_dec < m_bal,
                "seed {seed}: up-weighting decoys must pull logits down: {m_dec} !< {m_bal}"
            );
        }
    }

    // ------------------------------------------------- dev env overrides
    //
    // These drive `MlpConfig::apply_overrides` with a fake lookup instead of
    // `std::env::set_var`. That is not a convenience: `set_var` is `unsafe` in
    // edition 2024 and process-global, so a test that set `TIMSSEEK_MLP_LR`
    // would be racing every other test in this binary. The lookup argument makes
    // the whole override path a pure function.

    /// Lookup over a fixed `(name, value)` table; everything else is unset.
    fn env_of<'a>(pairs: &'a [(&'a str, &'a str)]) -> impl Fn(&str) -> Option<String> + use<'a> {
        move |k: &str| {
            pairs
                .iter()
                .find(|(n, _)| *n == k)
                .map(|(_, v)| (*v).to_string())
        }
    }

    /// THE LOAD-BEARING PROPERTY: with nothing set, the config is bit-identical
    /// to the compiled default, so the override cannot perturb a production run
    /// that did not ask for it.
    ///
    /// The literal values are re-asserted here on purpose. This affordance adds
    /// a way to override the defaults and must not retune them; a diff that
    /// changes one has to change this test and say so.
    #[test]
    fn unset_env_changes_nothing() {
        let mut cfg = MlpConfig::compiled_default();
        let applied = cfg.apply_overrides(env_of(&[]));
        assert!(
            applied.is_empty(),
            "nothing was set, so nothing may be reported as applied: {applied:?}"
        );
        assert_eq!(cfg, MlpConfig::compiled_default());

        // Retuned 2026-07 from the unmeasured `[64, 32]` / 3e-4 / 256 / 30 to
        // the winning arm of the two-scale sweep documented on
        // `compiled_default`. `lr` and `batch_size` are COUPLED there: a diff
        // that touches one of these two literals without the other is a bug.
        assert_eq!(cfg.hidden, vec![32, 16]);
        assert_eq!(cfg.lr.to_bits(), 1.2e-3f32.to_bits());
        assert_eq!(cfg.weight_decay.to_bits(), 1e-4f32.to_bits());
        assert_eq!(cfg.epochs, 60);
        assert_eq!(cfg.batch_size, 1024);
        // `TIMSSEEK_MLP_LOSS` can SELECT the focal loss but must never make it
        // the default: every number on `compiled_default` was measured on plain
        // BCE, and `Focal` is an unswept experiment. Unit variant, so this
        // equality is as exact as the `to_bits` pinning above.
        assert_eq!(cfg.loss, MlpLoss::Bce);
        assert_eq!(cfg.seed, 0x2545_F491_4F6C_DD1D);
        assert_eq!(cfg.early_stopping_patience, Some(5));
        // The rectifier every number on `compiled_default` was measured with.
        // `Square` is an experiment and must never become the default by
        // accident — see [`Square`] for the divergence hazard.
        assert_eq!(cfg.activation, Activation::LeakyRelu);
    }

    /// An unrelated variable is not a read, and one variable moves ONE field.
    #[test]
    fn one_variable_moves_one_field() {
        let base = MlpConfig::compiled_default();
        let mut cfg = base.clone();
        let applied = cfg.apply_overrides(env_of(&[
            ("TIMSSEEK_MLP_LR", "1e-3"),
            ("TIMSSEEK_MLP_NOT_A_KNOB", "7"),
        ]));
        assert_eq!(applied, vec!["TIMSSEEK_MLP_LR=\"1e-3\""]);
        assert_eq!(cfg.lr.to_bits(), 1e-3f32.to_bits());
        assert_eq!(
            MlpConfig { lr: base.lr, ..cfg },
            base,
            "overriding `lr` moved some other field too"
        );
    }

    /// Every documented variable is wired to its field, and the audit trail
    /// lists them all. A knob that parses but is never assigned would otherwise
    /// produce exactly the silent no-op sweep this whole affordance is meant to
    /// prevent.
    #[test]
    fn every_documented_variable_reaches_its_field() {
        let mut cfg = MlpConfig::compiled_default();
        let applied = cfg.apply_overrides(env_of(&[
            ("TIMSSEEK_MLP_EPOCHS", "77"),
            ("TIMSSEEK_MLP_LR", "2.5e-3"),
            ("TIMSSEEK_MLP_BATCH_SIZE", "512"),
            ("TIMSSEEK_MLP_HIDDEN", "128,64,16"),
            ("TIMSSEEK_MLP_WEIGHT_DECAY", "0"),
            ("TIMSSEEK_MLP_PATIENCE", "9"),
            ("TIMSSEEK_MLP_ACTIVATION", "square"),
            ("TIMSSEEK_MLP_SEED", "12345"),
            ("TIMSSEEK_MLP_LOSS", "focal:2:0.25"),
        ]));
        assert_eq!(cfg.epochs, 77);
        assert_eq!(cfg.lr.to_bits(), 2.5e-3f32.to_bits());
        assert_eq!(cfg.batch_size, 512);
        assert_eq!(cfg.hidden, vec![128, 64, 16]);
        assert_eq!(cfg.weight_decay, 0.0);
        assert_eq!(cfg.early_stopping_patience, Some(9));
        assert_eq!(cfg.activation, Activation::Square);
        assert_eq!(cfg.seed, 12345);
        assert_eq!(
            cfg.loss,
            MlpLoss::Focal {
                gamma: 2.0,
                alpha: 0.25
            }
        );
        assert_eq!(applied.len(), 9, "every read must be reported: {applied:?}");
    }

    /// The explicit spellings for the two fields that have an "off": `hidden`
    /// can be empty (a bare linear model) and `patience` is an `Option`.
    #[test]
    fn the_off_spellings_parse() {
        let mut cfg = MlpConfig::compiled_default();
        cfg.apply_overrides(env_of(&[
            ("TIMSSEEK_MLP_HIDDEN", "none"),
            ("TIMSSEEK_MLP_PATIENCE", "off"),
        ]));
        assert!(cfg.hidden.is_empty());
        assert_eq!(cfg.early_stopping_patience, None);

        let mut cfg = MlpConfig::compiled_default();
        cfg.apply_overrides(env_of(&[
            ("TIMSSEEK_MLP_HIDDEN", "NONE"),
            ("TIMSSEEK_MLP_PATIENCE", "None"),
        ]));
        assert!(cfg.hidden.is_empty(), "the spelling is case-insensitive");
        assert_eq!(cfg.early_stopping_patience, None);

        // `hidden = []` must still build a usable net, or the spelling is a trap.
        let mut rng = seeded();
        let mut net = Mlp::feedforward(4, &cfg.hidden, cfg.activation, &mut rng);
        let x = filled(3, 4, &mut rng);
        assert_eq!(net.forward(&x, false).rows(), 3);
    }

    /// Whitespace around a value survives a shell that added it.
    #[test]
    fn values_are_trimmed() {
        let mut cfg = MlpConfig::compiled_default();
        cfg.apply_overrides(env_of(&[
            ("TIMSSEEK_MLP_EPOCHS", " 40 "),
            ("TIMSSEEK_MLP_HIDDEN", " 32 , 8 "),
        ]));
        assert_eq!(cfg.epochs, 40);
        assert_eq!(cfg.hidden, vec![32, 8]);
    }

    /// Hex, so the compiled default seed can be restated verbatim.
    #[test]
    fn seed_takes_decimal_or_hex() {
        for spelling in [
            "0x2545F4914F6CDD1D",
            "0X2545f4914f6cdd1d",
            "2685821657736338717",
        ] {
            let mut cfg = MlpConfig::compiled_default();
            cfg.seed = 0;
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_SEED", spelling)]));
            assert_eq!(
                cfg.seed, 0x2545_F491_4F6C_DD1D,
                "{spelling} must round-trip to the default seed"
            );
        }
    }

    /// A malformed value must ABORT, not fall back to the default. Each of these
    /// would otherwise produce a sweep row labelled with a value that never
    /// trained. One case per rejection rule.
    macro_rules! rejects {
        ($name:ident, $var:literal, $val:literal) => {
            #[test]
            #[should_panic(expected = "is not usable")]
            fn $name() {
                let mut cfg = MlpConfig::compiled_default();
                cfg.apply_overrides(env_of(&[($var, $val)]));
            }
        };
    }

    rejects!(rejects_non_numeric_lr, "TIMSSEEK_MLP_LR", "banana");
    rejects!(rejects_zero_lr, "TIMSSEEK_MLP_LR", "0");
    rejects!(rejects_negative_lr, "TIMSSEEK_MLP_LR", "-1e-3");
    rejects!(rejects_nan_lr, "TIMSSEEK_MLP_LR", "NaN");
    rejects!(rejects_zero_epochs, "TIMSSEEK_MLP_EPOCHS", "0");
    rejects!(rejects_float_epochs, "TIMSSEEK_MLP_EPOCHS", "30.5");
    rejects!(rejects_zero_batch_size, "TIMSSEEK_MLP_BATCH_SIZE", "0");
    rejects!(
        rejects_negative_weight_decay,
        "TIMSSEEK_MLP_WEIGHT_DECAY",
        "-1"
    );
    rejects!(rejects_empty_hidden, "TIMSSEEK_MLP_HIDDEN", "");
    rejects!(rejects_zero_width_hidden, "TIMSSEEK_MLP_HIDDEN", "64,0");
    rejects!(rejects_non_numeric_hidden, "TIMSSEEK_MLP_HIDDEN", "64x32");
    rejects!(
        rejects_trailing_comma_hidden,
        "TIMSSEEK_MLP_HIDDEN",
        "64,32,"
    );
    // 0 is neither a patience nor a synonym for `none` — see `parse_patience`.
    rejects!(rejects_zero_patience, "TIMSSEEK_MLP_PATIENCE", "0");
    rejects!(rejects_negative_patience, "TIMSSEEK_MLP_PATIENCE", "-1");
    rejects!(rejects_non_numeric_seed, "TIMSSEEK_MLP_SEED", "0xZZ");
    // A near-miss spelling is a rejection, not a silent fall back to the
    // rectifier: a sweep arm labelled `square` that trained a ReLU is exactly
    // the outcome this abort exists to prevent.
    rejects!(
        rejects_unknown_activation,
        "TIMSSEEK_MLP_ACTIVATION",
        "relu"
    );
    rejects!(rejects_empty_activation, "TIMSSEEK_MLP_ACTIVATION", "");
    rejects!(rejects_unknown_loss, "TIMSSEEK_MLP_LOSS", "hinge");
    rejects!(rejects_empty_loss, "TIMSSEEK_MLP_LOSS", "");
    // Both focal parameters are spelled or neither is: an invented `gamma` is
    // exactly the mislabelled sweep row the abort exists for.
    rejects!(rejects_bare_focal, "TIMSSEEK_MLP_LOSS", "focal");
    rejects!(rejects_focal_without_alpha, "TIMSSEEK_MLP_LOSS", "focal:2");
    rejects!(
        rejects_focal_with_extra_field,
        "TIMSSEEK_MLP_LOSS",
        "focal:2:0.25:9"
    );
    rejects!(
        rejects_non_numeric_focal_gamma,
        "TIMSSEEK_MLP_LOSS",
        "focal:two:0.25"
    );
    rejects!(
        rejects_negative_focal_gamma,
        "TIMSSEEK_MLP_LOSS",
        "focal:-1:0.25"
    );
    rejects!(
        rejects_non_finite_focal_gamma,
        "TIMSSEEK_MLP_LOSS",
        "focal:inf:0.25"
    );
    // `alpha` is applied as `2*alpha` / `2*(1-alpha)`, so an endpoint zeroes one
    // class out entirely — a rejection, not a boundary case.
    rejects!(rejects_zero_focal_alpha, "TIMSSEEK_MLP_LOSS", "focal:2:0");
    rejects!(rejects_one_focal_alpha, "TIMSSEEK_MLP_LOSS", "focal:2:1");
    rejects!(
        rejects_negative_focal_alpha,
        "TIMSSEEK_MLP_LOSS",
        "focal:2:-0.25"
    );
    rejects!(
        rejects_above_one_focal_alpha,
        "TIMSSEEK_MLP_LOSS",
        "focal:2:1.5"
    );
    rejects!(
        rejects_non_finite_focal_alpha,
        "TIMSSEEK_MLP_LOSS",
        "focal:2:nan"
    );

    /// The accepted spellings of the loss. `focal:0:0.5` is the neutral point —
    /// see [`MlpLoss::Focal`] for why it is bit-identical to `Bce` rather than
    /// merely proportional to it, which is what makes an `alpha` sweep's rows
    /// comparable.
    #[test]
    fn loss_spellings_parse() {
        for spelling in ["bce", "BCE", " Bce "] {
            let mut cfg = MlpConfig::compiled_default();
            cfg.loss = MlpLoss::Focal {
                gamma: 3.0,
                alpha: 0.9,
            };
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_LOSS", spelling)]));
            assert_eq!(cfg.loss, MlpLoss::Bce, "{spelling:?}");
        }

        for (spelling, gamma, alpha) in [
            ("focal:2:0.25", 2.0f32, 0.25f32),
            ("FOCAL:2:0.25", 2.0, 0.25),
            (" focal : 2 : 0.25 ", 2.0, 0.25),
            ("focal:2.0:2.5e-1", 2.0, 0.25),
            // The neutral point, and the extremes of each range that ARE allowed.
            ("focal:0:0.5", 0.0, 0.5),
            ("focal:5:0.99", 5.0, 0.99),
        ] {
            let mut cfg = MlpConfig::compiled_default();
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_LOSS", spelling)]));
            assert_eq!(cfg.loss, MlpLoss::Focal { gamma, alpha }, "{spelling:?}");
        }
    }

    /// Overriding the loss moves ONLY the loss, and the audit trail quotes the
    /// value the operator actually set — the string a sweep table is labelled by.
    #[test]
    fn loss_override_moves_only_the_loss() {
        let base = MlpConfig::compiled_default();
        let mut cfg = base.clone();
        let applied = cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_LOSS", "focal:2:0.25")]));
        assert_eq!(applied, vec!["TIMSSEEK_MLP_LOSS=\"focal:2:0.25\""]);
        assert_eq!(
            MlpConfig {
                loss: base.loss,
                ..cfg
            },
            base,
            "overriding `loss` moved some other field too"
        );
    }

    /// The loss abort names the variable, the offending value and the grammar,
    /// including which of the two numbers has which range — a `focal:2:1` that
    /// said only "malformed" would leave the operator re-reading the source.
    #[test]
    fn the_loss_abort_names_the_variable_the_value_and_the_format() {
        let e = std::panic::catch_unwind(|| {
            let mut cfg = MlpConfig::compiled_default();
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_LOSS", "focal:2:1")]));
        })
        .expect_err("an out-of-range alpha must abort");
        let msg = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or("<not a string payload>");
        for needle in [
            "TIMSSEEK_MLP_LOSS",
            "focal:2:1",
            "focal:GAMMA:ALPHA",
            "GAMMA >= 0",
            "(0, 1)",
        ] {
            assert!(
                msg.contains(needle),
                "abort message lacks {needle:?}: {msg}"
            );
        }
    }

    /// The tolerant spellings of the activation, which a shell script is likely
    /// to produce. Everything outside this set aborts.
    #[test]
    fn activation_spellings_parse() {
        for spelling in ["square", "SQUARE", " Square "] {
            let mut cfg = MlpConfig::compiled_default();
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_ACTIVATION", spelling)]));
            assert_eq!(cfg.activation, Activation::Square, "{spelling:?}");
        }
        for spelling in ["leaky_relu", "leaky-relu", "LEAKY_RELU"] {
            let mut cfg = MlpConfig::compiled_default();
            cfg.activation = Activation::Square;
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_ACTIVATION", spelling)]));
            assert_eq!(cfg.activation, Activation::LeakyRelu, "{spelling:?}");
        }
    }

    /// The configured activation reaches the BUILT NET, and each variant builds
    /// the layer it names. Asserted against hand-assembled stacks rather than
    /// against each other: the hand-built arm draws from the RNG in the same
    /// order, so agreement is bit-exact and a plumbing mistake (either variant
    /// wired to the other's layer, or the choice dropped on the floor) fails.
    ///
    /// This is also the DETERMINISM check for the new field. Neither activation
    /// draws from the RNG, so the `LeakyRelu` stack is bit-identical to what
    /// `feedforward` produced before the choice existed — which is what the
    /// hand-built `Relu::leaky(_, LEAKY_SLOPE)` arm restates.
    #[test]
    fn the_configured_activation_reaches_the_built_net() {
        let x = filled(4, 3, &mut StdRng::seed_from_u64(9));
        let logits = |m: &mut Mlp| -> Vec<u32> {
            let out = m.forward(&x, false);
            (0..out.rows()).map(|i| out.row(i)[0].to_bits()).collect()
        };

        for (act, hand) in [
            (
                Activation::LeakyRelu,
                Box::new(|d: usize| -> Box<dyn Layer> { Box::new(Relu::leaky(d, LEAKY_SLOPE)) })
                    as Box<dyn Fn(usize) -> Box<dyn Layer>>,
            ),
            (
                Activation::Square,
                Box::new(|d: usize| -> Box<dyn Layer> { Box::new(Square::new(d)) }),
            ),
        ] {
            let mut r1 = StdRng::seed_from_u64(5);
            let mut built = Mlp::feedforward(3, &[8, 4], act, &mut r1);

            let mut r2 = StdRng::seed_from_u64(5);
            let mut expected = Mlp::new(vec![
                Box::new(Linear::new(3, 8, &mut r2)),
                hand(8),
                Box::new(Linear::new(8, 4, &mut r2)),
                hand(4),
                Box::new(Linear::new(4, 1, &mut r2)),
            ]);

            assert_eq!(
                logits(&mut built),
                logits(&mut expected),
                "{act:?} did not build the stack it names"
            );
        }

        // Non-vacuity: the two stacks are not the same net, so the equalities
        // above are not both true of everything.
        let mut r = StdRng::seed_from_u64(5);
        let mut leaky = Mlp::feedforward(3, &[8, 4], Activation::LeakyRelu, &mut r);
        let mut r = StdRng::seed_from_u64(5);
        let mut square = Mlp::feedforward(3, &[8, 4], Activation::Square, &mut r);
        assert_ne!(logits(&mut leaky), logits(&mut square));
    }

    /// The abort names the variable, the offending value AND the expected
    /// format. A panic that says only "invalid value" would leave an operator
    /// guessing which of nine variables their sweep script mangled.
    #[test]
    fn the_abort_names_the_variable_the_value_and_the_format() {
        let e = std::panic::catch_unwind(|| {
            let mut cfg = MlpConfig::compiled_default();
            cfg.apply_overrides(env_of(&[("TIMSSEEK_MLP_LR", "banana")]));
        })
        .expect_err("a malformed value must abort");
        let msg = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or("<not a string payload>");
        for needle in ["TIMSSEEK_MLP_LR", "banana", "finite positive float"] {
            assert!(
                msg.contains(needle),
                "abort message lacks {needle:?}: {msg}"
            );
        }
    }

    /// Under `cfg(test)` — which is how this assertion is compiled —
    /// `MlpConfig::default()` MUST NOT read the environment, so the many
    /// fixtures built as `..MlpConfig::default()` cannot be retuned by a stray
    /// variable in the developer's shell.
    ///
    /// A tripwire, not a proof: it is tautological while the `cfg(test)` arm
    /// stands, and starts failing the moment someone removes that arm while a
    /// `TIMSSEEK_MLP_*` is exported.
    #[test]
    fn default_ignores_the_environment_in_tests() {
        assert_eq!(MlpConfig::default(), MlpConfig::compiled_default());
    }
}
