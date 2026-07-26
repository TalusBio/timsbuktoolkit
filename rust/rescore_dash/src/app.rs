//! Dashboard state: which tab, which feature, which transforms, and the key
//! bindings that move between them. Rendering lives in [`crate::ui`].

use crate::curves;
use crate::stats::{
    self,
    FeatureSummary,
};
use crate::transform::{
    XTransform,
    YTransform,
};
use crate::view::RescoreView;
use ratatui::crossterm::event::{
    KeyCode,
    KeyEvent,
    KeyModifiers,
};
use ratatui::widgets::TableState;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Overview,
    Fdr,
    Calibration,
    Features,
}

impl Tab {
    pub const ALL: [Tab; 4] = [Self::Overview, Self::Fdr, Self::Calibration, Self::Features];

    pub fn title(self) -> &'static str {
        match self {
            Self::Overview => "Overview",
            Self::Fdr => "FDR",
            Self::Calibration => "Calibration",
            Self::Features => "Features",
        }
    }

    fn shift(self, delta: isize) -> Self {
        let i = Self::ALL.iter().position(|t| *t == self).unwrap_or(0) as isize;
        let n = Self::ALL.len() as isize;
        Self::ALL[((i + delta).rem_euclid(n)) as usize]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortKey {
    Name,
    TargetMean,
    DecoyMean,
    Auc,
    CohensD,
    NanFrac,
    Gain,
}

impl SortKey {
    const CYCLE: [SortKey; 7] = [
        Self::Name,
        Self::TargetMean,
        Self::DecoyMean,
        Self::Auc,
        Self::CohensD,
        Self::NanFrac,
        Self::Gain,
    ];

    pub fn next(self) -> Self {
        let i = Self::CYCLE.iter().position(|k| *k == self).unwrap_or(0);
        Self::CYCLE[(i + 1) % Self::CYCLE.len()]
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Name => "name",
            Self::TargetMean => "target mean",
            Self::DecoyMean => "decoy mean",
            Self::Auc => "AUC",
            Self::CohensD => "|d|",
            Self::NanFrac => "NaN%",
            Self::Gain => "gain",
        }
    }
}

pub enum Flow {
    Continue,
    Quit,
}

pub struct App<'a> {
    view: &'a RescoreView<'a>,
    tab: Tab,
    x: XTransform,
    y: YTransform,
    clip: bool,
    sort: SortKey,
    /// Direction, whose "forward" reading differs by key: for `SortKey::Name`
    /// `true` is ascending alphabetical (the natural reading order for text),
    /// while for every numeric key `true` is largest-first (the natural
    /// reading order for "what stands out"). `false` reverses either one.
    sort_desc: bool,
    filter: String,
    filter_editing: bool,
    /// Feature indices passing the filter, in sort order.
    visible: Vec<usize>,
    /// Position within `visible`, not a feature index.
    cursor: usize,
    cache: HashMap<usize, FeatureSummary>,
    /// Targets-passing-vs-q-value curve, cached: it depends only on
    /// `view.qvalue`/`view.is_target`, neither of which ever changes, so it is
    /// worth computing once rather than re-sorting on every keystroke.
    qvalue_curve: Option<Vec<(f64, f64)>>,
    /// Decoy/target PP curve, cached for the same reason.
    pp_curve: Option<Vec<(f64, f64)>>,
    /// Exact rank-based AUC of the discriminant score, cached: it depends only
    /// on `view.score`/`view.is_target` and sorts its input, so it is worth
    /// computing once rather than on every Overview frame.
    auc: Option<f64>,
    /// Scroll/selection state for the features table, kept on `App` rather
    /// than rebuilt per frame: `TableState::offset` is what lets a `Table`
    /// scroll, and a fresh `TableState::default()` every render recomputes it
    /// from `0`, pinning the viewport instead of tracking the cursor.
    table_state: TableState,
}

impl<'a> App<'a> {
    pub fn new(view: &'a RescoreView<'a>) -> Self {
        let mut app = Self {
            view,
            tab: Tab::Overview,
            x: XTransform::Linear,
            y: YTransform::Density,
            clip: true,
            sort: SortKey::Name,
            sort_desc: true,
            filter: String::new(),
            filter_editing: false,
            visible: Vec::new(),
            cursor: 0,
            cache: HashMap::new(),
            qvalue_curve: None,
            pp_curve: None,
            auc: None,
            table_state: TableState::default(),
        };
        app.refresh_visible();
        app
    }

    pub fn view(&self) -> &'a RescoreView<'a> {
        self.view
    }

    pub fn tab(&self) -> Tab {
        self.tab
    }

    pub fn x(&self) -> XTransform {
        self.x
    }

    pub fn y(&self) -> YTransform {
        self.y
    }

    pub fn clip(&self) -> bool {
        self.clip
    }

    pub fn sort_key(&self) -> SortKey {
        self.sort
    }

    pub fn filter(&self) -> &str {
        &self.filter
    }

    pub fn filter_editing(&self) -> bool {
        self.filter_editing
    }

    pub fn visible(&self) -> &[usize] {
        &self.visible
    }

    /// Cursor position within [`Self::visible`], for the table widget.
    pub fn cursor(&self) -> usize {
        self.cursor
    }

    /// Feature index under the cursor, or `None` when the filter matches
    /// nothing.
    pub fn selected_feature(&self) -> Option<usize> {
        self.visible.get(self.cursor).copied()
    }

    /// Summary for feature `j`, computed on first access and cached. Filling
    /// all ~131 up front would cost a full matrix walk before the first frame.
    ///
    /// `feature_column` is `Clone` (it is a cheap strided iterator over the
    /// borrowed matrix), so `summarize` can walk it directly rather than
    /// paying for an intermediate `Vec` per column.
    pub fn summary(&mut self, j: usize) -> FeatureSummary {
        if let Some(s) = self.cache.get(&j) {
            return *s;
        }
        let s = stats::summarize(self.view.feature_column(j), &self.view.is_target);
        self.cache.insert(j, s);
        s
    }

    /// Targets-passing-vs-q-value curve, computed on first access and
    /// cached, since `curves::qvalue_curve` sorts its input and neither
    /// `view.qvalue` nor `view.is_target` ever changes. Fixed at 100 points,
    /// matching the FDR panel's original per-frame call.
    pub fn qvalue_curve(&mut self) -> Vec<(f64, f64)> {
        self.qvalue_curve
            .get_or_insert_with(|| {
                curves::qvalue_curve(&self.view.qvalue, &self.view.is_target, 100)
            })
            .clone()
    }

    /// Decoy/target PP curve, computed on first access and cached. Fixed at
    /// 200 points, matching the calibration panel's original per-frame call.
    pub fn pp_curve(&mut self) -> Vec<(f64, f64)> {
        self.pp_curve
            .get_or_insert_with(|| curves::pp_curve(&self.view.score, &self.view.is_target, 200))
            .clone()
    }

    /// Exact rank-based AUC of the discriminant score, computed on first
    /// access and cached: `stats::auc_exact` sorts its input, so recomputing
    /// it on every Overview frame is a full sort of every row per keystroke.
    pub fn auc(&mut self) -> f64 {
        *self.auc.get_or_insert_with(|| {
            stats::auc_exact(
                self.view.score.iter().map(|&s| s as f64),
                &self.view.is_target,
            )
        })
    }

    /// The features table's `TableState`, synced to the current cursor before
    /// each render. Returning the same stored `TableState` (rather than a
    /// fresh `TableState::default()`) is what lets ratatui's `Table` keep its
    /// scroll offset between frames.
    pub fn table_state(&mut self) -> &mut TableState {
        let selected = if self.visible.is_empty() {
            None
        } else {
            Some(self.cursor)
        };
        self.table_state.select(selected);
        &mut self.table_state
    }

    pub fn handle_key(&mut self, key: KeyEvent) -> Flow {
        // Raw mode disables the terminal's own SIGINT handling, so Ctrl-C
        // arrives as a plain key event rather than a signal. Handled before
        // anything else — including the filter-editing dispatch below — so
        // it always quits rather than being typed into the filter box or
        // toggling clip (`KeyCode::Char('c')` with no modifiers).
        if is_ctrl_c(key) {
            return Flow::Quit;
        }
        if self.filter_editing {
            self.handle_filter_key(key);
            return Flow::Continue;
        }
        // CONTROL/ALT combinations are reserved (Ctrl-C above, and headroom
        // for future bindings); only a bare or Shifted character triggers a
        // command, so e.g. Ctrl-X does not cycle the x-transform.
        let plain = !key
            .modifiers
            .intersects(KeyModifiers::CONTROL | KeyModifiers::ALT);
        match key.code {
            KeyCode::Char('q') | KeyCode::Esc if plain => return Flow::Quit,
            KeyCode::Char('l') | KeyCode::Right | KeyCode::Tab if plain => {
                self.tab = self.tab.shift(1)
            }
            KeyCode::Char('h') | KeyCode::Left | KeyCode::BackTab if plain => {
                self.tab = self.tab.shift(-1)
            }
            KeyCode::Char('j') | KeyCode::Down if plain => {
                self.cursor = (self.cursor + 1).min(self.visible.len().saturating_sub(1));
            }
            KeyCode::Char('k') | KeyCode::Up if plain => {
                self.cursor = self.cursor.saturating_sub(1)
            }
            KeyCode::Char('x') if plain => self.x = self.x.next(),
            KeyCode::Char('X') if plain => self.x = self.x.prev(),
            KeyCode::Char('y') if plain => self.y = self.y.next(),
            KeyCode::Char('Y') if plain => self.y = self.y.prev(),
            KeyCode::Char('c') if plain => self.clip = !self.clip,
            KeyCode::Char('s') if plain => {
                self.sort = self.sort.next();
                self.refresh_visible();
            }
            KeyCode::Char('S') if plain => {
                self.sort_desc = !self.sort_desc;
                self.refresh_visible();
            }
            KeyCode::Char('/') if plain => {
                self.filter_editing = true;
                self.filter.clear();
                // Without this, the help line shows an empty filter while the
                // table still shows whatever the previous filter left
                // visible, until the user types the first character.
                self.refresh_visible();
            }
            _ => {}
        }
        Flow::Continue
    }

    fn handle_filter_key(&mut self, key: KeyEvent) {
        let plain = !key
            .modifiers
            .intersects(KeyModifiers::CONTROL | KeyModifiers::ALT);
        match key.code {
            KeyCode::Esc => {
                self.filter.clear();
                self.filter_editing = false;
                self.refresh_visible();
            }
            KeyCode::Enter => {
                self.filter_editing = false;
                self.refresh_visible();
            }
            KeyCode::Backspace => {
                self.filter.pop();
            }
            KeyCode::Char(c) if plain => self.filter.push(c),
            _ => {}
        }
    }

    /// Rebuild the visible feature list from the filter and sort key, keeping
    /// the cursor on the same feature when it survives.
    fn refresh_visible(&mut self) {
        let previous = self.selected_feature();
        let needle = self.filter.to_lowercase();
        let mut idx: Vec<usize> = (0..self.view.n_features())
            .filter(|&j| {
                needle.is_empty() || self.view.feature_names[j].to_lowercase().contains(&needle)
            })
            .collect();

        if self.sort == SortKey::Name {
            idx.sort_by(|&a, &b| self.view.feature_names[a].cmp(&self.view.feature_names[b]));
            if !self.sort_desc {
                idx.reverse();
            }
        } else {
            // Numeric sorts need the summaries; this is the one place the cache
            // is filled eagerly, and only for the filtered subset. Built as a
            // plain loop (rather than `.map` inside the closure) because
            // `self.summary(j)` takes `&mut self` and must not overlap with any
            // other borrow of `self` in the same expression.
            let mut keyed: Vec<(usize, f64)> = Vec::with_capacity(idx.len());
            for &j in &idx {
                let s = self.summary(j);
                let k = match self.sort {
                    SortKey::TargetMean => s.target_mean,
                    SortKey::DecoyMean => s.decoy_mean,
                    SortKey::Auc => s.auc,
                    SortKey::CohensD => s.cohens_d,
                    SortKey::NanFrac => s.nan_frac,
                    SortKey::Gain => self.view.mean_gain(&self.view.feature_names[j]) as f64,
                    SortKey::Name => 0.0,
                };
                keyed.push((j, k));
            }
            // `sort_desc` picks the finite-value order; NaN sorts last either
            // way. That rule has to live in the comparator itself, not in a
            // value substitution followed by a blanket `idx.reverse()` below
            // — a later reverse of the whole vector would undo a NaN-last
            // placement exactly when `sort_desc` is false.
            keyed.sort_by(|a, b| cmp_nan_last(a.1, b.1, self.sort_desc));
            idx = keyed.into_iter().map(|(j, _)| j).collect();
        }

        self.visible = idx;
        self.cursor = previous
            .and_then(|j| self.visible.iter().position(|&v| v == j))
            .unwrap_or(0)
            .min(self.visible.len().saturating_sub(1));
    }
}

/// Open the dashboard and block until the user quits.
///
/// Skips silently (with a warning) when stdout is not a terminal, so a piped
/// or containerized run is unaffected. Terminal setup failures are warn-only
/// too: if `try_init` fails (e.g. no usable controlling terminal even though
/// stdout is a tty — `setsid`, some container/CI pty setups), the terminal is
/// restored and `run` returns `Ok(())` rather than propagating the error or
/// panicking (unlike `ratatui::init`, which is `try_init().expect(...)` and
/// would unwind through the caller, past a warn-only `if let Err`, after
/// Phase 6 has already written output).
///
/// A panic inside the event loop is caught by [`catch_panics`] and turned
/// into a warning too — but see that function's doc comment: this guard is
/// inert in a release build, so a release-build panic in the loop still
/// aborts the process. What `run` actually guarantees unconditionally is that
/// terminal setup never panics or fails the caller, and that the terminal is
/// restored on every path that returns (a panic that reaches `catch_unwind`
/// at all triggers `try_init`'s panic hook first, which itself restores the
/// terminal, so the explicit `ratatui::restore()` below is belt-and-braces —
/// harmless if the hook already ran, load-bearing if it didn't).
pub fn run(view: &RescoreView<'_>) -> std::io::Result<()> {
    use std::io::IsTerminal;

    if !std::io::stdout().is_terminal() {
        tracing::warn!("rescore dashboard requested but stdout is not a terminal; skipping");
        return Ok(());
    }
    if let Err(e) = view.validate() {
        tracing::warn!("rescore dashboard input rejected: {e}");
        return Ok(());
    }

    let mut terminal = match ratatui::try_init() {
        Ok(terminal) => terminal,
        Err(e) => {
            // `try_init` can fail after partially succeeding (e.g. raw mode
            // enabled, then `EnterAlternateScreen` fails), so restore
            // unconditionally rather than assuming there is nothing to undo.
            // `restore()` swallows its own errors (it can only print them),
            // so this can't itself introduce a new failure path.
            tracing::warn!("rescore dashboard failed to initialize the terminal: {e}");
            ratatui::restore();
            return Ok(());
        }
    };
    let result = catch_panics(|| event_loop(&mut terminal, view));
    ratatui::restore();
    result
}

/// Runs `f`, converting a panic into a `tracing::warn!` + `Ok(())` rather than
/// letting it unwind — in builds where panics unwind at all. This workspace's
/// top-level `Cargo.toml` sets `[profile.release] panic = "abort"`, so in the
/// shipped release binary a panic inside `f` aborts the process immediately;
/// `catch_unwind` never gets a chance to run, and this function's recovery
/// does not apply there. It is real in dev/test builds (`panic = "unwind"`,
/// the default profile.dev), which is what this module's tests exercise.
/// Kept anyway because unwinding is strictly better than nothing when it is
/// available, and because a future decision to drop `panic = "abort"` should
/// not require re-adding this. Factored out of [`run`] so the recovery
/// behavior is testable without a real terminal: [`event_loop`] itself needs
/// one, but this wrapper does not.
fn catch_panics(f: impl FnOnce() -> std::io::Result<()>) -> std::io::Result<()> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(result) => result,
        Err(_) => {
            tracing::warn!(
                "rescore dashboard panicked; recovering so the search run is not aborted"
            );
            Ok(())
        }
    }
}

/// Whether `key` is Ctrl-C: raw mode disables the terminal's own `ISIG`
/// handling, so Ctrl-C arrives as this key event rather than delivering
/// `SIGINT`.
fn is_ctrl_c(key: KeyEvent) -> bool {
    key.code == KeyCode::Char('c') && key.modifiers.contains(KeyModifiers::CONTROL)
}

fn event_loop<B: ratatui::backend::Backend>(
    terminal: &mut ratatui::Terminal<B>,
    view: &RescoreView<'_>,
) -> std::io::Result<()> {
    use ratatui::crossterm::event::{
        self,
        Event,
        KeyEventKind,
    };

    let mut app = App::new(view);
    loop {
        terminal.draw(|f| crate::ui::draw(f, &mut app))?;
        // Only key *presses*: on Windows crossterm also emits releases, which
        // would double every keystroke.
        if let Event::Key(key) = event::read()?
            && key.kind == KeyEventKind::Press
            && matches!(app.handle_key(key), Flow::Quit)
        {
            return Ok(());
        }
    }
}

/// Orders two sort-key values with `NaN` always last, independent of `desc`.
/// `desc` breaks ties between two finite values only: `true` is largest
/// first, `false` is smallest first.
fn cmp_nan_last(a: f64, b: f64, desc: bool) -> std::cmp::Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) if desc => b.partial_cmp(&a).unwrap_or(std::cmp::Ordering::Equal),
        (false, false) => a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };
    use std::sync::Arc;

    fn key(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE)
    }

    fn code(c: KeyCode) -> KeyEvent {
        KeyEvent::new(c, KeyModifiers::NONE)
    }

    fn ctrl(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL)
    }

    static MATRIX: [f64; 6] = [1.0, 10.0, 2.0, 20.0, 3.0, 30.0];

    fn view() -> RescoreView<'static> {
        RescoreView {
            feature_names: vec![Arc::from("alpha_score"), Arc::from("beta_count")],
            features: &MATRIX,
            is_target: vec![true, false, true],
            score: vec![3.0, 1.0, 2.0],
            qvalue: vec![0.001, 0.9, 0.02],
            importance: Vec::new(),
        }
    }

    #[test]
    fn h_and_l_move_between_tabs() {
        let v = view();
        let mut app = App::new(&v);
        assert_eq!(app.tab(), Tab::Overview);
        app.handle_key(key('l'));
        assert_eq!(app.tab(), Tab::Fdr);
        app.handle_key(key('h'));
        assert_eq!(app.tab(), Tab::Overview);
        app.handle_key(key('h'));
        assert_eq!(app.tab(), Tab::Features, "wraps backwards");
    }

    #[test]
    fn j_and_k_move_the_feature_selection_and_clamp() {
        let v = view();
        let mut app = App::new(&v);
        assert_eq!(app.selected_feature(), Some(0));
        app.handle_key(key('j'));
        assert_eq!(app.selected_feature(), Some(1));
        app.handle_key(key('j'));
        assert_eq!(app.selected_feature(), Some(1), "clamps at the end");
        app.handle_key(key('k'));
        app.handle_key(key('k'));
        assert_eq!(app.selected_feature(), Some(0), "clamps at the start");
    }

    #[test]
    fn x_and_y_cycle_transforms_in_both_directions() {
        let v = view();
        let mut app = App::new(&v);
        let x0 = app.x();
        app.handle_key(key('x'));
        assert_ne!(app.x(), x0);
        app.handle_key(KeyEvent::new(KeyCode::Char('X'), KeyModifiers::SHIFT));
        assert_eq!(app.x(), x0);

        let y0 = app.y();
        app.handle_key(key('y'));
        assert_ne!(app.y(), y0);
        app.handle_key(KeyEvent::new(KeyCode::Char('Y'), KeyModifiers::SHIFT));
        assert_eq!(app.y(), y0);
    }

    #[test]
    fn c_toggles_clipping_and_q_quits() {
        let v = view();
        let mut app = App::new(&v);
        assert!(app.clip(), "percentile clip is the default");
        app.handle_key(key('c'));
        assert!(!app.clip());
        assert!(matches!(app.handle_key(key('q')), Flow::Quit));
    }

    #[test]
    fn ctrl_c_quits_from_normal_mode() {
        let v = view();
        let mut app = App::new(&v);
        assert!(
            matches!(app.handle_key(ctrl('c')), Flow::Quit),
            "Ctrl-C must quit rather than toggling clip"
        );
    }

    #[test]
    fn ctrl_c_quits_while_filter_editing() {
        let v = view();
        let mut app = App::new(&v);
        app.handle_key(key('/'));
        assert!(app.filter_editing());
        assert!(
            matches!(app.handle_key(ctrl('c')), Flow::Quit),
            "Ctrl-C must quit even mid-filter, not type a literal 'c'"
        );
    }

    #[test]
    fn ctrl_x_does_not_cycle_the_x_transform() {
        let v = view();
        let mut app = App::new(&v);
        let x0 = app.x();
        assert!(matches!(app.handle_key(ctrl('x')), Flow::Continue));
        assert_eq!(app.x(), x0, "Ctrl-X is reserved, not a transform cycle");
    }

    #[test]
    fn slash_filters_feature_names_and_esc_clears() {
        let v = view();
        let mut app = App::new(&v);
        assert_eq!(app.visible().len(), 2);

        app.handle_key(key('/'));
        assert!(app.filter_editing());
        for c in "beta".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert!(!app.filter_editing());
        assert_eq!(app.visible(), &[1], "only beta_count matches");
        assert_eq!(
            app.selected_feature(),
            Some(1),
            "selection follows the filter"
        );

        app.handle_key(key('/'));
        app.handle_key(code(KeyCode::Esc));
        assert_eq!(app.visible().len(), 2, "esc clears the filter");
    }

    /// Reopening the filter box clears `self.filter` immediately; `visible`
    /// must refresh in lockstep so the table does not keep showing the
    /// previous (now-stale) filtered set while the help line already shows an
    /// empty query.
    #[test]
    fn reopening_the_filter_refreshes_visible_immediately() {
        let v = view();
        let mut app = App::new(&v);
        app.handle_key(key('/'));
        for c in "beta".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert_eq!(app.visible(), &[1], "filtered down to beta_count");

        app.handle_key(key('/'));
        assert_eq!(app.filter(), "", "filter text cleared on reopen");
        assert_eq!(
            app.visible().len(),
            2,
            "visible must refresh to match the cleared filter, not stay stale"
        );
    }

    #[test]
    fn typing_while_filtering_does_not_trigger_commands() {
        let v = view();
        let mut app = App::new(&v);
        let tab0 = app.tab();
        app.handle_key(key('/'));
        app.handle_key(key('l'));
        app.handle_key(key('q'));
        assert_eq!(app.tab(), tab0, "keys are literal text while filtering");
        assert_eq!(app.filter(), "lq");
    }

    #[test]
    fn a_filter_matching_nothing_leaves_no_selection() {
        let v = view();
        let mut app = App::new(&v);
        app.handle_key(key('/'));
        for c in "zzz".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert!(app.visible().is_empty());
        assert_eq!(app.selected_feature(), None);
        // Must not panic.
        app.handle_key(key('j'));
    }

    #[test]
    fn s_cycles_the_sort_key_and_reorders() {
        let v = view();
        let mut app = App::new(&v);
        let first = app.visible().to_vec();
        // Sort by target mean descending: beta_count (mean 20) before alpha_score (mean 2).
        while app.sort_key() != SortKey::TargetMean {
            app.handle_key(key('s'));
        }
        assert_ne!(app.visible(), first.as_slice());
        assert_eq!(app.visible()[0], 1);
    }

    #[test]
    fn summary_is_cached_and_correct() {
        let v = view();
        let mut app = App::new(&v);
        let s = app.summary(0);
        assert_eq!(s.target_mean, 2.0, "rows 0 and 2 are targets: (1+3)/2");
        assert_eq!(s.decoy_mean, 2.0);
        assert_eq!(app.summary(0), s, "second call returns the cached value");
    }

    /// The FDR and PP curves sort their inputs, so they must be computed once
    /// per session rather than once per frame. `qvalue_curve`/`pp_curve` don't
    /// depend on any transform, clip, filter, or selection state, so a second
    /// call must both return the same data and leave the cache field filled
    /// rather than recomputing.
    #[test]
    fn qvalue_curve_is_cached() {
        let v = view();
        let mut app = App::new(&v);
        assert!(
            app.qvalue_curve.is_none(),
            "not computed until first access"
        );

        let first = app.qvalue_curve();
        assert!(app.qvalue_curve.is_some(), "cached after the first access");

        let second = app.qvalue_curve();
        assert_eq!(first, second, "cached data is returned, not recomputed");
    }

    #[test]
    fn pp_curve_is_cached() {
        let v = view();
        let mut app = App::new(&v);
        assert!(app.pp_curve.is_none(), "not computed until first access");

        let first = app.pp_curve();
        assert!(app.pp_curve.is_some(), "cached after the first access");

        let second = app.pp_curve();
        assert_eq!(first, second, "cached data is returned, not recomputed");
    }

    /// `auc_exact` sorts its input, so it must be computed once per session
    /// rather than once per Overview frame, exactly like the FDR/PP curves
    /// above.
    #[test]
    fn auc_is_cached() {
        let v = view();
        let mut app = App::new(&v);
        assert!(app.auc.is_none(), "not computed until first access");

        let first = app.auc();
        assert!(app.auc.is_some(), "cached after the first access");

        let second = app.auc();
        assert_eq!(first, second, "cached data is returned, not recomputed");
    }

    /// A fresh `TableState::default()` every frame resets `offset` to `0`, so
    /// the viewport can never scroll. `table_state()` must hand out the same
    /// stored state across calls, with `select()` following the cursor.
    #[test]
    fn table_state_is_shared_across_accesses_and_tracks_the_cursor() {
        let v = view();
        let mut app = App::new(&v);
        assert_eq!(app.table_state().selected(), Some(0));

        app.handle_key(key('j'));
        assert_eq!(
            app.table_state().selected(),
            Some(1),
            "table_state's selection follows the cursor"
        );
    }

    /// `nan_feat` is NaN in every row, so `stats::summarize` gives it a NaN
    /// AUC (no finite value contributes to either class): a real, naturally
    /// occurring case, not a synthetic key. It must sort last regardless of
    /// direction.
    #[test]
    fn nan_summary_keys_sort_last_in_both_directions() {
        static MATRIX: [f64; 12] = [
            1.0,
            10.0,
            f64::NAN,
            2.0,
            20.0,
            f64::NAN,
            3.0,
            30.0,
            f64::NAN,
            4.0,
            40.0,
            f64::NAN,
        ];
        let v = RescoreView {
            feature_names: vec![Arc::from("low"), Arc::from("high"), Arc::from("nan_feat")],
            features: &MATRIX,
            is_target: vec![true, false, true, false],
            score: vec![0.0; 4],
            qvalue: vec![0.0; 4],
            importance: Vec::new(),
        };
        let mut app = App::new(&v);
        while app.sort_key() != SortKey::Auc {
            app.handle_key(key('s'));
        }

        assert_eq!(
            app.visible().last().copied(),
            Some(2),
            "descending: NaN AUC still sorts last"
        );

        app.handle_key(key('S'));
        assert_eq!(
            app.visible().last().copied(),
            Some(2),
            "ascending: NaN AUC still sorts last"
        );
    }

    /// The test harness's stdout is not a terminal (captured, or at least not
    /// a real tty in CI), so `run` must take the early-return path rather
    /// than touching a real terminal, and it must return `Ok`, not panic.
    #[test]
    fn run_returns_ok_and_does_not_panic_when_stdout_is_not_a_terminal() {
        let v = view();
        assert!(super::run(&v).is_ok());
    }

    /// `catch_panics` is the piece of `run` that keeps a panic inside the
    /// event loop from aborting the whole search run — in builds where
    /// panics unwind, which is what this test runs under (`panic = "unwind"`,
    /// the default dev/test profile). `[profile.release]` sets `panic =
    /// "abort"` workspace-wide, where this guard cannot run at all; see
    /// `catch_panics`'s doc comment. Exercised directly here since, unlike
    /// `event_loop`, it needs no real terminal.
    #[test]
    fn catch_panics_converts_a_panic_into_ok() {
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {})); // silence the default panic printout
        let result = catch_panics(|| panic!("simulated event-loop panic"));
        std::panic::set_hook(prev_hook);
        assert!(result.is_ok(), "a caught panic must not fail the run");
    }

    #[test]
    fn catch_panics_passes_through_a_successful_result() {
        let result = catch_panics(|| Ok(()));
        assert!(result.is_ok());
    }

    #[test]
    fn catch_panics_passes_through_an_io_error() {
        let result = catch_panics(|| Err(std::io::Error::other("simulated io error")));
        assert!(result.is_err(), "a real Err must not be swallowed");
    }
}
