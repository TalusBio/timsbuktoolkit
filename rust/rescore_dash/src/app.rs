//! Dashboard state: which tab, which feature, which transforms, and the key
//! bindings that move between them.
//!
//! [`App`] is a [`Dashboard`] — everything materialized, see
//! [`crate::precompute`] — plus the handful of fields a keystroke can change.
//! Nothing here computes over the data; a key press moves an index.
//! Rendering lives in [`crate::ui`].

use crate::cycle;
use crate::precompute::{
    DEFAULT_Q_ZOOM,
    Dashboard,
    FeatureColumn,
};
use crate::transform::{
    Axis,
    YTransform,
};
use ratatui::crossterm::event::{
    KeyCode,
    KeyEvent,
    KeyModifiers,
};
use ratatui::widgets::TableState;

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
        cycle::step(&Self::ALL, self, delta)
    }
}

pub enum Flow {
    Continue,
    Quit,
}

/// Plot points, reused across frames.
///
/// `ratatui::Dataset::data` borrows a `&[(f64, f64)]`, so the points have to
/// outlive the widget; they live here rather than in a vector collected inside
/// the draw call.
#[derive(Default)]
pub(crate) struct Scratch {
    pub(crate) target: Vec<(f64, f64)>,
    pub(crate) decoy: Vec<(f64, f64)>,
}

pub struct App {
    pub(crate) dash: Dashboard,
    pub(crate) scratch: Scratch,
    tab: Tab,
    x: Axis,
    y: YTransform,
    clip: bool,
    /// Index into the dashboard's FDR zoom levels, not a q-value. Kept as an
    /// index so `z` is a bounded cycle rather than free-form zooming into a
    /// range no curve was gridded for.
    q_zoom: usize,
    sort: FeatureColumn,
    /// Direction, whose "forward" reading differs by key: for `FeatureColumn::Name`
    /// `true` is ascending alphabetical (the natural reading order for text),
    /// while for every numeric key `true` is largest-first (the natural
    /// reading order for "what stands out"). `false` reverses either one.
    sort_desc: bool,
    filter: String,
    filter_editing: bool,
    /// Feature indices passing the filter, in sort order.
    pub(crate) visible: Vec<usize>,
    /// Position within `visible`, not a feature index.
    pub(crate) cursor: usize,
    /// Scroll/selection state for the features table, kept on `App` rather
    /// than rebuilt per frame: `TableState::offset` is what lets a `Table`
    /// scroll, and a fresh `TableState::default()` every render recomputes it
    /// from `0`, pinning the viewport instead of tracking the cursor.
    ///
    /// Reachable directly from [`crate::ui`] because the table widget borrows
    /// `dash` while the state needs `&mut`; a whole-`App` accessor would make
    /// those two borrows overlap.
    pub(crate) table_state: TableState,
}

impl App {
    pub fn new(dash: Dashboard) -> Self {
        let mut app = Self {
            dash,
            scratch: Scratch::default(),
            tab: Tab::Overview,
            x: Axis::ALL[0],
            y: YTransform::Density,
            clip: true,
            q_zoom: DEFAULT_Q_ZOOM,
            sort: FeatureColumn::Name,
            sort_desc: true,
            filter: String::new(),
            filter_editing: false,
            visible: Vec::new(),
            cursor: 0,
            table_state: TableState::default(),
        };
        app.refresh_visible();
        app
    }

    pub fn dashboard(&self) -> &Dashboard {
        &self.dash
    }

    pub fn tab(&self) -> Tab {
        self.tab
    }

    pub fn x(&self) -> Axis {
        self.x
    }

    pub fn y(&self) -> YTransform {
        self.y
    }

    pub fn clip(&self) -> bool {
        self.clip
    }

    /// Which of the dashboard's FDR zoom levels the curve panel draws.
    pub fn q_zoom(&self) -> usize {
        self.q_zoom
    }

    pub fn sort_key(&self) -> FeatureColumn {
        self.sort
    }

    pub fn filter(&self) -> &str {
        &self.filter
    }

    pub fn filter_editing(&self) -> bool {
        self.filter_editing
    }

    /// Feature index under the cursor, or `None` when the filter matches
    /// nothing.
    pub fn selected_feature(&self) -> Option<usize> {
        self.visible.get(self.cursor).copied()
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
            KeyCode::Char('z') if plain => self.q_zoom = (self.q_zoom + 1) % self.dash.n_q_zooms(),
            KeyCode::Char('Z') if plain => {
                let n = self.dash.n_q_zooms();
                self.q_zoom = (self.q_zoom + n - 1) % n;
            }
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

    /// The value feature `j` sorts on under the current key. Every one of these
    /// is a precomputed array read.
    fn sort_value(&self, j: usize) -> f64 {
        // `None` is the name column, which the caller sorts as text.
        self.dash.feature_value(j, self.sort).unwrap_or(0.0)
    }

    /// Rebuild the visible feature list from the filter and sort key, keeping
    /// the cursor on the same feature when it survives.
    fn refresh_visible(&mut self) {
        let previous = self.selected_feature();
        let needle = self.filter.to_lowercase();
        let mut idx: Vec<usize> = (0..self.dash.n_features())
            .filter(|&j| {
                needle.is_empty() || self.dash.feature_names[j].to_lowercase().contains(&needle)
            })
            .collect();

        if self.sort == FeatureColumn::Name {
            idx.sort_by(|&a, &b| self.dash.feature_names[a].cmp(&self.dash.feature_names[b]));
            if !self.sort_desc {
                idx.reverse();
            }
        } else {
            // `sort_desc` picks the finite-value order; NaN sorts last either
            // way. That rule has to live in the comparator itself, not in a
            // value substitution followed by a blanket `idx.reverse()` — a
            // later reverse of the whole vector would undo a NaN-last
            // placement exactly when `sort_desc` is false.
            let desc = self.sort_desc;
            idx.sort_by(|&a, &b| cmp_nan_last(self.sort_value(a), self.sort_value(b), desc));
        }

        self.visible = idx;
        self.cursor = previous
            .and_then(|j| self.visible.iter().position(|&v| v == j))
            .unwrap_or(0)
            .min(self.visible.len().saturating_sub(1));
    }
}

/// Whether a dashboard could be shown at all.
///
/// Worth asking *before* [`Dashboard::build`], not after: building is the
/// expensive step, and there is no point paying for it to discover that stdout
/// is a pipe.
pub fn available() -> bool {
    use std::io::IsTerminal;
    std::io::stdout().is_terminal()
}

/// Open the dashboard and block until the user quits.
///
/// Blocks until the user quits.
///
/// Never fails the caller and never panics on setup: a non-terminal stdout or
/// a failed `try_init` warns, restores the terminal and returns `Ok`. This is
/// deliberately unlike `ratatui::init`, which is `try_init().expect(..)` and
/// would unwind past the caller's warn-only `if let Err` — after the results
/// have already been written to disk.
///
/// The terminal is restored on every returning path. [`catch_panics`] covers
/// the event loop, but only in a debug build; see its doc.
pub fn run(dash: Dashboard) -> std::io::Result<()> {
    if !available() {
        tracing::warn!("rescore dashboard requested but stdout is not a terminal; skipping");
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
    let mut app = App::new(dash);
    let result = catch_panics(|| event_loop(&mut terminal, &mut app));
    ratatui::restore();
    result
}

/// Runs `f`, converting a panic into a `tracing::warn!` + `Ok(())`.
///
/// Only where panics unwind. The workspace sets `[profile.release] panic =
/// "abort"`, so in the shipped binary a panic inside `f` aborts before
/// `catch_unwind` sees it; this recovery is real in dev/test builds only.
/// Separate from [`run`] because it needs no terminal, so it is testable.
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
    app: &mut App,
) -> std::io::Result<()> {
    use ratatui::crossterm::event::{
        self,
        Event,
        KeyEventKind,
    };

    loop {
        terminal.draw(|f| crate::ui::draw(f, app))?;
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
pub(crate) mod tests {
    use super::*;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };

    fn key(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE)
    }

    fn code(c: KeyCode) -> KeyEvent {
        KeyEvent::new(c, KeyModifiers::NONE)
    }

    fn ctrl(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL)
    }

    /// Build a dashboard from column-major data, so a test reads as a list of
    /// columns rather than an interleaved matrix literal. Shared with the `ui`
    /// tests.
    pub(crate) fn dashboard(names: &[&str], columns: &[&[f64]]) -> Dashboard {
        crate::precompute::tests::fixture_from_columns(names, columns).build()
    }

    fn app() -> App {
        App::new(dashboard(
            &["alpha_score", "beta_count"],
            &[&[1.0, 2.0, 3.0, 4.0], &[10.0, 20.0, 30.0, 40.0]],
        ))
    }

    #[test]
    fn h_and_l_move_between_tabs() {
        let mut app = app();
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
        let mut app = app();
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
        let mut app = app();
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
        let mut app = app();
        assert!(app.clip(), "percentile clip is the default");
        app.handle_key(key('c'));
        assert!(!app.clip());
        assert!(matches!(app.handle_key(key('q')), Flow::Quit));
    }

    #[test]
    fn z_cycles_the_fdr_zoom_in_both_directions() {
        let mut app = app();
        let n = app.dashboard().n_q_zooms();
        assert!(n > 1, "a single zoom level makes the key pointless");
        assert_ne!(
            app.dashboard().q_curve(app.q_zoom()).zoom,
            1.0,
            "the default view must be zoomed in; the full (0, 1] curve is the \
             one that shows the least"
        );

        let start = app.q_zoom();
        let mut seen = vec![app.dashboard().q_curve(start).zoom];
        for _ in 1..n {
            app.handle_key(key('z'));
            seen.push(app.dashboard().q_curve(app.q_zoom()).zoom);
        }
        app.handle_key(key('z'));
        assert_eq!(app.q_zoom(), start, "z must wrap, not run off the end");

        seen.sort_by(f64::total_cmp);
        seen.dedup();
        assert_eq!(seen.len(), n, "every zoom level must be reachable by z");

        app.handle_key(key('Z'));
        assert_eq!(
            app.q_zoom(),
            (start + n - 1) % n,
            "Z must step the other way"
        );
    }

    #[test]
    fn ctrl_c_quits_from_normal_mode() {
        let mut app = app();
        assert!(
            matches!(app.handle_key(ctrl('c')), Flow::Quit),
            "Ctrl-C must quit rather than toggling clip"
        );
    }

    #[test]
    fn ctrl_c_quits_while_filter_editing() {
        let mut app = app();
        app.handle_key(key('/'));
        assert!(app.filter_editing());
        assert!(
            matches!(app.handle_key(ctrl('c')), Flow::Quit),
            "Ctrl-C must quit even mid-filter, not type a literal 'c'"
        );
    }

    #[test]
    fn ctrl_x_does_not_cycle_the_x_transform() {
        let mut app = app();
        let x0 = app.x();
        assert!(matches!(app.handle_key(ctrl('x')), Flow::Continue));
        assert_eq!(app.x(), x0, "Ctrl-X is reserved, not a transform cycle");
    }

    #[test]
    fn slash_filters_feature_names_and_esc_clears() {
        let mut app = app();
        assert_eq!(app.visible.len(), 2);

        app.handle_key(key('/'));
        assert!(app.filter_editing());
        for c in "beta".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert!(!app.filter_editing());
        assert_eq!(app.visible, &[1], "only beta_count matches");
        assert_eq!(
            app.selected_feature(),
            Some(1),
            "selection follows the filter"
        );

        app.handle_key(key('/'));
        app.handle_key(code(KeyCode::Esc));
        assert_eq!(app.visible.len(), 2, "esc clears the filter");
    }

    /// Reopening the filter box clears `self.filter` immediately; `visible`
    /// must refresh in lockstep so the table does not keep showing the
    /// previous (now-stale) filtered set while the help line already shows an
    /// empty query.
    #[test]
    fn reopening_the_filter_refreshes_visible_immediately() {
        let mut app = app();
        app.handle_key(key('/'));
        for c in "beta".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert_eq!(app.visible, &[1], "filtered down to beta_count");

        app.handle_key(key('/'));
        assert_eq!(app.filter(), "", "filter text cleared on reopen");
        assert_eq!(
            app.visible.len(),
            2,
            "visible must refresh to match the cleared filter, not stay stale"
        );
    }

    #[test]
    fn typing_while_filtering_does_not_trigger_commands() {
        let mut app = app();
        let tab0 = app.tab();
        app.handle_key(key('/'));
        app.handle_key(key('l'));
        app.handle_key(key('q'));
        assert_eq!(app.tab(), tab0, "keys are literal text while filtering");
        assert_eq!(app.filter(), "lq");
    }

    #[test]
    fn a_filter_matching_nothing_leaves_no_selection() {
        let mut app = app();
        app.handle_key(key('/'));
        for c in "zzz".chars() {
            app.handle_key(key(c));
        }
        app.handle_key(code(KeyCode::Enter));
        assert!(app.visible.is_empty());
        assert_eq!(app.selected_feature(), None);
        // Must not panic.
        app.handle_key(key('j'));
    }

    #[test]
    fn s_cycles_the_sort_key_and_reorders() {
        let mut app = app();
        let first = app.visible.to_vec();
        // Sort by target mean descending: beta_count (mean 20) before
        // alpha_score (mean 2).
        while app.sort_key() != FeatureColumn::TargetMean {
            app.handle_key(key('s'));
        }
        assert_ne!(app.visible, first.as_slice());
        assert_eq!(app.visible[0], 1);
    }

    /// Gain is a column-aligned array read, so sorting on it must pick up the
    /// row's own number rather than a neighbour's.
    #[test]
    fn sorting_by_gain_reads_the_column_aligned_array() {
        let mut app = app();
        while app.sort_key() != FeatureColumn::Gain {
            app.handle_key(key('s'));
        }
        // The fixture sets gain[j] = j, so descending puts feature 1 first.
        assert_eq!(app.visible, &[1, 0]);
    }

    /// A wholly NaN column gets a NaN AUC — a real, naturally occurring case,
    /// not a synthetic key. It must sort last regardless of direction.
    #[test]
    fn nan_sort_keys_sort_last_in_both_directions() {
        let mut app = App::new(dashboard(
            &["low", "high", "nan_feat"],
            &[
                &[1.0, 2.0, 3.0, 4.0],
                &[10.0, 20.0, 30.0, 40.0],
                &[f64::NAN; 4],
            ],
        ));
        while app.sort_key() != FeatureColumn::Auc {
            app.handle_key(key('s'));
        }
        assert_eq!(
            app.visible.last().copied(),
            Some(2),
            "descending: NaN AUC still sorts last"
        );

        app.handle_key(key('S'));
        assert_eq!(
            app.visible.last().copied(),
            Some(2),
            "ascending: NaN AUC still sorts last"
        );
    }

    /// With no terminal, `run` must take the early-return path and return
    /// `Ok` rather than panicking or failing the search run.
    ///
    /// Skipped under a runner that does give the harness a tty — the point is
    /// the no-terminal path, and actually opening one here would take over the
    /// screen mid-test.
    #[test]
    fn run_returns_ok_and_does_not_panic_when_stdout_is_not_a_terminal() {
        if available() {
            return;
        }
        assert!(
            super::run(dashboard(&["a"], &[&[1.0, 2.0, 3.0, 4.0]])).is_ok(),
            "a non-terminal stdout must skip, not fail the run"
        );
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
    fn catch_panics_passes_through_an_io_error() {
        let result = catch_panics(|| Err(std::io::Error::other("simulated io error")));
        assert!(result.is_err(), "a real Err must not be swallowed");
    }
}
