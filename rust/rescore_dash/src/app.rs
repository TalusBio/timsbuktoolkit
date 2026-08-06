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
pub(crate) enum Tab {
    Overview,
    Fdr,
    Calibration,
    Features,
}

impl Tab {
    pub(crate) const ALL: [Tab; 4] = [Self::Overview, Self::Fdr, Self::Calibration, Self::Features];

    pub(crate) fn title(self) -> &'static str {
        match self {
            Self::Overview => "Overview",
            Self::Fdr => "FDR",
            Self::Calibration => "Calibration",
            Self::Features => "Features",
        }
    }

    pub(crate) fn index(self) -> usize {
        cycle::index_of(&Self::ALL, self)
    }

    fn shift(self, delta: isize) -> Self {
        cycle::step(&Self::ALL, self, delta)
    }
}

pub(crate) enum Flow {
    Continue,
    Quit,
}

pub(crate) struct App {
    pub(crate) dash: Dashboard,
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
    /// Scroll state for the features table, kept on `App` so `offset` survives
    /// a redraw. Its `selected` is written from `cursor` every frame and must
    /// never be read back; only `offset()` is.
    pub(crate) table_state: TableState,
}

impl App {
    pub(crate) fn new(dash: Dashboard) -> Self {
        let mut app = Self {
            dash,
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

    pub(crate) fn tab(&self) -> Tab {
        self.tab
    }

    pub(crate) fn x(&self) -> Axis {
        self.x
    }

    pub(crate) fn y(&self) -> YTransform {
        self.y
    }

    pub(crate) fn clip(&self) -> bool {
        self.clip
    }

    /// Which of the dashboard's FDR zoom levels the curve panel draws.
    pub(crate) fn q_zoom(&self) -> usize {
        self.q_zoom
    }

    /// The sort key and the order it puts rows in, e.g. `"AUC high-low"`.
    ///
    /// Spelled out rather than shown as an arrow because `sort_desc` reads
    /// differently per key: largest-first for a number, A-first for the name.
    pub(crate) fn sort_summary(&self) -> String {
        let order = match (self.sort, self.sort_desc) {
            (FeatureColumn::Name, true) => "a-z",
            (FeatureColumn::Name, false) => "z-a",
            (_, true) => "high-low",
            (_, false) => "low-high",
        };
        format!("{} {order}", self.sort.label())
    }

    pub(crate) fn filter(&self) -> &str {
        &self.filter
    }

    pub(crate) fn filter_editing(&self) -> bool {
        self.filter_editing
    }

    /// Feature index under the cursor, or `None` when the filter matches
    /// nothing.
    pub(crate) fn selected_feature(&self) -> Option<usize> {
        self.visible.get(self.cursor).copied()
    }

    pub(crate) fn handle_key(&mut self, key: KeyEvent) -> Flow {
        // Handled before anything else — including the filter-editing dispatch
        // below — so Ctrl-C always quits rather than being typed into the filter
        // box or toggling clip.
        if is_ctrl_c(key) {
            return Flow::Quit;
        }
        // CONTROL/ALT combinations are reserved (Ctrl-C above, and headroom
        // for future bindings); only a bare or Shifted character triggers a
        // command, so e.g. Ctrl-X does not cycle the x-transform.
        let plain = !key
            .modifiers
            .intersects(KeyModifiers::CONTROL | KeyModifiers::ALT);
        if self.filter_editing {
            self.handle_filter_key(key, plain);
            return Flow::Continue;
        }
        match key.code {
            KeyCode::Char('q') if plain => return Flow::Quit,
            KeyCode::Char('l') | KeyCode::Right | KeyCode::Tab if plain => {
                self.tab = self.tab.shift(1)
            }
            KeyCode::Char('h') | KeyCode::Left | KeyCode::BackTab if plain => {
                self.tab = self.tab.shift(-1)
            }
            KeyCode::Char('j') | KeyCode::Down if plain => {
                self.cursor += 1;
                self.clamp_cursor();
            }
            KeyCode::Char('k') | KeyCode::Up if plain => {
                self.cursor = self.cursor.saturating_sub(1)
            }
            KeyCode::Char('x') if plain => self.x = self.x.next(),
            KeyCode::Char('X') if plain => self.x = self.x.prev(),
            KeyCode::Char('y') if plain => self.y = self.y.next(),
            KeyCode::Char('Y') if plain => self.y = self.y.prev(),
            KeyCode::Char('c') if plain => self.clip = !self.clip,
            KeyCode::Char('z') if plain => self.cycle_q_zoom(1),
            KeyCode::Char('Z') if plain => self.cycle_q_zoom(-1),
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
                self.refresh_visible();
            }
            _ => {}
        }
        Flow::Continue
    }

    fn handle_filter_key(&mut self, key: KeyEvent, plain: bool) {
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
                self.refresh_visible();
            }
            KeyCode::Char(c) if plain => {
                self.filter.push(c);
                self.refresh_visible();
            }
            _ => {}
        }
    }

    /// Move `q_zoom` `delta` levels through the dashboard's zoom list, wrapping.
    fn cycle_q_zoom(&mut self, delta: isize) {
        let n = self.dash.n_q_zooms() as isize;
        self.q_zoom = (self.q_zoom as isize + delta).rem_euclid(n) as usize;
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
            let value = |j: usize| self.dash.feature_value(j, self.sort).unwrap_or(f64::NAN);
            idx.sort_by(|&a, &b| cmp_nan_last(value(a), value(b), desc));
        }

        self.visible = idx;
        self.cursor = previous
            .and_then(|j| self.visible.iter().position(|&v| v == j))
            .unwrap_or(0);
        self.clamp_cursor();
    }

    fn clamp_cursor(&mut self) {
        self.cursor = self.cursor.min(self.visible.len().saturating_sub(1));
    }
}

/// Whether a dashboard could be shown at all.
fn available() -> bool {
    use std::io::IsTerminal;
    std::io::stdout().is_terminal()
}

/// Open the dashboard and block until the user quits.
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

/// `Error = std::io::Error` pins the backend to one whose draw failures are
/// already `io::Error`, matching `event::read`. `ratatui::try_init` returns
/// such a backend.
fn event_loop<B: ratatui::backend::Backend<Error = std::io::Error>>(
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
    use crate::precompute::tests::dashboard;
    use ratatui::crossterm::event::{
        KeyCode,
        KeyEvent,
        KeyModifiers,
    };

    pub(crate) fn key(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE)
    }

    pub(crate) fn code(c: KeyCode) -> KeyEvent {
        KeyEvent::new(c, KeyModifiers::NONE)
    }

    fn ctrl(c: char) -> KeyEvent {
        KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL)
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
    fn c_toggles_clipping() {
        let mut app = app();
        assert!(app.clip(), "percentile clip is the default");
        app.handle_key(key('c'));
        assert!(!app.clip());
        app.handle_key(key('c'));
        assert!(app.clip());
    }

    #[test]
    fn z_cycles_the_fdr_zoom_in_both_directions() {
        let mut app = app();
        let n = app.dash.n_q_zooms();
        assert!(n > 1, "a single zoom level makes the key pointless");

        let start = app.q_zoom();
        for _ in 0..n {
            app.handle_key(key('z'));
        }
        assert_eq!(app.q_zoom(), start, "z must wrap, not run off the end");

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
        assert_eq!(app.filter(), "", "reopening clears the filter text");
        assert_eq!(app.visible.len(), 2, "and refreshes visible with it");
        app.handle_key(code(KeyCode::Esc));
        assert_eq!(app.visible.len(), 2, "esc clears the filter");
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

    /// Gain is a column-aligned array read, so sorting on it must pick up the
    /// row's own number rather than a neighbour's.
    #[test]
    fn sorting_by_gain_reads_the_column_aligned_array() {
        let mut app = app();
        while !app.sort_summary().starts_with("gain ") {
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
        while !app.sort_summary().starts_with("AUC ") {
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

    /// A panic must not fail the run, and a real `Err` must not be swallowed by
    /// the same guard.
    #[test]
    fn catch_panics_converts_a_panic_into_ok_but_passes_errors_through() {
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {})); // silence the default panic printout
        let caught = catch_panics(|| panic!("simulated event-loop panic"));
        std::panic::set_hook(prev_hook);
        assert!(caught.is_ok(), "a caught panic must not fail the run");

        let passed = catch_panics(|| Err(std::io::Error::other("simulated io error")));
        assert!(passed.is_err(), "a real Err must not be swallowed");
    }
}
