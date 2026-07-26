//! Dashboard state: which tab, which feature, which transforms, and the key
//! bindings that move between them. Rendering lives in [`crate::ui`].

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
};
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
    pub fn summary(&mut self, j: usize) -> FeatureSummary {
        if let Some(s) = self.cache.get(&j) {
            return *s;
        }
        let column: Vec<f64> = self.view.feature_column(j).collect();
        let s = stats::summarize(column.iter().copied(), &self.view.is_target);
        self.cache.insert(j, s);
        s
    }

    pub fn handle_key(&mut self, key: KeyEvent) -> Flow {
        if self.filter_editing {
            self.handle_filter_key(key);
            return Flow::Continue;
        }
        match key.code {
            KeyCode::Char('q') | KeyCode::Esc => return Flow::Quit,
            KeyCode::Char('l') | KeyCode::Right | KeyCode::Tab => self.tab = self.tab.shift(1),
            KeyCode::Char('h') | KeyCode::Left | KeyCode::BackTab => self.tab = self.tab.shift(-1),
            KeyCode::Char('j') | KeyCode::Down => {
                self.cursor = (self.cursor + 1).min(self.visible.len().saturating_sub(1));
            }
            KeyCode::Char('k') | KeyCode::Up => self.cursor = self.cursor.saturating_sub(1),
            KeyCode::Char('x') => self.x = self.x.next(),
            KeyCode::Char('X') => self.x = self.x.prev(),
            KeyCode::Char('y') => self.y = self.y.next(),
            KeyCode::Char('Y') => self.y = self.y.prev(),
            KeyCode::Char('c') => self.clip = !self.clip,
            KeyCode::Char('s') => {
                self.sort = self.sort.next();
                self.refresh_visible();
            }
            KeyCode::Char('S') => {
                self.sort_desc = !self.sort_desc;
                self.refresh_visible();
            }
            KeyCode::Char('/') => {
                self.filter_editing = true;
                self.filter.clear();
            }
            _ => {}
        }
        Flow::Continue
    }

    fn handle_filter_key(&mut self, key: KeyEvent) {
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
            KeyCode::Char(c) => self.filter.push(c),
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
}
