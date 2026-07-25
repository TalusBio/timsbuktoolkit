//! Thin wrapper over the real `timsseek` apex finder.
//!
//! Runs BOTH scoring stages that the production pipeline performs on an
//! extraction:
//!   * Pass 1 -- quick apex location (`suggest_apex`), the Phase-1 "prescore".
//!   * Pass 2 -- full scoring (`score_at`), the Phase-3 calibrated/secondary
//!     score. Its `joint_apex_cycle` may differ from Pass 1's peak pick.
//!
//! No reimplementation: `TraceScorer`, `ElutionTraces`, `ApexLocation` and
//! `ApexBlocks` are the exact types the CLI pipeline uses.

use timsseek::scoring::apex_finding::{
    ApexBlocks,
    ApexLocation,
    ElutionTraces,
    Extraction,
    TraceScorer,
};
use timsseek::scoring::blocks::apex_evidence::ApexEvidence;

/// The staged pass results for one extraction.
pub struct Passes {
    /// Pass 1: quick apex location (prescore).
    pub pass1: ApexLocation,
    /// The window-global apex evidence Pass 2 was scored against. Returned
    /// rather than kept internal so a caller sweeping cycles reuses this one
    /// instead of recomputing it per cycle.
    pub evidence: ApexEvidence,
    /// Pass 2: full scored apex (may sit at a different cycle than pass1).
    pub pass2: ApexBlocks,
}

/// Everything the UI needs after scoring one synthetic extraction.
pub struct ScoreResult {
    /// The per-cycle intermediate traces (cosine/scribe/lazyscore/...).
    pub traces: ElutionTraces,
    /// Pass 1: quick apex location (prescore).
    pub pass1: ApexLocation,
    /// Pass 2: full scored apex (may sit at a different cycle than pass1).
    pub pass2: ApexBlocks,
    /// `main_score` from Pass 2 evaluated at EVERY cycle -- the score landscape
    /// of this one extraction. Only the 11 apex features vary with cycle: the
    /// apex evidence is window-global and identical at every entry.
    pub landscape: Vec<f32>,
}

/// Run the two staged passes on a CALLER-OWNED scorer, returning the apex
/// results without cloning the traces. This is the allocation-free hot path:
/// reuse one `TraceScorer` across many extractions, exactly like production's
/// `ScoringWorker`. Used by the sensitivity bench.
pub fn run_with(
    scorer: &mut TraceScorer,
    extraction: &Extraction<String>,
    rt_mapper: &dyn Fn(usize) -> u32,
) -> Result<Passes, String> {
    // Stage A: per-cycle traces (resizes the scorer's buffers as needed).
    scorer
        .compute_traces(extraction)
        .map_err(|e| format!("compute_traces failed: {e:?}"))?;

    let cycle_offset = extraction.chromatograms.cycle_offset();

    // Stage B (Pass 1): quick apex suggestion.
    let pass1 = scorer
        .suggest_apex(rt_mapper, cycle_offset)
        .map_err(|e| format!("suggest_apex failed: {e:?}"))?;

    // Stage B' : window-global apex evidence (cycle-invariant).
    let evidence = scorer.compute_apex_evidence(extraction);

    // Stage C (Pass 2): full score at the suggested apex.
    let pass2 = scorer
        .score_at(extraction, &evidence, pass1.apex_cycle, &pass1, rt_mapper)
        .map_err(|e| format!("score_at failed: {e:?}"))?;

    Ok(Passes {
        pass1,
        evidence,
        pass2,
    })
}

/// Run the two staged passes over `extraction`, allocating a fresh scorer and
/// cloning the traces for display. Convenience for the interactive app.
///
/// `rt_mapper` maps a global cycle index to retention time in ms, exactly as
/// the pipeline supplies it.
pub fn run(
    extraction: &Extraction<String>,
    rt_mapper: &dyn Fn(usize) -> u32,
) -> Result<ScoreResult, String> {
    let n_cycles = extraction.chromatograms.num_cycles();
    let max_frags = extraction.chromatograms.fragments.num_ions().max(1);
    let mut scorer = TraceScorer::new(n_cycles, max_frags);
    let Passes {
        pass1,
        evidence,
        pass2,
    } = run_with(&mut scorer, extraction, rt_mapper)?;
    // Traces and the window-global apex evidence are already computed, so each
    // extra `score_at` is just the cycle-local feature block — the sweep is
    // O(cycles * window), not O(cycles^2). Cycles that fail to score give 0.
    let landscape = (0..n_cycles)
        .map(|c| {
            scorer
                .score_at(extraction, &evidence, c, &pass1, rt_mapper)
                .map(|b| b.primary.main_score)
                .unwrap_or(0.0)
        })
        .collect();
    Ok(ScoreResult {
        traces: scorer.traces().clone(),
        pass1,
        pass2,
        landscape,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::{
        SimParams,
        build,
    };

    #[test]
    fn clean_peak_passes_agree() {
        let mut p = SimParams::default();
        p.noise_floor = 0.0;
        p.random_peaks.enabled = false;
        p.apex_cycle = 30.0;
        let data = build(&p);
        let map = p.rt_mapper();
        let res = run(&data.extraction, &map).unwrap();
        // apex_profile peak (pass1) must be near the configured apex.
        assert!((res.pass1.apex_cycle as i64 - 30).abs() <= 2);
        // pass2 joint apex agrees with pass1 on a clean peak.
        assert!((res.pass2.joint_apex_cycle as i64 - res.pass1.apex_cycle as i64).abs() <= 3);
    }

    #[test]
    fn landscape_peaks_at_the_scored_apex() {
        let mut p = SimParams::default();
        p.noise_floor = 0.05;
        p.random_peaks.enabled = false;
        p.apex_cycle = 30.0;
        let data = build(&p);
        let res = run(&data.extraction, &p.rt_mapper()).unwrap();

        assert_eq!(res.landscape.len(), p.n_cycles);
        // The entry at the scored apex IS the reported main score.
        assert_eq!(
            res.landscape[res.pass2.joint_apex_cycle],
            res.pass2.primary.main_score
        );
        // ... and on a clean peak that entry is the landscape maximum.
        let best = res
            .landscape
            .iter()
            .cloned()
            .fold(f32::NEG_INFINITY, f32::max);
        assert_eq!(best, res.pass2.primary.main_score);
    }

    #[test]
    fn empty_transition_lowers_score() {
        // With realistic noise (>0) an absent fragment's row is noise, so it
        // survives filter_zero_intensity_ions and stays expected -> it should
        // score lower than the fully-observed case (coverage/cosine penalty).
        let mut base = SimParams::default();
        base.noise_floor = 0.1;
        base.random_peaks.enabled = false;
        base.apex_cycle = 30.0;

        let clean = build(&base);
        let map = base.rt_mapper();
        let clean_score = run(&clean.extraction, &map)
            .unwrap()
            .pass2
            .primary
            .main_score;

        let mut dropped = base.clone();
        dropped.real_fragments[0].obs_scale = 0.0; // top fragment goes missing
        let d = build(&dropped);
        let dropped_score = run(&d.extraction, &map).unwrap().pass2.primary.main_score;

        assert!(
            dropped_score < clean_score,
            "dropped={dropped_score:e} clean={clean_score:e}"
        );
    }
}
