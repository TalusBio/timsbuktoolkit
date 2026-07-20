use crate::IonAnnot;
use crate::models::capabilities::LibCapabilities;

#[derive(Debug, Clone)]
pub struct ModDefinition {
    pub token: String, // verbatim, e.g. "[UNIMOD:4]"
    pub mono_delta: f64,
    pub cs_delta: (i16, i16),
}

#[derive(Debug, Clone)]
pub struct QueryCollection {
    pub caps: LibCapabilities,
    // per-target scalars, len = n_targets (no `id` column - library_id = target index)
    pub precursor_mz: Vec<f64>,
    pub charge: Vec<u8>,
    pub rt_seconds: Vec<f32>,
    pub mobility: Vec<f32>,
    // CSR prefix offsets (n+1)
    pub frag_off: Vec<u32>,
    pub seq_strip_off: Vec<u32>,
    pub seq_mod_off: Vec<u32>,
    pub mod_off: Vec<u32>,
    // fragment arenas (len = total fragments)
    pub frag_labels: Vec<IonAnnot>,
    pub frag_mzs: Vec<f64>,
    // sequences
    pub seq_strip_blob: String,
    pub seq_mod_blob: String,
    // structured mods
    pub mods: Vec<(u8, u16)>,
    pub mod_registry: Vec<ModDefinition>,
}

impl QueryCollection {
    pub fn with_capabilities(caps: LibCapabilities) -> Self {
        Self {
            caps,
            precursor_mz: Vec::new(),
            charge: Vec::new(),
            rt_seconds: Vec::new(),
            mobility: Vec::new(),
            frag_off: vec![0],
            seq_strip_off: vec![0],
            seq_mod_off: vec![0],
            mod_off: vec![0],
            frag_labels: Vec::new(),
            frag_mzs: Vec::new(),
            seq_strip_blob: String::new(),
            seq_mod_blob: String::new(),
            mods: Vec::new(),
            mod_registry: Vec::new(),
        }
    }

    /// Append one target. `mods` are (position, registry_idx) pairs; the caller
    /// is responsible for having registered the mod tokens in `mod_registry`.
    pub fn push_target(
        &mut self,
        precursor_mz: f64,
        charge: u8,
        rt_seconds: f32,
        mobility: f32,
        frags: &[(IonAnnot, f64)],
        seq_strip: &str,
        seq_mod: &str,
        mods: &[(u8, u16)],
    ) {
        self.precursor_mz.push(precursor_mz);
        self.charge.push(charge);
        self.rt_seconds.push(rt_seconds);
        self.mobility.push(mobility);
        for (lab, mz) in frags {
            self.frag_labels.push(*lab);
            self.frag_mzs.push(*mz);
        }
        // CSR offsets are u32 (documented ~40x headroom over real arena sizes).
        // If a monster library ever exceeds that, fail loud instead of wrapping
        // an offset and silently corrupting every range past the overflow.
        self.frag_off.push(
            u32::try_from(self.frag_labels.len()).expect("fragment arena exceeds u32 offset range"),
        );
        self.seq_strip_blob.push_str(seq_strip);
        self.seq_strip_off.push(
            u32::try_from(self.seq_strip_blob.len())
                .expect("stripped-seq blob exceeds u32 offset range"),
        );
        self.seq_mod_blob.push_str(seq_mod);
        self.seq_mod_off.push(
            u32::try_from(self.seq_mod_blob.len())
                .expect("modified-seq blob exceeds u32 offset range"),
        );
        self.mods.extend_from_slice(mods);
        self.mod_off
            .push(u32::try_from(self.mods.len()).expect("mods arena exceeds u32 offset range"));
    }

    pub fn n_targets(&self) -> usize {
        self.charge.len()
    }

    pub fn frag_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.frag_off[tgt] as usize..self.frag_off[tgt + 1] as usize
    }

    pub fn seq_strip_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.seq_strip_off[tgt] as usize..self.seq_strip_off[tgt + 1] as usize
    }

    pub fn seq_mod_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.seq_mod_off[tgt] as usize..self.seq_mod_off[tgt + 1] as usize
    }

    pub fn mod_range(&self, tgt: usize) -> std::ops::Range<usize> {
        self.mod_off[tgt] as usize..self.mod_off[tgt + 1] as usize
    }

    /// Seal after build: release excess capacity on every arena.
    pub fn seal(&mut self) {
        self.precursor_mz.shrink_to_fit();
        self.charge.shrink_to_fit();
        self.rt_seconds.shrink_to_fit();
        self.mobility.shrink_to_fit();
        self.frag_off.shrink_to_fit();
        self.seq_strip_off.shrink_to_fit();
        self.seq_mod_off.shrink_to_fit();
        self.mod_off.shrink_to_fit();
        self.frag_labels.shrink_to_fit();
        self.frag_mzs.shrink_to_fit();
        self.seq_strip_blob.shrink_to_fit();
        self.seq_mod_blob.shrink_to_fit();
        self.mods.shrink_to_fit();
        self.mod_registry.shrink_to_fit();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IonAnnot;

    #[test]
    fn csr_ranges_recover_per_target_fragments() {
        let mut c = QueryCollection::with_capabilities(LibCapabilities::default_diann());
        // target 0: 2 frags; target 1: 3 frags
        c.push_target(
            // precursor_mz
            500.0,
            // charge
            2,
            // rt
            1.0,
            // mob
            0.8,
            // frags
            &[
                (IonAnnot::try_from("y3").unwrap(), 300.0),
                (IonAnnot::try_from("y4").unwrap(), 400.0),
            ],
            // strip
            "PEPTIDEK",
            // modified
            "PEPTIDEK",
            // mods
            &[],
        );
        c.push_target(
            600.0,
            3,
            2.0,
            0.9,
            &[
                (IonAnnot::try_from("y2").unwrap(), 200.0),
                (IonAnnot::try_from("y5").unwrap(), 500.0),
                (IonAnnot::try_from("y6").unwrap(), 600.0),
            ],
            "SAMPLERK",
            "SAMPLERK",
            &[],
        );
        c.seal();
        assert_eq!(c.n_targets(), 2);
        assert_eq!(c.frag_range(0), 0..2);
        assert_eq!(c.frag_range(1), 2..5);
        assert_eq!(c.frag_mzs[c.frag_range(1)][2], 600.0);
        assert_eq!(&c.seq_strip_blob[c.seq_strip_range(0)], "PEPTIDEK");
        assert_eq!(&c.seq_strip_blob[c.seq_strip_range(1)], "SAMPLERK");
    }
}
