//! The one set of modification ontologies this process builds.
//!
//! Two crates need it and neither owns it: the mzSpecLib reader parses each
//! analyte's ProForma through mzannotate, which takes an `&Ontologies`, and
//! timsseek's fallback sequence parser needs the same indexes. It lives here
//! because timsseek depends on timsquery and not the reverse, so this is the
//! lowest crate both can reach.
//!
//! Build the indexes anywhere else and the process holds two copies.
//! `ast-grep --lang rust -p 'mzcv::CVIndex::init_static()'` should report only
//! this file.

use std::sync::OnceLock;

/// The modification ontologies, built once on first use.
///
/// GNOme is omitted. It adds 191,529 entries and 26.4 MB, taking
/// initialization from about 48 ms to 2.6 s, and a GNO modification already
/// produces no usable parse downstream: `count_carbon_sulphur_in_sequence`
/// rejects it during formula counting. PSI-MOD, XL-MOD, Unimod and RESID are
/// loaded because formula counts need them.
///
/// Note this is PSI-**MOD**, the protein-modification vocabulary, and not
/// PSI-**MS**. Nothing here carries `MS:` terms, which is why the reader's
/// `spectrum origin type` hierarchy is vendored separately in
/// [`serde::psims_origin_type`](crate::serde). The names are close enough to
/// mislead.
///
/// Lazily built, and which libraries pay for it differs by format. A DIA-NN
/// library whose sequences all match timsseek's fast byte-walk parser never
/// reaches here. An mzSpecLib library always does: mzannotate takes
/// `&Ontologies` to parse an analyte at all.
pub fn ontologies() -> &'static mzcore::ontology::Ontologies {
    static ONTOLOGIES: OnceLock<mzcore::ontology::Ontologies> = OnceLock::new();
    ONTOLOGIES.get_or_init(|| {
        let mut ontologies = mzcore::ontology::Ontologies::empty();
        *ontologies.unimod_mut() = mzcv::CVIndex::init_static();
        *ontologies.psimod_mut() = mzcv::CVIndex::init_static();
        *ontologies.xlmod_mut() = mzcv::CVIndex::init_static();
        *ontologies.resid_mut() = mzcv::CVIndex::init_static();
        ontologies
    })
}
