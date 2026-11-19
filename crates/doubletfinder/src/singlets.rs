//! The singlets as a new 10x directory (milestone M6b; docs/PLAN.md A3 and B3): the counts
//! DoubletFinder read, without the cells it calls doublets. This is the workflow's
//! DoubletFinder step at the level of the counts: `doubletfinder.R` keeps the singlets with
//! `subset()`, which leaves their stored entries as they are (explicit zeros included), and
//! its later SCTransform and clustering add assays without touching the RNA counts.
//!
//! The run frees the counts when it builds the final run's input, which keeps its memory peak
//! down; the command reads them again once the calls are made. [`singlets`] keeps the singlets'
//! columns in memory ([`scio::Counts`], the type tools hand each other: `scrust doubletfinder
//! --recluster` clusters them without reading `singlets/` back, M15), and [`write_singlets`]
//! reads, keeps and writes in one call.

use crate::run::Output;
use crate::DfError;
use scio::tenx::Written;
use scio::Counts;
use std::path::Path;

/// The singlets' directory in a run's output.
pub const DIR: &str = "singlets";

/// The counts of the cells `output` calls singlets, in input order, each column exactly as stored;
/// the feature lines as they are. `counts` must be the run's input (its cells and genes).
pub fn singlets(counts: Counts, output: &Output) -> Result<Counts, DfError> {
    if counts.barcodes != output.barcodes {
        return Err(DfError::Input("the counts' barcodes are not the run's".into()));
    }
    if counts.n_features() != output.n_genes || output.doublet.len() != counts.n_cells() {
        return Err(DfError::Input(format!(
            "the matrix is {} x {}; the run had {} genes and {} cells",
            counts.n_features(),
            counts.n_cells(),
            output.n_genes,
            output.doublet.len()
        )));
    }
    let keep: Vec<bool> = output.doublet.iter().map(|&d| !d).collect();
    Ok(counts.keep_cells(&keep))
}

/// Writes the singlets as `<out>/singlets/` (`matrix.mtx.gz`, `barcodes.tsv.gz`,
/// `features.tsv.gz`); `software` goes in the matrix's `%metadata_json` line. The files do not
/// depend on `threads`.
pub fn write_counts(singlets: &Counts, out: &Path, software: &str, threads: usize) -> Result<Written, DfError> {
    Ok(singlets.write(&out.join(DIR), software, threads)?)
}

/// Reads the counts in `counts_dir` again, keeps the singlets ([`singlets`]) and writes them
/// ([`write_counts`]).
pub fn write_singlets(counts_dir: &Path, output: &Output, out: &Path, software: &str, threads: usize) -> Result<Written, DfError> {
    let counts = Counts::read(counts_dir)?;
    let s = singlets(counts, output).map_err(|e| DfError::Input(format!("{}: {e}", counts_dir.display())))?;
    write_counts(&s, out, software, threads)
}
