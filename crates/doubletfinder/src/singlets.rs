//! The singlets as a new 10x directory (milestone M6b; docs/PLAN.md A3 and B3): the counts
//! DoubletFinder read, without the cells it calls doublets. This is the workflow's
//! DoubletFinder step at the level of the counts: `doubletfinder.R` keeps the singlets with
//! `subset()`, which leaves their stored entries as they are (explicit zeros included), and
//! its later SCTransform and clustering add assays without touching the RNA counts.
//!
//! The run frees the counts when it builds the final run's input, which keeps its memory peak
//! down; [`write_singlets`] reads them again from the input directory once the calls are made.

use crate::run::Output;
use crate::DfError;
use scio::tenx::{self, Written};
use std::path::Path;

/// The singlets' directory in a run's output.
pub const DIR: &str = "singlets";

/// Writes `<out>/singlets/` (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`): the
/// columns of the counts in `counts_dir` that `output` calls singlets, in input order, each
/// exactly as stored; the input's feature lines, as they are. `software` goes in the matrix's
/// `%metadata_json` line. The files do not depend on `threads`.
pub fn write_singlets(counts_dir: &Path, output: &Output, out: &Path, software: &str, threads: usize) -> Result<Written, DfError> {
    let barcodes = tenx::read_barcodes(counts_dir)?;
    if barcodes != output.barcodes {
        return Err(DfError::Input(format!("{}: the barcodes are not the run's", counts_dir.display())));
    }
    let features = tenx::read_lines(&tenx::find(counts_dir, "features.tsv")?)?;
    let m = scio::sparse::read_mtx(counts_dir)?;
    if m.nrow != features.len() || m.ncol != barcodes.len() || m.nrow != output.n_genes || output.doublet.len() != m.ncol {
        return Err(DfError::Input(format!(
            "{}: the matrix is {} x {}, with {} features and {} barcodes; the run had {} genes and {} cells",
            counts_dir.display(),
            m.nrow,
            m.ncol,
            features.len(),
            barcodes.len(),
            output.n_genes,
            output.doublet.len()
        )));
    }
    let cols: Vec<usize> = (0..m.ncol).filter(|&j| !output.doublet[j]).collect();
    let kept: Vec<String> = cols.iter().map(|&j| barcodes[j].clone()).collect();
    drop(barcodes);
    Ok(tenx::write_columns(&out.join(DIR), &m, &cols, &kept, &features, software, threads)?)
}
