//! What `scrust doubletfinder` reads (docs/PLAN.md A3 and B3): a 10x Matrix Market directory
//! (`matrix.mtx(.gz)`, `barcodes.tsv(.gz)`, `features.tsv(.gz)`) and the cells' cluster labels.
//!
//! The workflow's DoubletFinder reads its Seurat object's RNA counts, the cell names and
//! `seurat_clusters`, a factor whose level order the homotypic adjustment depends on
//! (`modelHomotypic` tabulates it). The export writes the labels as a TSV (`barcode`,
//! `seurat_clusters`, in the matrix's column order) and the levels, in order, one per line
//! (`seurat_clusters_levels.txt`).

use crate::run::Input;
use crate::DfError;
use scio::{tenx, CscMatrix, Table};
use std::path::{Path, PathBuf};

/// Where the factor's levels came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LevelsSource {
    File(PathBuf),
    /// No levels file: the order of first appearance.
    FirstAppearance,
}

impl LevelsSource {
    pub fn describe(&self) -> String {
        match self {
            LevelsSource::File(p) => p.display().to_string(),
            LevelsSource::FirstAppearance => "order of first appearance".into(),
        }
    }
}

/// The cluster labels of `barcodes` from a TSV with a header: the `barcode` column (else the
/// first) and the `seurat_clusters` column (else the second). The barcodes must be the
/// matrix's, in its order.
pub fn read_clusters(path: &Path, barcodes: &[String]) -> Result<Vec<String>, DfError> {
    let t = Table::read(path).map_err(|e| DfError::Input(format!("{}: {e}", path.display())))?;
    if t.header.len() < 2 {
        return Err(DfError::Input(format!("{}: expected a barcode and a cluster column", path.display())));
    }
    let bc = t.col("barcode").unwrap_or(0);
    let cl = t.col("seurat_clusters").unwrap_or(1);
    if t.rows.len() != barcodes.len() || t.rows.iter().zip(barcodes).any(|(r, b)| r.get(bc) != Some(b)) {
        return Err(DfError::Input(format!("{}: barcodes differ from the matrix's (order or set)", path.display())));
    }
    Ok(t.rows.iter().map(|r| r.get(cl).cloned().unwrap_or_default()).collect())
}

/// The factor's levels: from `levels` if given, else `seurat_clusters_levels.txt` next to the
/// clusters file if there is one, else the labels' order of first appearance. (A label that is
/// not a level is NA in R's factor, which `modelHomotypic`'s `table()` leaves out.)
pub fn read_levels(levels: Option<&Path>, clusters_path: &Path, labels: &[String]) -> Result<(Vec<String>, LevelsSource), DfError> {
    let beside = clusters_path.parent().map(|d| d.join("seurat_clusters_levels.txt")).filter(|p| p.is_file());
    let (lv, source) = match levels.map(Path::to_path_buf).or(beside) {
        Some(p) => {
            let text = std::fs::read_to_string(&p).map_err(|e| DfError::Input(format!("{}: {e}", p.display())))?;
            (text.lines().map(str::to_owned).collect::<Vec<String>>(), LevelsSource::File(p))
        }
        None => {
            let mut lv: Vec<String> = Vec::new();
            for l in labels {
                if !lv.contains(l) {
                    lv.push(l.clone());
                }
            }
            (lv, LevelsSource::FirstAppearance)
        }
    };
    Ok((lv, source))
}

/// The cells, their labels and levels, and the number of genes, from a 10x directory and the
/// clusters file.
pub fn read_input(counts_dir: &Path, clusters: &Path, levels: Option<&Path>) -> Result<(Input, LevelsSource), DfError> {
    let barcodes = tenx::read_barcodes(counts_dir).map_err(DfError::Io)?;
    let n_genes = tenx::read_features(counts_dir).map_err(DfError::Io)?.len();
    let labels = read_clusters(clusters, &barcodes)?;
    let (cluster_levels, source) = read_levels(levels, clusters, &labels)?;
    Ok((Input { barcodes, clusters: labels, cluster_levels, n_genes }, source))
}

/// The counts matrix (genes x cells, every stored entry kept), checked against the barcodes and
/// the number of features.
pub fn read_counts(counts_dir: &Path, input: &Input) -> Result<CscMatrix, DfError> {
    let m = scio::sparse::read_mtx(counts_dir).map_err(DfError::Io)?;
    if m.nrow != input.n_genes || m.ncol != input.barcodes.len() {
        return Err(DfError::Input(format!(
            "{}: the matrix is {} x {}, with {} features and {} barcodes",
            counts_dir.display(),
            m.nrow,
            m.ncol,
            input.n_genes,
            input.barcodes.len()
        )));
    }
    Ok(m)
}
