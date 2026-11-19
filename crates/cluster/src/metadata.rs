//! The per-cell metadata the clustering's `metadata.tsv` starts from (docs/PLAN.md C3,
//! "Outputs"):
//!
//! - **the counts' own columns**, as `CreateSeuratObject()` computes them (SeuratObject 5.4.0's
//!   `.CalcN.default`): `nCount_RNA = Matrix::colSums(counts)`, which sums each column's stored
//!   values in stored order in double (Matrix 1.7-6's `Csparse_colsum`), and
//!   `nFeature_RNA = Matrix::colSums(counts > 0)`, the stored values above 0 (explicit zeros do
//!   not count); and, when asked for (milestone M15b), `percent.mt`, as the workflow's scripts
//!   compute it with `PercentageFeatureSet()` before SCTransform ([`crate::covariates`]);
//! - **the columns of `--metadata`**, carried over by barcode.
//!
//! The staged route needs the counts' columns too, without the counts: `scrust sctransform`
//! writes them in [`FILE`], and `scrust pca` and `scrust neighbors` carry the file forward
//! ([`carry`]), so `scrust cluster --sct`, `--pca` and `--snn` write the same `metadata.tsv` as
//! `scrust cluster --counts`.

use scio::{fmt17, CscMatrix, Table};
use std::io;
use std::path::Path;

/// The file the counts' columns go to in a stage's output directory.
pub const FILE: &str = "cell_metadata.tsv";

/// `nCount_RNA` and `nFeature_RNA` of every cell, in the counts' column order, and `percent.mt`
/// when it was asked for (M15b).
#[derive(Debug, Clone, PartialEq)]
pub struct RnaMetadata {
    pub n_count: Vec<f64>,
    pub n_feature: Vec<i64>,
    pub percent_mt: Option<Vec<f64>>,
}

impl RnaMetadata {
    /// The columns `CreateSeuratObject()` computes from `counts` (genes x cells).
    pub fn of(counts: &CscMatrix) -> RnaMetadata {
        let mut n_count = Vec::with_capacity(counts.ncol);
        let mut n_feature = Vec::with_capacity(counts.ncol);
        for j in 0..counts.ncol {
            let (_, x) = counts.col(j);
            let mut s = 0.0f64;
            for &v in x {
                s += v;
            }
            n_count.push(s);
            n_feature.push(x.iter().filter(|&&v| v > 0.0).count() as i64);
        }
        RnaMetadata { n_count, n_feature, percent_mt: None }
    }

    /// The column texts as R's reference scripts write them (doubles with 17 significant digits,
    /// integers as they are).
    pub fn texts(&self) -> (Vec<String>, Vec<String>) {
        (self.n_count.iter().map(|&v| fmt17(v)).collect(), self.n_feature.iter().map(|v| v.to_string()).collect())
    }

    /// `percent.mt`'s texts (17 significant digits, `NaN` for a cell with no counts), when computed.
    pub fn percent_mt_texts(&self) -> Option<Vec<String>> {
        self.percent_mt.as_ref().map(|v| v.iter().map(|&x| fmt17(x)).collect())
    }

    /// The columns [`FILE`] holds: 3 (`barcode`, `nCount_RNA`, `nFeature_RNA`), or 4 with `percent.mt`.
    pub fn n_columns(&self) -> usize {
        if self.percent_mt.is_some() {
            4
        } else {
            3
        }
    }

    /// Writes [`FILE`]'s layout: `barcode`, `nCount_RNA`, `nFeature_RNA`, and `percent.mt` when
    /// computed.
    pub fn write(&self, path: &Path, barcodes: &[String]) -> io::Result<()> {
        assert_eq!(barcodes.len(), self.n_count.len());
        let (a, b) = self.texts();
        let pm = self.percent_mt_texts();
        let mut s = String::from(if pm.is_some() { "barcode\tnCount_RNA\tnFeature_RNA\tpercent.mt\n" } else { "barcode\tnCount_RNA\tnFeature_RNA\n" });
        for (i, bc) in barcodes.iter().enumerate() {
            s.push_str(bc);
            s.push('\t');
            s.push_str(&a[i]);
            s.push('\t');
            s.push_str(&b[i]);
            if let Some(p) = &pm {
                s.push('\t');
                s.push_str(&p[i]);
            }
            s.push('\n');
        }
        std::fs::write(path, s).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
    }

    /// Reads [`FILE`]'s layout: the barcodes and the columns (`percent.mt` when the file has it).
    pub fn read(path: &Path) -> Result<(Vec<String>, RnaMetadata), String> {
        let t = Table::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let with_pm = t.header == ["barcode", "nCount_RNA", "nFeature_RNA", "percent.mt"];
        if !with_pm && t.header != ["barcode", "nCount_RNA", "nFeature_RNA"] {
            return Err(format!("{}: expected the columns barcode, nCount_RNA, nFeature_RNA (and percent.mt)", path.display()));
        }
        let mut barcodes = Vec::with_capacity(t.rows.len());
        let mut m = RnaMetadata { n_count: Vec::with_capacity(t.rows.len()), n_feature: Vec::with_capacity(t.rows.len()), percent_mt: with_pm.then(|| Vec::with_capacity(t.rows.len())) };
        for (i, r) in t.rows.iter().enumerate() {
            let bad = || format!("{}: row {}", path.display(), i + 1);
            barcodes.push(r[0].clone());
            m.n_count.push(scio::parse17(r.get(1).ok_or_else(bad)?).ok_or_else(bad)?);
            m.n_feature.push(r.get(2).ok_or_else(bad)?.parse().map_err(|_| bad())?);
            if let Some(p) = m.percent_mt.as_mut() {
                p.push(scio::parse17(r.get(3).ok_or_else(bad)?).ok_or_else(bad)?);
            }
        }
        Ok((barcodes, m))
    }
}

/// Copies [`FILE`] from a stage's input directory to its output directory, when the input has
/// it (an R run's directory has not; then nothing is copied).
pub fn carry(from_dir: &Path, to_dir: &Path) -> io::Result<()> {
    let src = from_dir.join(FILE);
    if src.is_file() {
        std::fs::copy(&src, to_dir.join(FILE))?;
    }
    Ok(())
}

/// Columns carried over from `--metadata` (a TSV with a `barcode` column, else the first
/// column holds the barcodes): for each of `barcodes`, the row with that barcode, or `NA` where
/// the table has none. Returns the column names and, per column, the cells' values.
pub fn extra_columns(path: &Path, barcodes: &[String], taken: &[&str]) -> Result<(Vec<String>, Vec<Vec<String>>), String> {
    let bytes = scio::tenx::read_maybe_gz(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let text = String::from_utf8(bytes).map_err(|_| format!("{}: not text", path.display()))?;
    let t = Table::parse(&text).ok_or_else(|| format!("{}: empty", path.display()))?;
    let bc = t.col("barcode").unwrap_or(0);
    let mut index = std::collections::HashMap::with_capacity(t.rows.len());
    for (i, r) in t.rows.iter().enumerate() {
        let b = r.get(bc).ok_or_else(|| format!("{}: row {} has no barcode", path.display(), i + 1))?;
        if index.insert(b.as_str(), i).is_some() {
            return Err(format!("{}: barcode {b} appears twice", path.display()));
        }
    }
    let mut names = Vec::new();
    let mut cols = Vec::new();
    for (j, h) in t.header.iter().enumerate() {
        if j == bc {
            continue;
        }
        if taken.contains(&h.as_str()) || names.contains(h) {
            return Err(format!("{}: column {h} is one the clustering writes, or appears twice", path.display()));
        }
        if h.is_empty() || h.contains(['\t', '\n']) {
            return Err(format!("{}: column {} has no usable name", path.display(), j + 1));
        }
        names.push(h.clone());
        cols.push(barcodes.iter().map(|b| index.get(b.as_str()).and_then(|&i| t.rows[i].get(j)).cloned().unwrap_or_else(|| "NA".into())).collect());
    }
    Ok((names, cols))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_columns_as_seurat_computes_them() {
        // 3 genes x 2 cells; an explicit zero in cell 1
        let m = CscMatrix { nrow: 3, ncol: 2, p: vec![0, 3, 4], i: vec![0, 1, 2, 1], x: vec![2.0, 0.0, 5.0, 0.5] };
        let r = RnaMetadata::of(&m);
        assert_eq!(r.n_count, vec![7.0, 0.5]);
        assert_eq!(r.n_feature, vec![2, 1]);
        assert_eq!(r.texts().0, vec!["7", "0.5"]);
    }

    #[test]
    fn extra_columns_by_barcode() {
        let dir = std::env::temp_dir().join(format!("scrust-meta-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("m.tsv");
        std::fs::write(&p, "barcode\tpercent.mt\tsample\nB\t1.5\tx\nA\t2\ty\n").unwrap();
        let (names, cols) = extra_columns(&p, &["A".into(), "C".into(), "B".into()], &["nCount_RNA"]).unwrap();
        assert_eq!(names, vec!["percent.mt", "sample"]);
        assert_eq!(cols[0], vec!["2", "NA", "1.5"]);
        std::fs::write(&p, "barcode\tnCount_RNA\nA\t1\n").unwrap();
        assert!(extra_columns(&p, &["A".into()], &["nCount_RNA"]).is_err());
        std::fs::write(&p, "barcode\tx\nA\t1\nA\t2\n").unwrap();
        assert!(extra_columns(&p, &["A".into()], &[]).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
