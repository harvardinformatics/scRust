//! Cluster labels with their levels: R's factor, as one tool hands it to the next in memory
//! (docs/PLAN.md A3, "In-memory hand-offs"). The clustering chain returns `seurat_clusters` as
//! one; DoubletFinder's homotypic adjustment reads one (milestone M15, clustering on the fly).
//!
//! On disk the same labels are `seurat_clusters.tsv` (barcode and label, in the cells' order)
//! and `seurat_clusters_levels.txt` (the levels in order, one per line): the layout
//! `scrust cluster` writes and `scrust doubletfinder --clusters` reads.

use std::fmt;

/// A factor: one label per cell, as codes into its levels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Factor {
    /// The column's name (`seurat_clusters`, `SCT_snn_res.0.8`).
    pub name: String,
    /// The levels, in the factor's order.
    pub levels: Vec<String>,
    /// 0-based codes into `levels`, one per cell.
    pub codes: Vec<u32>,
}

/// A label that is not one of the factor's levels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NotALevel {
    pub cell: usize,
    pub label: String,
}

impl fmt::Display for NotALevel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "cell {}: label '{}' is not one of the levels", self.cell + 1, self.label)
    }
}

impl std::error::Error for NotALevel {}

impl Factor {
    /// The factor of `labels` with the given levels; every label must be one of them.
    pub fn from_labels(name: &str, labels: &[String], levels: Vec<String>) -> Result<Factor, NotALevel> {
        // a level that repeats maps to its first position, as R's match() does
        let mut index: std::collections::HashMap<&str, u32> = std::collections::HashMap::with_capacity(levels.len());
        for (k, l) in levels.iter().enumerate() {
            index.entry(l.as_str()).or_insert(k as u32);
        }
        let mut codes = Vec::with_capacity(labels.len());
        for (cell, l) in labels.iter().enumerate() {
            match index.get(l.as_str()) {
                Some(&c) => codes.push(c),
                None => return Err(NotALevel { cell, label: l.clone() }),
            }
        }
        Ok(Factor { name: name.to_string(), levels, codes })
    }

    /// The number of cells.
    pub fn len(&self) -> usize {
        self.codes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    /// The cells' labels as text.
    pub fn labels(&self) -> Vec<&str> {
        self.codes.iter().map(|&c| self.levels[c as usize].as_str()).collect()
    }

    /// The cells' labels as owned strings.
    pub fn label_strings(&self) -> Vec<String> {
        self.codes.iter().map(|&c| self.levels[c as usize].clone()).collect()
    }

    /// `seurat_clusters.tsv`'s text: a header (`barcode`, the factor's name), then each cell's
    /// barcode and label, in the cells' order.
    pub fn table_text(&self, barcodes: &[String]) -> String {
        assert_eq!(barcodes.len(), self.codes.len(), "{} barcodes for {} labels", barcodes.len(), self.codes.len());
        let mut s = String::with_capacity(barcodes.len() * 24);
        s.push_str("barcode\t");
        s.push_str(&self.name);
        s.push('\n');
        for (b, &c) in barcodes.iter().zip(&self.codes) {
            s.push_str(b);
            s.push('\t');
            s.push_str(&self.levels[c as usize]);
            s.push('\n');
        }
        s
    }

    /// The levels file's text: one level per line, as R's `writeLines` writes them.
    pub fn levels_text(&self) -> String {
        let mut s = String::new();
        for l in &self.levels {
            s.push_str(l);
            s.push('\n');
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_codes_and_text() {
        let labels: Vec<String> = ["1", "0", "10", "1"].iter().map(|s| s.to_string()).collect();
        let lv: Vec<String> = ["0", "1", "10"].iter().map(|s| s.to_string()).collect();
        let f = Factor::from_labels("seurat_clusters", &labels, lv).unwrap();
        assert_eq!(f.codes, vec![1, 0, 2, 1]);
        assert_eq!(f.labels(), vec!["1", "0", "10", "1"]);
        assert_eq!(f.label_strings(), labels);
        let bc: Vec<String> = ["a", "b", "c", "d"].iter().map(|s| s.to_string()).collect();
        assert_eq!(f.table_text(&bc), "barcode\tseurat_clusters\na\t1\nb\t0\nc\t10\nd\t1\n");
        assert_eq!(f.levels_text(), "0\n1\n10\n");
        let e = Factor::from_labels("x", &labels, vec!["0".into(), "1".into()]).unwrap_err();
        assert_eq!((e.cell, e.label.as_str()), (2, "10"));
    }
}
