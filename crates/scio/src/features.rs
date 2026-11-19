//! Features as Seurat names them, and Seurat's `PercentageFeatureSet()` (milestone M15b: the
//! workflow's `percent.mt`, docs/PLAN.md C3).
//!
//! The workflow computes, before SCTransform,
//! `obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "(?i)^mt-")`, on the object's
//! row names:
//!
//! - **Seurat's row names** ([`seurat_names`]): `Read10X()` names the rows by the features file's
//!   column 2 made unique (`make.unique()`, `unique.features = TRUE`), and SeuratObject's
//!   `CreateSeuratObject()` replaces `_` and `|` with `-`. The R loader names them the same way
//!   (`loaders/R/scrustio/R/names.R`);
//! - **the pattern** ([`matches_mt_pattern`]): R's `grep()` with `perl = FALSE` is the TRE engine,
//!   where `(?i)` makes the rest case-insensitive and `^` anchors at the start: a name matches
//!   when it begins with `m` or `M`, then `t` or `T`, then `-` (checked in the pinned env on
//!   names such as `MT-ND1`, `mT-Cytb`, `MT--`, ` mt-a`, `mt_a` and a fullwidth `ｍ`:
//!   `validation/R/regress_vectors.R`). `make.unique()` only appends `.1`, `.2`, ... to a name,
//!   so it never changes whether a name matches; the `_` and `|` replacement does (`Mt_Atp6`
//!   becomes `Mt-Atp6`);
//! - **the percentage** ([`Counts::percent_feature_set`]): Seurat 5.5.1's
//!   `layer.sums / object[[]][cells, "nCount_RNA"] * 100`, where `layer.sums` is
//!   `Matrix::colSums()` of the counts' rows for the features (Matrix 1.7-6's `Csparse_colsum`:
//!   each column's stored values in stored order, added in double) and `nCount_RNA` is
//!   `CreateSeuratObject()`'s `colSums()` of the whole column, the same way. The rows are taken in
//!   the counts' order, as `grep()` returns them (for integer counts the order cannot matter:
//!   every partial sum is exact). A cell with no counts gets `0/0 * 100`, NaN, as in R; a feature
//!   set that matches nothing gives zeros.

use crate::counts::Counts;
use crate::tenx::Feature;
use std::collections::HashSet;

/// The workflow's pattern for mitochondrial features (`PercentageFeatureSet(pattern = )`).
pub const MT_PATTERN: &str = "(?i)^mt-";

/// R's `make.unique(names, sep = ".")` (`do_makeunique` in R 4.5.3's `src/main/unique.c`): a name
/// seen before becomes `name.<cnt>`, trying `cnt` from the first unused value for that name
/// upwards until the result is not among the names so far (the original names and those already
/// made). The first occurrence keeps its name.
pub fn make_unique(names: &[String]) -> Vec<String> {
    let n = names.len();
    let mut out: Vec<String> = names.to_vec();
    if n < 2 {
        return out;
    }
    let mut used: HashSet<String> = names.iter().cloned().collect();
    // index of each name's first occurrence, and the next count to try for it
    let mut first: std::collections::HashMap<&str, usize> = std::collections::HashMap::with_capacity(n);
    let mut cnts = vec![1usize; n];
    // R's one-element buffer `newx`, reused from name to name
    let mut newx = String::new();
    for (i, s) in names.iter().enumerate() {
        let Some(&dp) = first.get(s.as_str()) else {
            first.insert(s.as_str(), i);
            continue;
        };
        // try appending 1, 2, 3, ..., n - 1 until it is not already in use
        let mut cnt = cnts[dp];
        while cnt < n {
            newx = format!("{s}.{cnt}");
            if !used.contains(&newx) {
                break;
            }
            cnt += 1;
        }
        used.insert(newx.clone());
        out[i] = newx.clone();
        cnts[dp] = cnt + 1;
    }
    out
}

/// Seurat's row names for features named `names` (the features file's column 2, in row order):
/// `make.unique()` as `Read10X()` applies it, then `_` and `|` replaced with `-` as
/// `CreateSeuratObject()` does.
pub fn seurat_names(names: &[String]) -> Vec<String> {
    make_unique(names).into_iter().map(|s| s.replace(['_', '|'], "-")).collect()
}

/// Whether a Seurat row name matches the workflow's [`MT_PATTERN`], `"(?i)^mt-"`, as R's
/// `grep()` (TRE) matches it: it begins with `m` or `M`, `t` or `T`, then `-`.
pub fn matches_mt_pattern(name: &str) -> bool {
    let b = name.as_bytes();
    b.len() >= 3 && (b[0] == b'm' || b[0] == b'M') && (b[1] == b't' || b[1] == b'T') && b[2] == b'-'
}

/// The rows (0-based, in the counts' order) whose Seurat name matches [`MT_PATTERN`].
pub fn mt_pattern_rows(features: &[Feature]) -> Vec<usize> {
    let names: Vec<String> = features.iter().map(|f| f.name.clone()).collect();
    seurat_names(&names).iter().enumerate().filter(|(_, s)| matches_mt_pattern(s)).map(|(i, _)| i).collect()
}

/// Which rows a list of feature ids (the features file's column 1) selects: the rows whose id is
/// listed, in the counts' order, and the listed ids found among none of the features, in the
/// list's order (each once).
pub fn id_rows(features: &[Feature], ids: &[String]) -> (Vec<usize>, Vec<String>) {
    let listed: HashSet<&str> = ids.iter().map(String::as_str).collect();
    let rows: Vec<usize> = features.iter().enumerate().filter(|(_, f)| listed.contains(f.id.as_str())).map(|(i, _)| i).collect();
    let present: HashSet<&str> = features.iter().map(|f| f.id.as_str()).collect();
    let mut seen = HashSet::new();
    let missing = ids.iter().filter(|id| !present.contains(id.as_str()) && seen.insert(id.as_str())).cloned().collect();
    (rows, missing)
}

/// A list of feature ids as a file holds it: one id per line, each line's surrounding whitespace
/// (spaces, tabs, a carriage return) ignored, blank lines and lines starting with `#` skipped.
/// Returns the ids in the file's order, repeats included.
pub fn parse_id_list(text: &str) -> Vec<String> {
    text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')).map(str::to_string).collect()
}

impl Counts {
    /// Seurat's row names for these counts ([`seurat_names`] of the features' names).
    pub fn seurat_names(&self) -> Vec<String> {
        let names: Vec<String> = self.features().into_iter().map(|f| f.name).collect();
        seurat_names(&names)
    }

    /// `PercentageFeatureSet()` for the features in `rows` (0-based, ascending): per cell, their
    /// counts summed in stored order, divided by the cell's `nCount_RNA` (all its counts, summed
    /// in stored order), times 100, in R's order of operations.
    pub fn percent_feature_set(&self, rows: &[usize]) -> Vec<f64> {
        let m = &self.matrix;
        let mut set = vec![false; m.nrow];
        for &r in rows {
            set[r] = true;
        }
        (0..m.ncol)
            .map(|j| {
                let (i, x) = m.col(j);
                let mut s = 0.0f64;
                let mut total = 0.0f64;
                for (&r, &v) in i.iter().zip(x) {
                    total += v;
                    if set[r as usize] {
                        s += v;
                    }
                }
                s / total * 100.0
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::CscMatrix;

    fn strs(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn make_unique_as_r() {
        // R: make.unique(c("a", "a", "b", "a", "a.1", "b")): "a.1" is an original name, so the
        // first repeat of "a" becomes "a.2" (tests/fixtures/clustering/percent_mt/make_unique.tsv
        // has R's results for more cases)
        assert_eq!(make_unique(&strs(&["a", "a", "b", "a", "a.1", "b"])), strs(&["a", "a.2", "b", "a.3", "a.1", "b.1"]));
        // R: make.unique(c("a", "a.1", "a")): "a.1" is taken, so the duplicate becomes "a.2"
        assert_eq!(make_unique(&strs(&["a", "a.1", "a"])), strs(&["a", "a.1", "a.2"]));
        assert_eq!(make_unique(&strs(&["x"])), strs(&["x"]));
    }

    #[test]
    fn seurat_names_and_the_pattern() {
        let names = strs(&["mt-Co1", "MT-ND2", "Mt_Atp6", "mT|Cytb", "Mtor", "xmt-a", "mt", "mt", "Gene_1", "Gene_1", "MT-"]);
        let s = seurat_names(&names);
        assert_eq!(s, strs(&["mt-Co1", "MT-ND2", "Mt-Atp6", "mT-Cytb", "Mtor", "xmt-a", "mt", "mt.1", "Gene-1", "Gene-1.1", "MT-"]));
        let hits: Vec<bool> = s.iter().map(|n| matches_mt_pattern(n)).collect();
        assert_eq!(hits, vec![true, true, true, true, false, false, false, false, false, false, true]);
        assert!(!matches_mt_pattern("ｍt-a") && !matches_mt_pattern(" mt-a") && matches_mt_pattern("MT--"));
    }

    #[test]
    fn id_lists() {
        let f: Vec<Feature> = ["A\ta\tG", "B\tb\tG", "C\tc\tG"].iter().map(|l| Feature::parse(l)).collect();
        let ids = parse_id_list("# mt\nC\n\n  A \r\nC\nZ\nZ\n#B\n");
        assert_eq!(ids, strs(&["C", "A", "C", "Z", "Z"]));
        assert_eq!(id_rows(&f, &ids), (vec![0, 2], strs(&["Z"])));
    }

    #[test]
    fn percent_feature_set_as_r() {
        // 3 genes x 3 cells; the third cell has no counts
        let m = CscMatrix { nrow: 3, ncol: 3, p: vec![0, 3, 4, 4], i: vec![0, 1, 2, 1], x: vec![1.0, 2.0, 5.0, 4.0] };
        let c = Counts::new(m, strs(&["c1", "c2", "c3"]), strs(&["g1\tg1\tG", "g2\tg2\tG", "g3\tg3\tG"])).unwrap();
        let p = c.percent_feature_set(&[0, 2]);
        assert_eq!(p[0].to_bits(), (6.0f64 / 8.0 * 100.0).to_bits());
        assert_eq!(p[1], 0.0);
        assert!(p[2].is_nan());
        let none = c.percent_feature_set(&[]);
        assert_eq!((none[0], none[1]), (0.0, 0.0));
        assert!(none[2].is_nan());
    }
}
