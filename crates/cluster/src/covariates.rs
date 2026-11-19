//! The per-cell covariates of SCTransform's `vars.to.regress`, and the workflow's `percent.mt`
//! (milestone M15b; docs/PLAN.md C3, "vars.to.regress and percent.mt").
//!
//! **percent.mt** is `PercentageFeatureSet()` ([`scio::Counts::percent_feature_set`]) of the
//! mitochondrial features, which are either the workflow's ([`MtGenes::Pattern`]: `"(?i)^mt-"` on
//! Seurat's row names, as every workflow script computes it before SCTransform) or the features
//! whose id (the features file's column 1) a list names ([`MtGenes::Ids`], `--mt-genes`). Listed
//! ids that no feature has are reported and skipped; a list that matches no feature is an error.
//! A pattern that matches no feature gives zeros, with a warning, and is refused as a covariate
//! (`--vars-to-regress percent.mt`), naming the remedy (`--mt-genes`).
//!
//! **The covariates** are named as Seurat names them (`vars.to.regress`): `percent.mt`,
//! `nCount_RNA`, `nFeature_RNA` (the counts' own columns, as `CreateSeuratObject()` computes
//! them), or a numeric column of `--metadata` ([`covariate_values`]). Seurat builds a formula from
//! the names (`RegressOutMatrix`: `paste("GENE ~", paste(vars, collapse = "+"))`) and takes
//! `ScaleData`'s rows named like a covariate as covariates too, so the port refuses what R would
//! not regress as written: an unknown name, a repeated one, a name that is not syntactic in R
//! (`make.names(x) != x`, reserved words included), `GENE` (the formula's response), and a name
//! that is also a feature's. Values are read as R's `read.table()` reads them (`R_strtod`,
//! [`rcompat::as_numeric`]), so the covariates are the bits R regresses; a value that is missing,
//! not a number or not finite is an error (R drops the cell from `lm()`'s model frame, or stops,
//! and `qr.resid()` then fails).
//!
//! The values reach the chain in memory ([`CellCovariates`]); the commands read the files.

use crate::metadata::RnaMetadata;
use scio::features::{id_rows, mt_pattern_rows, parse_id_list, MT_PATTERN};
use scio::{Counts, Json};
use std::collections::HashSet;

/// The covariate `PercentageFeatureSet()` makes.
pub const PERCENT_MT: &str = "percent.mt";

/// The covariates the counts themselves give.
pub const FROM_COUNTS: [&str; 3] = [PERCENT_MT, "nCount_RNA", "nFeature_RNA"];

/// How many of the mitochondrial features `summary.json` lists by name, and how many of the ids
/// not found.
pub const FEATURES_LISTED: usize = 50;
pub const NOT_FOUND_LISTED: usize = 10;

/// Which features `percent.mt` counts.
#[derive(Debug, Clone, PartialEq)]
pub enum MtGenes {
    /// The workflow's: Seurat's row names that match `"(?i)^mt-"`.
    Pattern,
    /// The features whose id (column 1 of `features.tsv`) is listed: `--mt-genes`' ids, in the
    /// file's order (repeats included), with the file's path and SHA-256 for the record.
    Ids { ids: Vec<String>, file: String, sha256: String },
}

impl MtGenes {
    /// `--mt-genes`' file: one id per line; blank lines and lines starting with `#` are skipped,
    /// and each line's surrounding whitespace is ignored ([`scio::features::parse_id_list`]).
    pub fn from_file_text(text: &str, file: &str, sha256: &str) -> MtGenes {
        MtGenes::Ids { ids: parse_id_list(text), file: file.to_string(), sha256: sha256.to_string() }
    }
}

/// `percent.mt`, with what defined it.
#[derive(Debug, Clone, PartialEq)]
pub struct PercentMt {
    /// Per cell, in the counts' column order.
    pub values: Vec<f64>,
    pub definition: MtGenes,
    /// The features counted (rows, 0-based, in the counts' order) and their Seurat names.
    pub rows: Vec<usize>,
    pub features: Vec<String>,
    /// With a list: the distinct ids listed, and those no feature has (in the list's order).
    pub ids_listed: usize,
    pub ids_not_found: Vec<String>,
}

impl PercentMt {
    /// `PercentageFeatureSet()` of the features `mt` selects. A list that selects none is an error.
    pub fn compute(counts: &Counts, mt: &MtGenes) -> Result<PercentMt, String> {
        let feats = counts.features();
        let (rows, ids_listed, ids_not_found) = match mt {
            MtGenes::Pattern => (mt_pattern_rows(&feats), 0, Vec::new()),
            MtGenes::Ids { ids, file, .. } => {
                let distinct: HashSet<&str> = ids.iter().map(String::as_str).collect();
                let (rows, missing) = id_rows(&feats, ids);
                if rows.is_empty() {
                    return Err(format!(
                        "--mt-genes {file}: none of its {} id(s) is among the counts' features (column 1 of features.tsv){}",
                        distinct.len(),
                        if missing.is_empty() { String::new() } else { format!("; first: {}", missing.iter().take(NOT_FOUND_LISTED).cloned().collect::<Vec<_>>().join(", ")) }
                    ));
                }
                (rows, distinct.len(), missing)
            }
        };
        let names = counts.seurat_names();
        let features = rows.iter().map(|&r| names[r].clone()).collect();
        Ok(PercentMt { values: counts.percent_feature_set(&rows), definition: mt.clone(), rows, features, ids_listed, ids_not_found })
    }

    /// The warnings to print and record: ids not found, no feature matched, no counts at all.
    pub fn notes(&self) -> Vec<String> {
        let mut v = Vec::new();
        if let MtGenes::Ids { file, .. } = &self.definition {
            if !self.ids_not_found.is_empty() {
                v.push(format!(
                    "--mt-genes {file}: {} of its {} id(s) are not among the counts' features and are skipped (first: {})",
                    self.ids_not_found.len(),
                    self.ids_listed,
                    self.ids_not_found.iter().take(NOT_FOUND_LISTED).cloned().collect::<Vec<_>>().join(", ")
                ));
            }
        }
        if self.rows.is_empty() {
            v.push(format!("percent.mt: no feature's name matches \"{MT_PATTERN}\", so percent.mt is 0 in every cell; give the mitochondrial genes' ids with --mt-genes"));
        } else if self.values.iter().all(|&x| x == 0.0) {
            v.push(format!("percent.mt: the {} mitochondrial feature(s) have no counts, so percent.mt is 0 in every cell", self.rows.len()));
        }
        v
    }

    /// `summary.json`'s record (`result.sctransform.percent_mt`).
    pub fn record_json(&self) -> Json {
        let mut m: Vec<(&str, Json)> = Vec::new();
        match &self.definition {
            MtGenes::Pattern => {
                m.push(("definition", Json::Str("pattern".into())));
                m.push(("pattern", Json::Str(MT_PATTERN.into())));
                m.push(("matched_on", Json::Str("Seurat's row names: features.tsv's column 2, made unique (Read10X), with _ and | replaced by - (CreateSeuratObject)".into())));
            }
            MtGenes::Ids { file, sha256, .. } => {
                m.push(("definition", Json::Str("ids".into())));
                m.push(("file", Json::Str(file.clone())));
                m.push(("sha256", Json::Str(sha256.clone())));
                m.push(("matched_on", Json::Str("features.tsv's column 1".into())));
                m.push(("ids_listed", Json::Int(self.ids_listed as i64)));
                m.push(("ids_not_found", Json::Int(self.ids_not_found.len() as i64)));
                m.push(("ids_not_found_first", Json::strs(&self.ids_not_found.iter().take(NOT_FOUND_LISTED).cloned().collect::<Vec<_>>())));
            }
        }
        m.push(("n_features", Json::Int(self.rows.len() as i64)));
        m.push(("features_first", Json::strs(&self.features.iter().take(FEATURES_LISTED).cloned().collect::<Vec<_>>())));
        m.push(("cells_zero", Json::Int(self.values.iter().filter(|&&x| x == 0.0).count() as i64)));
        m.push(("cells_nan", Json::Int(self.values.iter().filter(|x| x.is_nan()).count() as i64)));
        Json::obj(m)
    }
}

/// What the chain needs besides the counts and the parameters to compute `percent.mt` and the
/// covariates, handed over in memory (docs/PLAN.md A3): the list `--mt-genes` gives (the pattern
/// when there is none), and `--metadata`'s columns for the cells, as text.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CellCovariates {
    pub mt_genes: Option<MtGenes>,
    /// The columns' names and, per column, each cell's text (`NA` where the table has no row
    /// for a cell), in the counts' column order ([`crate::metadata::extra_columns`]).
    pub metadata: Option<(Vec<String>, Vec<Vec<String>>)>,
}

/// R's reserved words (`?Reserved`), which `make.names()` does not leave as they are.
const RESERVED: [&str; 19] = [
    "if", "else", "repeat", "while", "function", "for", "in", "next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", "NA_character_", "NA_complex_",
];

/// Whether `make.names(name) == name` in R, for ASCII names: letters, digits, `.` and `_`, starting
/// with a letter or with a `.` not followed by a digit, and not a reserved word (`...` and `..1`,
/// `..2`, ... included). Names with other characters are treated as not syntactic.
pub fn is_syntactic(name: &str) -> bool {
    let b = name.as_bytes();
    if b.is_empty() || !b.iter().all(|c| c.is_ascii_alphanumeric() || *c == b'.' || *c == b'_') {
        return false;
    }
    let starts = b[0].is_ascii_alphabetic() || (b[0] == b'.' && !b.get(1).is_some_and(|c| c.is_ascii_digit()));
    let dots = name == "..." || (name.starts_with("..") && name.len() > 2 && name[2..].bytes().all(|c| c.is_ascii_digit()));
    starts && !dots && !RESERVED.contains(&name)
}

/// Checks `vars.to.regress`' names before any data is read: each is syntactic in R, none is
/// `GENE`, none repeats.
pub fn check_names(vars: &[String]) -> Result<(), String> {
    let mut seen = HashSet::new();
    for v in vars {
        if !seen.insert(v.as_str()) {
            return Err(format!("vars.to.regress: '{v}' is given twice"));
        }
        if v == "GENE" {
            return Err("vars.to.regress: 'GENE' is the response Seurat's RegressOutMatrix names in its formula (GENE ~ ...); rename the column".into());
        }
        if !is_syntactic(v) {
            return Err(format!("vars.to.regress: '{v}' is not a syntactic name in R (letters, digits, . and _, starting with a letter or a . not followed by a digit; not a reserved word): Seurat builds the formula GENE ~ {v} from it, which would not regress the column as named; rename the column"));
        }
    }
    Ok(())
}

/// The covariates' values, in `vars`' order (the formula's): `percent.mt`, `nCount_RNA` and
/// `nFeature_RNA` from `rna` (`percent.mt` must have been computed), any other name from
/// `metadata`'s column of that name, read as R reads numbers. `feature_names` are Seurat's row
/// names, which no covariate may share.
pub fn covariate_values(vars: &[String], rna: &RnaMetadata, metadata: Option<&(Vec<String>, Vec<Vec<String>>)>, feature_names: &[String]) -> Result<Vec<Vec<f64>>, String> {
    check_names(vars)?;
    let features: HashSet<&str> = feature_names.iter().map(String::as_str).collect();
    let n = rna.n_count.len();
    let mut out = Vec::with_capacity(vars.len());
    for v in vars {
        if features.contains(v.as_str()) {
            return Err(format!("vars.to.regress: '{v}' is also a feature's name; Seurat's ScaleData would add that feature's residuals as a covariate too; rename the column"));
        }
        let vals: Vec<f64> = match v.as_str() {
            PERCENT_MT => rna.percent_mt.clone().ok_or("vars.to.regress: percent.mt was not computed")?,
            "nCount_RNA" => rna.n_count.clone(),
            "nFeature_RNA" => rna.n_feature.iter().map(|&x| x as f64).collect(),
            name => {
                let Some((names, cols)) = metadata else {
                    return Err(format!("vars.to.regress: '{name}' is not percent.mt, nCount_RNA or nFeature_RNA, and no --metadata table is given"));
                };
                let Some(k) = names.iter().position(|c| c == name) else {
                    return Err(format!("vars.to.regress: '{name}' is not percent.mt, nCount_RNA, nFeature_RNA or a column of --metadata (its columns: {})", names.join(", ")));
                };
                let mut vals = Vec::with_capacity(n);
                for (i, t) in cols[k].iter().enumerate() {
                    let x = rcompat::as_numeric(t);
                    if rcompat::is_na(x) || x.is_nan() {
                        return Err(format!("vars.to.regress: --metadata column '{name}' has no number for cell {} ('{t}'): a missing or non-numeric value", i + 1));
                    }
                    vals.push(x);
                }
                vals
            }
        };
        if let Some(i) = vals.iter().position(|x| !x.is_finite()) {
            return Err(format!("vars.to.regress: '{v}' is not finite in cell {} ({}); R's lm() stops on NA/NaN/Inf", i + 1, scio::fmt17(vals[i])));
        }
        if vals.len() != n {
            return Err(format!("vars.to.regress: '{v}' has {} values for {n} cells", vals.len()));
        }
        out.push(vals);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn syntactic_names_as_make_names() {
        for ok in ["percent.mt", "nCount_RNA", "x", ".x", "a1", "a_b.c", "T", "F", "pi", ".a.2", "..a"] {
            assert!(is_syntactic(ok), "{ok}");
        }
        for bad in ["", "1x", ".2way", "my-col", "a b", "TRUE", "NULL", "if", "in", "...", "..1", "..12", "_x", "é", "a+b"] {
            assert!(!is_syntactic(bad), "{bad}");
        }
    }

    #[test]
    fn names_and_values() {
        let rna = RnaMetadata { n_count: vec![10.0, 20.0, 30.0], n_feature: vec![1, 2, 3], percent_mt: Some(vec![0.5, 1.5, 0.0]) };
        let md = (vec!["score".to_string(), "text".to_string(), "gappy".to_string()], vec![
            vec!["0.1".into(), "1e-3".into(), " 2 ".into()],
            vec!["a".into(), "b".into(), "c".into()],
            vec!["1".into(), "NA".into(), "3".into()],
        ]);
        let v = |s: &[&str]| s.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let got = covariate_values(&v(&["percent.mt", "nFeature_RNA", "score"]), &rna, Some(&md), &v(&["Gapdh"])).unwrap();
        assert_eq!(got[0], vec![0.5, 1.5, 0.0]);
        assert_eq!(got[1], vec![1.0, 2.0, 3.0]);
        assert_eq!(got[2].iter().map(|x| x.to_bits()).collect::<Vec<_>>(), [rcompat::as_numeric("0.1"), rcompat::as_numeric("1e-3"), 2.0].iter().map(|x| x.to_bits()).collect::<Vec<_>>());
        let err = |vars: &[&str], needle: &str| {
            let e = covariate_values(&v(vars), &rna, Some(&md), &v(&["Gapdh"])).unwrap_err();
            assert!(e.contains(needle), "{vars:?}: {e}");
        };
        err(&["text"], "non-numeric");
        err(&["gappy"], "missing");
        err(&["nope"], "not percent.mt");
        err(&["score", "score"], "twice");
        err(&["GENE"], "response");
        err(&["Gapdh"], "feature");
        let mut no_pm = rna.clone();
        no_pm.percent_mt = None;
        assert!(covariate_values(&v(&["percent.mt"]), &no_pm, None, &[]).is_err());
        assert!(covariate_values(&v(&["score"]), &rna, None, &[]).unwrap_err().contains("no --metadata"));
    }
}
