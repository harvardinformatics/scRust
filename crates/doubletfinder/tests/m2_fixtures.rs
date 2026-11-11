//! The tier-T fixtures of milestone M2 (tests/fixtures/doubletfinder/, docs/PLAN.md B4 and B6).
//!
//! No DoubletFinder code yet (that starts at M3): these tests check that the committed fixtures
//! are complete and consistent, and that their text files carry exact doubles. The pANN values
//! in cells.tsv are written with 17 significant digits, so Rust's correctly rounded parser must
//! give back the very bits of the binary copy in stages/S5_pann_final.npy; and the direct
//! reference run must reproduce the saved workflow output (doublet set, singlet pANN, columns).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

const OBJECTS: [&str; 5] = [
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn object_dir(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/doubletfinder")
        .join(object)
}

/// A .npy array as written by validation/R/lib/io.R: format 1.0, little-endian `<f8` or `<i4`.
struct Npy {
    dtype: String,
    shape: Vec<usize>,
    data: Vec<u8>,
}

impl Npy {
    fn read(path: &Path) -> Npy {
        let bytes = fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert_eq!(&bytes[..6], b"\x93NUMPY", "{}: not a .npy file", path.display());
        let hlen = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let header = std::str::from_utf8(&bytes[10..10 + hlen]).unwrap();
        let field = |key: &str| {
            let start = header.find(key).unwrap() + key.len();
            header[start..].trim_start_matches([':', ' ', '\''])
        };
        let dtype = field("'descr'").split('\'').next().unwrap().to_string();
        let shape_txt = field("'shape'");
        let shape = shape_txt[1..shape_txt.find(')').unwrap()]
            .split(',')
            .filter_map(|s| s.trim().parse().ok())
            .collect();
        Npy { dtype, shape, data: bytes[10 + hlen..].to_vec() }
    }

    fn rows(&self) -> usize {
        self.shape.first().copied().unwrap_or(1)
    }

    fn cols(&self) -> usize {
        if self.shape.len() > 1 { self.shape[1] } else { 1 }
    }

    fn f64s(&self) -> Vec<f64> {
        assert_eq!(self.dtype, "<f8");
        self.data.as_chunks::<8>().0.iter().map(|c| f64::from_le_bytes(*c)).collect()
    }

    fn i32s(&self) -> Vec<i32> {
        assert_eq!(self.dtype, "<i4");
        self.data.as_chunks::<4>().0.iter().map(|c| i32::from_le_bytes(*c)).collect()
    }
}

/// Rows of a TSV file with a header line, as (header, rows).
fn tsv(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').map(str::to_owned).collect();
    let rows = lines.map(|l| l.split('\t').map(str::to_owned).collect()).collect();
    (header, rows)
}

/// The string value of `"key": "value"` in one of the fixtures' JSON files.
fn json_string(path: &Path, key: &str) -> String {
    let text = fs::read_to_string(path).unwrap();
    let pat = format!("\"{key}\": \"");
    let start = text.find(&pat).unwrap_or_else(|| panic!("{key} not in {}", path.display())) + pat.len();
    text[start..start + text[start..].find('"').unwrap()].to_string()
}

#[test]
fn final_pann_text_is_exact_and_matches_the_calls() {
    for object in OBJECTS {
        let dir = object_dir(object);
        let (header, cells) = tsv(&dir.join("reference/cells.tsv"));
        assert_eq!(header, ["barcode", "pANN", "classification"]);
        let pann = Npy::read(&dir.join("stages/S5_pann_final.npy")).f64s();
        let calls = Npy::read(&dir.join("stages/S7_doublet.npy")).i32s();
        assert_eq!(cells.len(), pann.len(), "{object}");
        for ((row, &p), &c) in cells.iter().zip(&pann).zip(&calls) {
            let parsed: f64 = row[1].parse().unwrap();
            assert_eq!(parsed.to_bits(), p.to_bits(), "{object} {}: {} vs {p:e}", row[0], row[1]);
            assert_eq!(row[2] == "Doublet", c == 1, "{object} {}", row[0]);
        }
    }
}

#[test]
fn reference_reproduces_the_saved_workflow_output() {
    for object in OBJECTS {
        let dir = object_dir(object);
        let (_, cells) = tsv(&dir.join("reference/cells.tsv"));
        let doublets: BTreeSet<&str> = cells
            .iter()
            .filter(|r| r[2] == "Doublet")
            .map(|r| r[0].as_str())
            .collect();
        let saved_doublets_text = fs::read_to_string(dir.join("saved/doublets.txt")).unwrap();
        let saved_doublets: BTreeSet<&str> = saved_doublets_text.lines().collect();
        assert_eq!(doublets, saved_doublets, "{object}: doublet sets differ");

        let (_, saved) = tsv(&dir.join("saved/cells.tsv"));
        let singlets: Vec<&Vec<String>> = cells.iter().filter(|r| r[2] == "Singlet").collect();
        assert_eq!(saved.len(), singlets.len(), "{object}");
        for (s, r) in saved.iter().zip(&singlets) {
            assert_eq!(s[0], r[0], "{object}: kept cells differ");
            assert_eq!(s[1], r[1], "{object} {}: singlet pANN differs", s[0]);
        }
        for key in ["pANN_column", "classification_column", "pK_label"] {
            assert_eq!(
                json_string(&dir.join("reference/summary.json"), key),
                json_string(&dir.join("saved/saved.json"), key),
                "{object}: {key}"
            );
        }
    }
}

#[test]
fn stage_index_matches_the_arrays() {
    // Files that stay in data/ (the size budget, tests/fixtures/README.md).
    let data_only = |f: &str| {
        f == "S6_bkde.npy"
            || ["S2d_y.", "S2d_scale_data.", "S4_positions."]
                .iter()
                .any(|p| f.starts_with(p))
    };
    for object in OBJECTS {
        let stages = object_dir(object).join("stages");
        let (header, rows) = tsv(&stages.join("index.tsv"));
        assert_eq!(header, ["file", "ctx", "row_start", "nrow", "ncol", "dtype", "stage"]);
        let mut contexts = BTreeSet::new();
        for row in &rows {
            let (file, ctx) = (&row[0], &row[1]);
            contexts.insert(ctx.split(':').next().unwrap().to_string());
            let path = stages.join(file);
            if !path.exists() {
                assert!(data_only(file), "{object}: {file} is missing");
                continue;
            }
            let a = Npy::read(&path);
            let (start, nrow, ncol): (usize, usize, usize) =
                (row[2].parse().unwrap(), row[3].parse().unwrap(), row[4].parse().unwrap());
            assert_eq!(a.dtype, row[5], "{object} {file}");
            assert!(start + nrow <= a.rows(), "{object} {file} {ctx}: rows out of range");
            assert_eq!(a.cols(), ncol, "{object} {file} {ctx}");
            assert!(row[6].starts_with('S'), "{object} {file}: stage {}", row[6]);
        }
        let want: BTreeSet<String> = ["final", "summary", "sweep_pN0.05", "sweep_pN0.1", "sweep_pN0.15",
            "sweep_pN0.2", "sweep_pN0.25", "sweep_pN0.3"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(contexts, want, "{object}: contexts");
    }
}

#[test]
fn rng_checkpoints_are_mersenne_twister_states() {
    for object in OBJECTS {
        let stages = object_dir(object).join("stages");
        let rng = Npy::read(&stages.join("rng.npy"));
        assert_eq!(rng.cols(), 626, "{object}");
        let (_, labels) = tsv(&stages.join("rng.tsv"));
        assert_eq!(labels.len(), rng.rows(), "{object}");
        // Fortran order: the first column holds .Random.seed[1] of every checkpoint, 10403 for
        // Mersenne-Twister with Inversion and Rejection sampling (R >= 3.6 defaults).
        let seeds = rng.i32s();
        assert!(seeds[..rng.rows()].iter().all(|&kind| kind == 10403), "{object}");
    }
}
