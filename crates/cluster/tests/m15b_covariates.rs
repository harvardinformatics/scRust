//! Milestone M15b: `percent.mt` against Seurat 5.5.1 (validation/R/regress_vectors.R; fixtures in
//! tests/fixtures/clustering/percent_mt/): a 10x directory with tricky feature names (case,
//! underscores and pipes Seurat turns into dashes, names Read10X makes unique, a mitochondrial
//! feature without counts, a cell without mitochondrial counts), with integer counts and with
//! non-integer ones, and R's results:
//!
//! - Seurat's row names (`Read10X()`'s `make.unique()`, then `CreateSeuratObject()`'s `_` and `|`
//!   to `-`), and `make.unique()` on its own (`make_unique.tsv`);
//! - `PercentageFeatureSet(obj, pattern = "(?i)^mt-")`: the features matched and the values, bit
//!   for bit;
//! - `PercentageFeatureSet(obj, features = <the features whose id the list names>)`, the list
//!   read as `--mt-genes` reads it (blank lines, comments, an id twice, ids not found);
//! - `nCount_RNA`, `PercentageFeatureSet()`'s denominator.
//!
//! No BLAS is involved: the checks are exact on any machine.

use cluster::covariates::{MtGenes, PercentMt};
use scio::features::make_unique;
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/percent_mt")
}

fn hexl(s: &str) -> Vec<f64> {
    s.split(',').map(|h| f64::from_bits(u64::from_str_radix(h, 16).unwrap())).collect()
}

fn list(s: &str) -> Vec<String> {
    if s == "-" {
        Vec::new()
    } else {
        s.split(',').map(str::to_owned).collect()
    }
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn make_unique_as_r() {
    let text = std::fs::read_to_string(dir().join("make_unique.tsv")).unwrap();
    let mut n = 0;
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(make_unique(&list(f[0])), list(f[1]), "make.unique({})", f[0]);
        n += 1;
    }
    assert_eq!(n, 5);
}

#[test]
fn percent_mt_as_seurat_computes_it() {
    let expected = std::fs::read_to_string(dir().join("expected.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = expected.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 6);
    let want_names: Vec<String> = std::fs::read_to_string(dir().join("seurat_names.txt")).unwrap().lines().map(str::to_owned).collect();
    let ids_file = dir().join("mt_genes.txt");
    let ids_text = std::fs::read_to_string(&ids_file).unwrap();
    for kind in ["integer", "real"] {
        let c = scio::Counts::read(&dir().join(kind)).unwrap();
        assert_eq!(c.seurat_names(), want_names, "{kind}: Seurat's row names");
        let rna = cluster::metadata::RnaMetadata::of(&c.matrix);
        for r in rows.iter().filter(|r| r[0] == kind) {
            match r[1] {
                "nCount_RNA" => assert_eq!(bits(&rna.n_count), bits(&hexl(r[3])), "{kind}: nCount_RNA"),
                "pattern" => {
                    let p = PercentMt::compute(&c, &MtGenes::Pattern).unwrap();
                    assert_eq!(p.features, list(r[2]), "{kind}: the features the pattern matches");
                    assert_eq!(bits(&p.values), bits(&hexl(r[3])), "{kind}: percent.mt by the pattern");
                    assert!(p.notes().is_empty(), "{:?}", p.notes());
                }
                "list" => {
                    let mt = MtGenes::from_file_text(&ids_text, &ids_file.display().to_string(), "-");
                    let p = PercentMt::compute(&c, &mt).unwrap();
                    assert_eq!(p.features, list(r[2]), "{kind}: the listed features");
                    assert_eq!(bits(&p.values), bits(&hexl(r[3])), "{kind}: percent.mt by the list");
                    assert_eq!((p.ids_listed, p.ids_not_found.clone()), (6, vec!["NOPE1".to_string(), "NOPE2".to_string()]));
                    assert_eq!(p.notes().len(), 1, "the ids not found are reported");
                }
                other => panic!("unknown row {other}"),
            }
        }
    }
    // a list that matches nothing is an error; a pattern that matches nothing gives zeros
    let c = scio::Counts::read(&dir().join("integer")).unwrap();
    let none = MtGenes::Ids { ids: vec!["NOPE".into()], file: "x".into(), sha256: "-".into() };
    assert!(PercentMt::compute(&c, &none).unwrap_err().contains("none of its 1 id(s)"));
}
