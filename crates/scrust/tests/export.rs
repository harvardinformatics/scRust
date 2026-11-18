//! Milestone M13 (docs/PLAN.md C3f, C5's K10): SCTransform's corrected counts and the export,
//! `scrust cluster --export` and `scrust sctransform --corrected-counts`, on the committed tier-T
//! objects:
//!
//! - **K10 end to end**: the corrected counts and their log1p against M8's portable instrumented
//!   runs (`reference/k10.tsv`: the SHA-256 of R's files' data lines, R's digests of the values,
//!   correct()'s median log_umi and its residuals beyond the clip range), `nCount_SCT` and
//!   `nFeature_SCT` against R's `metadata.tsv`. The products are OpenBLAS calls and the references
//!   portable runs, so the tests check the kernel first (Haswell, env.sh's default);
//! - **the export's layout**: `counts/` the input byte for byte, the manifest describing every
//!   file, the naming rules, and the 10x readers' own file lookup (Seurat's `Read10X`, scanpy's
//!   `read_10x_mtx` with and without `prefix=`) failing everywhere but `counts/`;
//! - **threads**: every file byte-identical at 1 and at 3 threads, the manifest's entries too.

use scio::manifest::Manifest;
use scio::{Json, Table};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const OBJECTS: [&str; 7] = [
    "filtered_seurat_tenx_test",
    "filtered_seurat_emptydrops_test",
    "soupx_tenx_test",
    "soupx_emptydrops_test",
    "decontx_raw_tenx_test",
    "decontx_raw_emptydrops_test",
    "cellbender_fromraw_test",
];

fn scrust(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_scrust")).args(args).output().expect("failed to run the scrust binary")
}

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
}

fn reference(object: &str) -> PathBuf {
    fixtures().join("clustering").join(object).join("reference")
}

/// The object's counts: the base objects' own, or M2's DoubletFinder input.
fn counts_dir(object: &str) -> PathBuf {
    let own = fixtures().join("clustering").join(object).join("input");
    if own.join("matrix.mtx.gz").exists() {
        own
    } else {
        fixtures().join("doubletfinder").join(object).join("input")
    }
}

fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("export").join(name);
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn k10(object: &str) -> BTreeMap<String, String> {
    let t = Table::read(&reference(object).join("k10.tsv")).unwrap();
    t.rows.iter().map(|r| (r[0].clone(), r[1].clone())).collect()
}

/// The SHA-256 of a Matrix Market file's data lines: the banner, the size line and the entries,
/// each ending in a newline, its comment lines left out (as validation/slurm/fixtures_m13.sh
/// digests R's).
fn data_lines_sha256(path: &Path) -> String {
    let text = scio::tenx::read_maybe_gz(path).unwrap();
    let mut h = scio::sha256::Sha256::new();
    for (i, line) in text.split(|&b| b == b'\n').enumerate() {
        if line.is_empty() || (i > 0 && line[0] == b'%') {
            continue;
        }
        h.update(line);
        h.update(b"\n");
    }
    h.hex()
}

fn gz_sha256(path: &Path) -> String {
    scio::sha256::hex(&scio::tenx::read_maybe_gz(path).unwrap())
}

fn column(t: &Table, name: &str) -> Vec<String> {
    let c = t.col(name).unwrap_or_else(|| panic!("no column {name}"));
    t.rows.iter().map(|r| r[c].clone()).collect()
}

fn gz_table(path: &Path) -> Table {
    Table::parse(&String::from_utf8(scio::tenx::read_maybe_gz(path).unwrap()).unwrap()).unwrap()
}

/// `scrust cluster --counts <object> --preset workflow --export` into `out`.
fn export(object: &str, out: &Path, threads: &str, extra: &[&str]) {
    let c = counts_dir(object).display().to_string();
    let o = out.display().to_string();
    let mut a = vec!["cluster", "--counts", c.as_str(), "--preset", "workflow", "--export", "--out", o.as_str(), "--threads", threads];
    a.extend_from_slice(extra);
    let r = scrust(&a);
    assert!(r.status.success(), "{object}: {}", String::from_utf8_lossy(&r.stderr));
}

/// What the corrected counts of `dir` (an export's `sct/`, or `scrust sctransform
/// --corrected-counts`' output) must hold for `object`.
fn check_k10(object: &str, dir: &Path, corrected: &Json) {
    let want = k10(object);
    let w = |k: &str| want.get(k).unwrap_or_else(|| panic!("k10.tsv: no {k}")).as_str();
    assert_eq!(data_lines_sha256(&dir.join("corrected_counts.mtx.gz")), w("corrected_counts_data_lines_sha256"), "{object}: corrected counts");
    assert_eq!(data_lines_sha256(&dir.join("log1p_corrected_counts.mtx.gz")), w("log1p_data_lines_sha256"), "{object}: log1p data");
    assert_eq!(gz_sha256(&dir.join("sct_features.tsv.gz")), w("sct_features_sha256"), "{object}: sct_features");
    assert_eq!(gz_sha256(&dir.join("sct_cells.tsv.gz")), w("sct_cells_sha256"), "{object}: sct_cells");
    let s = |k: &str| corrected.get(k).and_then(Json::as_str).unwrap_or_default().to_string();
    let n = |k: &str| corrected.get(k).and_then(Json::as_f64).unwrap_or(f64::NAN);
    assert_eq!(s("counts_x_sha256"), w("sct_counts_x_sha256"), "{object}: R's digest of the counts");
    assert_eq!(s("data_x_sha256"), w("sct_data_x_sha256"), "{object}: R's digest of the log1p data");
    assert_eq!(n("nnz") as usize, w("sct_nnz").parse::<usize>().unwrap());
    assert_eq!(n("n_features") as usize, w("sct_n_features").parse::<usize>().unwrap());
    assert_eq!(n("median_log_umi").to_bits(), w("correct_median_log_umi").parse::<f64>().unwrap().to_bits(), "{object}: the median log_umi");
    assert_eq!(n("residuals_beyond_clip") as u64, w("correct_residuals_beyond_clip").parse::<u64>().unwrap(), "{object}: the residuals beyond the clip range");
    // the comment lines say what the matrix is, and that it is not raw counts
    let text = String::from_utf8(scio::tenx::read_maybe_gz(&dir.join("corrected_counts.mtx.gz")).unwrap()).unwrap();
    let head: Vec<&str> = text.lines().take(12).collect();
    assert!(head[0] == "%%MatrixMarket matrix coordinate integer general" && head[1].starts_with("% scrust: SCTransform's CORRECTED counts") && head.iter().any(|l| l.contains("NOT raw counts")), "{head:?}");
    assert!(head.iter().any(|l| l.starts_with("% Made by scrust ")) && head.iter().filter(|l| l.starts_with("% Input sha256: ")).count() == 3, "{head:?}");
}

/// The files 10x readers would read in `dir`, if any: Seurat's `Read10X(dir)` needs
/// `barcodes.tsv(.gz)`, `features.tsv(.gz)` or `genes.tsv`, and `matrix.mtx(.gz)`; scanpy's
/// `read_10x_mtx(dir, prefix=p)` needs `<p>matrix.mtx.gz` (or `<p>matrix.mtx` with `<p>genes.tsv`).
fn tenx_readable(dir: &Path, prefix: &str) -> bool {
    let f = |n: &str| dir.join(format!("{prefix}{n}")).is_file();
    let matrix = f("matrix.mtx") || f("matrix.mtx.gz");
    let barcodes = f("barcodes.tsv") || f("barcodes.tsv.gz");
    let features = f("features.tsv") || f("features.tsv.gz") || f("genes.tsv");
    matrix && barcodes && features
}

/// The naming rules on a real export: only `counts/` is a 10x directory, and no prefix makes a
/// normalised file readable as one (scanpy's `prefix=`: every prefix of every file name).
fn check_naming(out: &Path) {
    let files = scio::manifest::list_files(out).unwrap();
    let mut dirs: Vec<PathBuf> = vec![out.to_path_buf()];
    dirs.extend(files.iter().filter_map(|f| Path::new(f).parent().filter(|p| !p.as_os_str().is_empty()).map(|p| out.join(p))));
    dirs.sort();
    dirs.dedup();
    for d in &dirs {
        let readable = tenx_readable(d, "");
        assert_eq!(readable, d.ends_with("counts"), "{}: a 10x reader {} it", d.display(), if readable { "reads" } else { "does not read" });
        // scanpy with prefix=: every prefix of every file name in this directory (and, at the top,
        // of the files under it, as "sct/..." prefixes)
        for f in &files {
            let name = f.as_str();
            for k in 0..=name.len() {
                let p = &name[..k];
                if !name.is_char_boundary(k) {
                    continue;
                }
                let target = d.join(format!("{p}matrix.mtx.gz"));
                if target.is_file() {
                    let rel = target.strip_prefix(out).unwrap().display().to_string();
                    assert!(rel.starts_with("counts/"), "prefix {p:?} on {} reads {rel}", d.display());
                }
            }
        }
    }
    for f in &files {
        if !f.starts_with("counts/") {
            assert!(!f.ends_with("matrix.mtx") && !f.ends_with("matrix.mtx.gz"), "{f}");
            let name = Path::new(f).file_name().unwrap().to_str().unwrap();
            assert!(!scio::manifest::TENX_NAMES.contains(&name), "{f}: a 10x name outside counts/");
        }
    }
    let m = Manifest::read(out).unwrap();
    assert!(scio::manifest::naming_problems(&m.entries).is_empty());
}

/// K10 end to end on the seven committed objects, the export's layout, and the manifest.
#[test]
fn the_export_reproduces_r_s_corrected_counts() {
    for object in OBJECTS {
        rcompat::blas::require_fixture_kernel(&reference(object), "scrust cluster --export (K10)");
        let out = tmp(object);
        export(object, &out, "2", &[]);
        let summary = Json::parse(&std::fs::read_to_string(out.join("summary.json")).unwrap()).unwrap();
        let corrected = summary.get("result").and_then(|r| r.get("sctransform")).and_then(|s| s.get("corrected")).cloned().expect("result.sctransform.corrected");
        check_k10(object, &out.join("sct"), &corrected);
        // nCount_SCT and nFeature_SCT against R's metadata.tsv, after nFeature_RNA
        let r_meta = gz_table(&reference(object).join("metadata.tsv.gz"));
        let c_meta = Table::read(&out.join("metadata.tsv")).unwrap();
        for col in ["nCount_SCT", "nFeature_SCT"] {
            assert_eq!(column(&c_meta, col), column(&r_meta, col), "{object}: {col}");
        }
        assert_eq!(&c_meta.header[..5], &["barcode", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT"]);
        // counts/: the input byte for byte; the manifest describes every file
        for f in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"] {
            assert!(std::fs::read(counts_dir(object).join(f)).unwrap() == std::fs::read(out.join("counts").join(f)).unwrap(), "{object}: counts/{f}");
        }
        let m = Manifest::read(&out).unwrap();
        assert!(m.check(&out, true).unwrap().is_empty(), "{object}: {:?}", m.check(&out, true).unwrap());
        assert_eq!(m.entry("sct/corrected_counts.mtx.gz").unwrap().scale, "corrected_counts");
        assert_eq!(m.entry("counts/matrix.mtx.gz").unwrap().scale, "raw_counts");
        let input = &m.provenance.inputs[0];
        assert_eq!(input.files.iter().map(|f| f.sha256.clone()).collect::<Vec<_>>(), ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"].iter().map(|f| m.entry(&format!("counts/{f}")).unwrap().sha256.clone().unwrap()).collect::<Vec<_>>());
        assert!(std::fs::read_to_string(out.join("README.txt")).unwrap().contains("THE ONLY RAW COUNTS ARE IN counts/"));
        check_naming(&out);
        let v = scrust(&["validate", "manifest", "--dir", &out.display().to_string()]);
        assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    }
}

/// Every file of an export at 1 and at 3 threads, and with the stage dumps, byte for byte; the
/// manifest's entries too (its provenance's command line names the thread count).
#[test]
fn the_export_does_not_depend_on_the_thread_count() {
    let object = "filtered_seurat_emptydrops_test";
    let (a, b, c) = (tmp("threads_1"), tmp("threads_3"), tmp("threads_dumps"));
    export(object, &a, "1", &[]);
    export(object, &b, "3", &[]);
    export(object, &c, "2", &["--dump-stages", "full", "--keep-stages"]);
    let files = scio::manifest::list_files(&a).unwrap();
    assert_eq!(files, scio::manifest::list_files(&b).unwrap());
    for f in &files {
        if f == "summary.json" {
            continue;
        }
        assert!(std::fs::read(a.join(f)).unwrap() == std::fs::read(b.join(f)).unwrap(), "{f} depends on the thread count");
        // and on the stage dumps: the normalised matrices' comment lines leave them out (the
        // README lists the directory's files, the stage dumps among them)
        if f != "README.txt" {
            assert!(std::fs::read(a.join(f)).unwrap() == std::fs::read(c.join(f)).unwrap(), "{f} depends on --dump-stages");
        }
    }
    let (ma, mb) = (Manifest::read(&a).unwrap(), Manifest::read(&b).unwrap());
    assert_eq!(ma.entries, mb.entries);
    assert_eq!(ma.provenance.command, mb.provenance.command);
    assert_ne!(ma.provenance.command_line, mb.provenance.command_line);
    let result = |d: &Path| scio::json::to_file_text(Json::parse(&std::fs::read_to_string(d.join("summary.json")).unwrap()).unwrap().get("result").unwrap());
    assert_eq!(result(&a), result(&b));
    // the stage dumps hold K10 in R's layout
    let st = scio::stages::StageDir::open(&c.join("stages")).unwrap();
    assert!(st.has("K10_genes.npy", "chain"));
    let rng = Table::read(&c.join("stages/rng.tsv")).unwrap();
    let labels = column(&rng, "label");
    let at = |l: &str| labels.iter().position(|x| x == l).unwrap();
    assert!(at("fits_exit") < at("correct_entry") && at("correct_entry") < at("vst_exit"), "R's order: {labels:?}");
}

/// The routes that cannot export refuse `--export`, with the reason.
#[test]
fn only_a_run_from_the_counts_exports() {
    let pca = reference("filtered_seurat_tenx_test").join("pca");
    let o = tmp("refused");
    let r = scrust(&["cluster", "--pca", &pca.display().to_string(), "--export", "--out", &o.display().to_string()]);
    assert!(!r.status.success());
    let e = String::from_utf8_lossy(&r.stderr);
    assert!(e.contains("--export needs --counts") && e.contains("vst's residuals"), "{e}");
    assert!(!o.join("summary.json").exists());
}

/// `scrust sctransform --corrected-counts` (seed 12345, the workflow's): the SCT part of the
/// export, R's corrected counts, its README and manifest; without the option, no such file.
#[test]
fn sctransform_writes_the_sct_part() {
    let object = "soupx_tenx_test";
    rcompat::blas::require_fixture_kernel(&reference(object), "scrust sctransform --corrected-counts (K10)");
    let c = counts_dir(object).display().to_string();
    let out = tmp("sctransform_corrected");
    let r = scrust(&["sctransform", "--counts", &c, "--seed", "12345", "--corrected-counts", "--out", &out.display().to_string(), "--threads", "2"]);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    let summary = Json::parse(&std::fs::read_to_string(out.join("summary.json")).unwrap()).unwrap();
    check_k10(object, &out, summary.get("result").unwrap().get("corrected").unwrap());
    let m = Manifest::read(&out).unwrap();
    assert!(m.check(&out, true).unwrap().is_empty(), "{:?}", m.check(&out, true).unwrap());
    assert!(!tenx_readable(&out, ""));
    assert!(out.join("README.txt").is_file());
    let plain = tmp("sctransform_plain");
    let r = scrust(&["sctransform", "--counts", &c, "--seed", "12345", "--out", &plain.display().to_string()]);
    assert!(r.status.success());
    assert!(!plain.join("corrected_counts.mtx.gz").exists() && !plain.join("README.txt").exists());
    // the files both write are the same bytes
    for f in ["cell_attr.tsv", "gene_attr.tsv", "model_pars.tsv", "model_pars_fit.tsv", "variable_features.tsv", "scale_data.npy", "cell_metadata.tsv"] {
        assert!(std::fs::read(out.join(f)).unwrap() == std::fs::read(plain.join(f)).unwrap(), "{f}");
    }
    let mp = Manifest::read(&plain).unwrap();
    assert!(mp.check(&plain, true).unwrap().is_empty());
    let r = scrust(&["sctransform", "--counts", &c, "--corrected-counts", "--from-stages", "x", "--context", "final", "--out", &tmp("refused_sct").display().to_string()]);
    assert!(!r.status.success() && String::from_utf8_lossy(&r.stderr).contains("--corrected-counts: from the counts only"));
}
