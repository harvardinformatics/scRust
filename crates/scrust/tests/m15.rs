//! Milestone M15 (docs/PLAN.md C3h; A3's in-memory hand-offs): clustering on the fly and the
//! singlets' re-clustering in `scrust doubletfinder`, and the in-memory rule's test.
//!
//! - **The in-memory rule:** the chain's entry point called in this process, its result written
//!   by the library's writer, gives the files `scrust cluster --counts` writes, byte for byte, and
//!   the same `result` section.
//! - **Clustering on the fly:** `scrust doubletfinder` given no `--clusters` writes DoubletFinder's
//!   files exactly as a run given the labels it computed; its own labels and levels are R's (under
//!   the fixtures' kernel); its manifest, nested outputs included, describes the directory.
//! - **Re-clustering the singlets in process** (`--recluster`, with and without the export) gives
//!   the two-command route's files (`scrust doubletfinder`, then `scrust cluster --counts
//!   <out>/singlets`), byte for byte, apart from the timings and the command line in its
//!   `summary.json`. The routes run in separate working directories with the same relative paths,
//!   so the paths they record are the same.
//!
//! All but the comparison with R's labels are Rust against Rust: any OpenBLAS kernel.

use cluster::chain::{self, Options};
use cluster::chain_output::Writer;
use cluster::ChainParams;
use scio::Json;
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn scrust_in(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_scrust")).current_dir(dir).args(args).output().expect("failed to run the scrust binary")
}

fn ok(dir: &Path, args: &[&str]) -> Output {
    let o = scrust_in(dir, args);
    assert!(o.status.success(), "{args:?}: {}", String::from_utf8_lossy(&o.stderr));
    o
}

fn fixture(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder").join(object)
}

fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m15").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Every file under `dir`, by relative path.
fn tree(dir: &Path) -> BTreeMap<String, PathBuf> {
    let mut out = BTreeMap::new();
    let mut todo = vec![dir.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                todo.push(p);
            } else {
                out.insert(p.strip_prefix(dir).unwrap().display().to_string(), p);
            }
        }
    }
    out
}

/// The members of `summary.json` that record how a run went, not what it computed.
const VOLATILE: [&str; 17] = [
    "start_time", "elapsed_s", "phases", "peak_rss_kib", "printed_elapsed", "sctransform_s", "pca_s", "neighbors_s", "cluster_s", "silhouette_s", "purity_s", "network_s", "resolutions_s",
    "optimizer_s", "groupsingletons_s", "export_counts_copy_s", "manifest_s",
];

/// `j` without the volatile members, and without `files.provenance.command_line`.
fn settled(j: &Json, path: &str) -> Json {
    match j {
        Json::Object(m) => Json::Object(
            m.iter()
                .filter(|(k, _)| !VOLATILE.contains(&k.as_str()) && format!("{path}.{k}") != ".files.provenance.command_line")
                .map(|(k, v)| (k.clone(), settled(v, &format!("{path}.{k}"))))
                .collect(),
        ),
        Json::List(v) => Json::List(v.iter().map(|x| settled(x, path)).collect()),
        Json::Vector(v) => Json::Vector(v.iter().map(|x| settled(x, path)).collect()),
        other => other.clone(),
    }
}

fn read_json(p: &Path) -> Json {
    Json::parse(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// The same files under `a` and `b`, byte for byte, but the `summary.json` files, which must
/// agree as JSON once the timings and the command line are set aside. `skip_a` and `skip_b`: paths
/// (or directories, ending in `/`) on one side only, left out.
fn assert_same_tree(a: &Path, b: &Path, skip_a: &[&str], skip_b: &[&str], skip_summary: &[&str]) {
    let keep = |rel: &String, skip: &[&str]| !skip.iter().any(|s| rel == s || (s.ends_with('/') && rel.starts_with(s)));
    let ta: BTreeMap<String, PathBuf> = tree(a).into_iter().filter(|(k, _)| keep(k, skip_a)).collect();
    let tb: BTreeMap<String, PathBuf> = tree(b).into_iter().filter(|(k, _)| keep(k, skip_b)).collect();
    assert_eq!(ta.keys().collect::<Vec<_>>(), tb.keys().collect::<Vec<_>>(), "{} and {} hold different files", a.display(), b.display());
    for (rel, pa) in &ta {
        let pb = &tb[rel];
        if rel.ends_with("summary.json") {
            if !skip_summary.contains(&rel.as_str()) {
                assert_eq!(settled(&read_json(pa), ""), settled(&read_json(pb), ""), "{rel} differs beyond the timings and the command line");
            }
        } else {
            assert!(std::fs::read(pa).unwrap() == std::fs::read(pb).unwrap(), "{rel} differs");
        }
    }
}

/// The in-memory rule (A3): the entry point's result, written by the library's writer, is
/// `scrust cluster --counts`' output, file for file, with the export, and the same `result`.
#[test]
fn the_entry_point_s_result_writes_the_command_s_files() {
    let input = fixture("decontx_raw_tenx_test").join("input");
    let d = tmp("in_memory_rule");
    let (cmd_dir, lib_dir) = (d.join("command"), d.join("library"));
    let c = input.display().to_string();
    ok(&d, &["cluster", "--counts", &c, "--preset", "workflow", "--export", "--out", cmd_dir.to_str().unwrap(), "--threads", "2"]);
    // in this process: read, call, write
    let counts = scio::Counts::read(&input).unwrap();
    let opts = Options { correct: true, keep_reductions: true, keep_sct: true, ..Options::default() };
    let r = chain::cluster_counts(Cow::Owned(counts), &ChainParams::workflow(), 3, opts, &mut chain::NoObserver).unwrap();
    let summary = read_json(&cmd_dir.join("summary.json"));
    let prov = summary.get("files").and_then(|f| f.get("provenance")).unwrap();
    std::fs::create_dir_all(&lib_dir).unwrap();
    let mut w = Writer::new(&lib_dir, r.cells.clone(), 3);
    w.write_stages = false;
    w.export = true;
    // the comment lines of the normalised matrices name the command and the input's checksums
    w.command = prov.get("command").and_then(|v| v.as_str()).unwrap().to_string();
    let Some(Json::List(inputs)) = prov.get("inputs") else { panic!("no provenance.inputs") };
    let Some(Json::List(files)) = inputs[0].get("files") else { panic!("no input files") };
    w.input_checksums = files.iter().map(|f| (f.get("sha256").and_then(|v| v.as_str()).unwrap().to_string(), input.join(f.get("name").and_then(|v| v.as_str()).unwrap()).display().to_string())).collect();
    w.write_clustered(&r).unwrap();
    w.finish_labels(&r.chain.clusters, r.chain.silhouette.as_ref(), r.chain.purity.as_ref()).unwrap();
    // every file the writer wrote is the command's, byte for byte; the command adds counts/ (the
    // input, copied), README.txt and summary.json around them
    let written = tree(&lib_dir);
    assert!(written.contains_key("sct/corrected_counts.mtx.gz") && written.contains_key("pca/embeddings.tsv") && written.contains_key("neighbors/snn.mtx.gz") && written.contains_key("metadata.tsv"));
    let all: Vec<String> = tree(&cmd_dir).into_keys().collect();
    let extra: Vec<&String> = all.iter().filter(|k| !written.contains_key(*k)).collect();
    assert_eq!(extra, ["README.txt", "counts/barcodes.tsv.gz", "counts/features.tsv.gz", "counts/matrix.mtx.gz", "summary.json"], "the command's files beyond the writer's");
    for (rel, p) in &written {
        assert!(std::fs::read(p).unwrap() == std::fs::read(cmd_dir.join(rel)).unwrap(), "{rel} differs");
    }
    // the result section, from the result alone, as the file's text holds it
    let as_text = |j: &Json| Json::parse(&scio::json::to_file_text(j)).unwrap();
    let mine = as_text(&cluster::chain_output::clustered_result_json(&r));
    let theirs = summary.get("result").unwrap();
    if let (Json::Object(a), Json::Object(b)) = (&mine, theirs) {
        let ka: Vec<&String> = a.iter().map(|(k, _)| k).collect();
        let kb: Vec<&String> = b.iter().map(|(k, _)| k).collect();
        assert_eq!(ka, kb, "result's members");
        for ((k, x), (_, y)) in a.iter().zip(b) {
            assert_eq!(x, y, "result.{k}");
        }
    }
    assert_eq!(&mine, theirs, "summary.json's result");
    assert_eq!(&as_text(&w.result_json(&r.chain.clusters)), theirs);
}

/// Clustering on the fly: DoubletFinder's files are those of a run given the labels the hook
/// computed; `clusters/` holds them in `--clusters`' layout; the manifest, with its nested
/// output, describes the directory.
#[test]
fn doubletfinder_clusters_its_own_input() {
    let f = fixture("soupx_emptydrops_test");
    let input = f.join("input").display().to_string();
    let d = tmp("on_the_fly");
    let o = ok(&d, &["doubletfinder", "--counts", &input, "--out", "hook", "--cluster-stats", "silhouette,purity", "--threads", "2"]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("clusters/: ") && err.contains("(on the fly)"), "{err}");
    let hook = d.join("hook");
    for file in ["seurat_clusters.tsv", "seurat_clusters_levels.txt", "metadata.tsv", "clusters.tsv", "cluster_ids.txt", "nclusters.txt", "summary.json"] {
        assert!(hook.join("clusters").join(file).is_file(), "clusters/{file}");
    }
    assert!(!hook.join("clusters/pca").exists() && !hook.join("clusters/neighbors").exists(), "only the labels are kept");
    let meta = std::fs::read_to_string(hook.join("clusters/metadata.tsv")).unwrap();
    assert!(meta.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tSCT_snn_res.0.8\tseurat_clusters\tsilhouette_width\tsilhouette_other\tneighborhood_purity\tpurity_maximum\n"), "{}", &meta[..120]);
    // the same run given those labels: DoubletFinder's files and result, byte for byte
    let labels = hook.join("clusters/seurat_clusters.tsv").display().to_string();
    ok(&d, &["doubletfinder", "--counts", &input, "--clusters", &labels, "--out", "given", "--threads", "3"]);
    let given = d.join("given");
    for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv", "singlets/matrix.mtx.gz", "singlets/barcodes.tsv.gz", "singlets/features.tsv.gz"] {
        assert!(std::fs::read(hook.join(file)).unwrap() == std::fs::read(given.join(file)).unwrap(), "{file}");
    }
    let (rh, rg) = (read_json(&hook.join("summary.json")), read_json(&given.join("summary.json")));
    let Some(Json::Object(h)) = rh.get("result") else { panic!("no result") };
    let Some(Json::Object(g)) = rg.get("result") else { panic!("no result") };
    let h_df: Vec<&(String, Json)> = h.iter().filter(|(k, _)| k != "clustering").collect();
    assert_eq!(h_df, g.iter().collect::<Vec<_>>(), "DoubletFinder's result");
    let clustering = h.iter().find(|(k, _)| k == "clustering").map(|(_, v)| v).expect("result.clustering");
    assert_eq!(Some(clustering), read_json(&hook.join("clusters/summary.json")).get("result"), "the clustering's numbers");
    assert!(clustering.get("purity").and_then(|p| p.get("computed")) == Some(&Json::Bool(true)));
    // the manifests: the nested output listed by its summary.json, every file described
    let v = scrust_in(&d, &["validate", "manifest", "--dir", hook.to_str().unwrap()]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    let m = scio::manifest::Manifest::read(&hook).unwrap();
    assert_eq!(m.entry("clusters/summary.json").map(|e| e.role.as_str()), Some(scio::manifest::NESTED_OUTPUT));
    assert!(m.entry("clusters/seurat_clusters.tsv").is_none(), "the nested output's files are its own manifest's");
    assert_eq!(m.provenance.inputs.iter().map(|i| i.option.as_str()).collect::<Vec<_>>(), ["counts"]);
    let inner = scio::manifest::Manifest::read(&hook.join("clusters")).unwrap();
    assert_eq!(inner.command, "cluster");
    assert_eq!(inner.entry("seurat_clusters.tsv").map(|e| e.role.as_str()), Some("cluster_labels"));
    assert!(inner.provenance.command.starts_with("scrust cluster --counts ") && inner.provenance.command.ends_with(" --preset workflow --stats silhouette,purity"), "{}", inner.provenance.command);
}

/// The labels clustering on the fly computes are R's, with R's levels: the workflow clustered
/// DoubletFinder's inputs with the same chain (K12, on tier T). BLAS-dependent: the fixtures'
/// kernel only.
#[test]
fn clustering_on_the_fly_gives_r_s_labels() {
    for object in ["soupx_tenx_test", "decontx_raw_emptydrops_test", "cellbender_fromraw_test"] {
        let f = fixture(object);
        rcompat::blas::require_fixture_kernel(&f, "clustering on the fly against R's labels");
        let d = tmp(&format!("on_the_fly_r_{object}"));
        ok(&d, &["doubletfinder", "--counts", &f.join("input").display().to_string(), "--out", "run", "--threads", "2"]);
        for file in ["seurat_clusters.tsv", "seurat_clusters_levels.txt"] {
            assert!(std::fs::read(d.join("run/clusters").join(file)).unwrap() == std::fs::read(f.join("input").join(file)).unwrap(), "{object}: {file}");
        }
        // and so DoubletFinder's files are those of the run given R's labels
        let reference = f.join("reference");
        for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
            if reference.join(file).is_file() {
                assert!(std::fs::read(d.join("run").join(file)).unwrap() == std::fs::read(reference.join(file)).unwrap(), "{object}: {file} differs from R's run");
            }
        }
    }
}

/// The singlets re-clustered in process are the two-command route's files, with and without the
/// export.
#[test]
fn reclustering_in_process_is_the_two_command_route() {
    let f = fixture("decontx_raw_emptydrops_test");
    let input = f.join("input").display().to_string();
    let clusters = f.join("input/seurat_clusters.tsv").display().to_string();
    for export in [false, true] {
        let d = tmp(if export { "recluster_export" } else { "recluster" });
        let (one, two) = (d.join("in_process"), d.join("two_commands"));
        std::fs::create_dir_all(&one).unwrap();
        std::fs::create_dir_all(&two).unwrap();
        let mut a = vec!["doubletfinder", "--counts", &input, "--clusters", &clusters, "--out", "run", "--threads", "2", "--recluster"];
        if export {
            a.push("--recluster-export");
        }
        ok(&one, &a);
        // the same thread count (run.threads is recorded); the files do not depend on it
        ok(&two, &["doubletfinder", "--counts", &input, "--clusters", &clusters, "--out", "run", "--threads", "3"]);
        let mut b = vec!["cluster", "--counts", "run/singlets", "--preset", "workflow", "--out", "run/recluster", "--threads", "2"];
        if export {
            b.push("--export");
        }
        ok(&two, &b);
        // recluster/ against the second command's output: every file, its summary.json included
        assert_same_tree(&one.join("run/recluster"), &two.join("run/recluster"), &[], &[], &[]);
        // DoubletFinder's own files too; its summary.json records the re-clustering
        assert_same_tree(&one.join("run"), &two.join("run"), &["recluster/"], &["recluster/"], &["summary.json"]);
        let s = read_json(&one.join("run/recluster/summary.json"));
        let p = s.get("files").and_then(|f| f.get("provenance")).unwrap();
        let cl: Vec<&str> = match p.get("command_line") {
            Some(Json::Vector(v) | Json::List(v)) => v.iter().filter_map(|x| x.as_str()).collect(),
            _ => panic!("no command_line"),
        };
        assert_eq!(&cl[..2], ["scrust", "doubletfinder"], "the command line as given");
        let want = format!("scrust cluster --counts run/singlets --preset workflow{}", if export { " --export" } else { "" });
        assert_eq!(p.get("command").and_then(|v| v.as_str()), Some(want.as_str()));
        let v = scrust_in(&one, &["validate", "manifest", "--dir", "run"]);
        assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
        let top = read_json(&one.join("run/summary.json"));
        assert_eq!(top.get("result").and_then(|r| r.get("recluster")).and_then(|r| r.get("n_cells")).and_then(|v| v.as_f64()), Some(scio::tenx::read_barcodes(&one.join("run/singlets")).unwrap().len() as f64));
    }
}

#[test]
fn the_new_options_refuse_what_they_cannot_do() {
    let f = fixture("soupx_tenx_test");
    let input = f.join("input").display().to_string();
    let clusters = f.join("input/seurat_clusters.tsv").display().to_string();
    let d = tmp("refusals");
    let refused = |extra: &[&str], needle: &str| {
        let mut a = vec!["doubletfinder", "--counts", input.as_str(), "--out", "x"];
        a.extend_from_slice(extra);
        let o = scrust_in(&d, &a);
        assert!(!o.status.success(), "{extra:?} succeeded");
        let e = String::from_utf8_lossy(&o.stderr);
        assert!(e.contains(needle), "{extra:?}: {e}");
    };
    refused(&["--clusters", &clusters, "--cluster-preset", "workflow"], "clustering on the fly is for runs without --clusters");
    refused(&["--clusters", &clusters, "--cluster-stats", "purity"], "clustering on the fly is for runs without --clusters");
    refused(&["--cluster-levels", "levels.txt"], "--cluster-levels goes with --clusters");
    refused(&["--from-stages", "stages"], "give --clusters");
    refused(&["--recluster-export"], "--recluster-export goes with --recluster");
    refused(&["--recluster-preset", "workflow"], "goes with --recluster");
    refused(&["--cluster-preset", "nonsense"], "unknown preset");
    refused(&["--cluster-stats", "entropy"], "entropy");
    refused(&["--recluster", "--recluster-params", "missing.toml"], "missing.toml");
    assert!(!d.join("x").exists(), "nothing written when refused");
}
