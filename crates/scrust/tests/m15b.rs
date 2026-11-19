//! Milestone M15b (docs/PLAN.md C3, "vars.to.regress and percent.mt"): SCTransform's
//! `vars.to.regress` and the workflow's `percent.mt`, end to end.
//!
//! - **Against R** (the fixtures' kernel only): `scrust cluster --counts` with the cases of
//!   validation/m15b_cases.tsv on two tier-T objects (a v5 object and the emptyDrops object's v3
//!   path) gives the labels and levels, the silhouette widths, the purities, `percent.mt`, the
//!   regressed column's bits, `scale.data` (R's digest) and the regression's rank and pivot of R's
//!   runs (validation/R/reference_clustering.R, portable profile; tests/fixtures/clustering/m15b/).
//!   On the test data, which has no mitochondrial gene, R's `percent.mt` by the pattern is 0 in
//!   every cell: R regresses it through `dqrdc2`'s rank-deficient path (rank 1), and scrust gives
//!   the same all-zero column through a list of features without counts.
//! - **The routes** (any kernel): `scrust sctransform` with the options, then `scrust pca`,
//!   `scrust neighbors` and the staged starts, write the files `scrust cluster --counts` writes
//!   with them, `percent.mt` carried in `cell_metadata.tsv`.
//! - **The in-memory rule:** the entry point called with the covariates in memory, its result
//!   written by the library's writer, gives the command's files.
//! - **In process:** `scrust doubletfinder` with the options in its parameter files (and
//!   `--recluster-metadata`) gives the two-command route's files, and DoubletFinder's own
//!   SCTransform does not regress.
//! - **The defaults:** without the options, nothing new is written; and the refusals.

use cluster::chain::{self, Options};
use cluster::chain_output::Writer;
use cluster::covariates::{CellCovariates, MtGenes};
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

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
}

/// The committed counts of a tier-T object.
fn input(object: &str) -> PathBuf {
    let df = fixtures().join("doubletfinder").join(object).join("input");
    if df.join("matrix.mtx.gz").is_file() {
        df
    } else {
        fixtures().join("clustering").join(object).join("input")
    }
}

fn m15b(object: &str) -> PathBuf {
    fixtures().join("clustering/m15b").join(object)
}

fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m15b").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn read_json(p: &Path) -> Json {
    Json::parse(&std::fs::read_to_string(p).unwrap()).unwrap()
}

fn table(text: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').map(String::from).collect();
    (header, lines.map(|l| l.split('\t').map(String::from).collect()).collect())
}

fn column(t: &(Vec<String>, Vec<Vec<String>>), name: &str) -> Option<Vec<String>> {
    let k = t.0.iter().position(|h| h == name)?;
    Some(t.1.iter().map(|r| r[k].clone()).collect())
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

/// `summary.json` without the manifest's entries (`files.entries`) and the files written
/// (`run.files`) under `dirs` (each ending in `/`).
fn without_entries(j: &Json, dirs: &[&str]) -> Json {
    let Json::Object(top) = j else { panic!("not an object") };
    let under = |p: Option<&str>| p.is_some_and(|p| dirs.iter().any(|d| p.starts_with(d)));
    let filtered = |s: &Json, key: &str| -> Json {
        let Json::Object(m) = s else { panic!("not an object") };
        Json::Object(
            m.iter()
                .map(|(k, v)| {
                    let keep = |xs: &Vec<Json>| -> Vec<Json> { xs.iter().filter(|x| !under(x.as_str().or_else(|| x.get("path").and_then(|p| p.as_str())))).cloned().collect() };
                    match v {
                        Json::List(xs) if k == key => (k.clone(), Json::List(keep(xs))),
                        Json::Vector(xs) if k == key => (k.clone(), Json::Vector(keep(xs))),
                        _ => (k.clone(), v.clone()),
                    }
                })
                .collect(),
        )
    };
    Json::Object(
        top.iter()
            .map(|(k, v)| match k.as_str() {
                "files" => (k.clone(), filtered(v, "entries")),
                "run" => (k.clone(), filtered(v, "files")),
                _ => (k.clone(), v.clone()),
            })
            .collect(),
    )
}

/// The same files under `a` and `b`, byte for byte, but `summary.json`, compared as JSON without
/// the timings and the command line (unless listed in `skip_summary`).
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

/// A case's options for scrust (the R run's, with an all-zero list where R's pattern finds no
/// feature: validation/m15b_cases.tsv's rust_mt_genes).
fn case_args(object: &str, case: &str) -> Vec<String> {
    let i = m15b(object).join("inputs");
    let p = |f: &str| i.join(f).display().to_string();
    let v = |s: &[&str]| s.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    match case {
        "pm" => v(&["--vars-to-regress", "percent.mt", "--mt-genes", &p("zero_genes.txt")]),
        "ids" => v(&["--vars-to-regress", "percent.mt", "--mt-genes", &p("mt_genes.txt")]),
        "ids_ncount" => v(&["--vars-to-regress", "percent.mt,nCount_RNA", "--mt-genes", &p("mt_genes.txt")]),
        "meta" => v(&["--vars-to-regress", "score", "--metadata", &p("metadata.tsv")]),
        _ => unreachable!(),
    }
}

/// R's runs with vars.to.regress (portable profile): the labels, levels, scores, percent.mt, the
/// regressed column's bits, scale.data's digest and the regression's rank and pivot.
#[test]
fn the_regression_against_r() {
    for object in ["soupx_tenx_test", "filtered_seurat_emptydrops_test"] {
        for case in ["pm", "ids", "ids_ncount", "meta"] {
            let r = m15b(object).join(case);
            rcompat::blas::require_fixture_kernel(&r, "vars.to.regress against R");
            let d = tmp(&format!("r_{object}_{case}"));
            let c = input(object).display().to_string();
            let mut a: Vec<String> = ["cluster", "--counts", &c, "--preset", "workflow", "--out", "run", "--threads", "3"].iter().map(|x| x.to_string()).collect();
            a.extend(case_args(object, case));
            ok(&d, &a.iter().map(String::as_str).collect::<Vec<_>>());
            let run = d.join("run");
            let what = format!("{object} {case}");
            let rm = table(&String::from_utf8(scio::gz::gunzip(&std::fs::read(r.join("metadata.tsv.gz")).unwrap()).unwrap()).unwrap());
            let cm = table(&std::fs::read_to_string(run.join("metadata.tsv")).unwrap());
            for col in ["nCount_RNA", "nFeature_RNA", "percent.mt", "SCT_snn_res.0.8", "seurat_clusters", "silhouette_width", "neighborhood_purity"] {
                match column(&rm, col) {
                    Some(want) => assert_eq!(column(&cm, col).as_ref(), Some(&want), "{what}: metadata.tsv {col}"),
                    None => assert!(col != "percent.mt" || case == "meta", "{what}: R's run has no {col}"),
                }
            }
            if case == "meta" {
                // R regressed the doubles read.delim() parsed; scrust carries the table's text and
                // parses it as R does
                let (want, got) = (column(&rm, "score").unwrap(), column(&cm, "score").unwrap());
                for (w, g) in want.iter().zip(&got) {
                    assert_eq!(scio::parse17(w).unwrap().to_bits(), rcompat::as_numeric(g).to_bits(), "{what}: score {g}");
                }
            }
            for f in ["seurat_clusters_levels.txt", "SCT_snn_res.0.8_levels.txt"] {
                assert_eq!(std::fs::read(r.join(f)).unwrap(), std::fs::read(run.join(f)).unwrap(), "{what}: {f}");
            }
            let (rs, cs) = (read_json(&r.join("summary.json")), read_json(&run.join("summary.json")));
            let (rr, cr) = (rs.get("result").unwrap(), cs.get("result").unwrap());
            assert_eq!(rr.get("scale_data_sha256"), cr.get("scale_data_sha256"), "{what}: scale.data (R's digest)");
            let (rg, cg) = (rr.get("regression").unwrap(), cr.get("sctransform").and_then(|s| s.get("regression")).unwrap());
            assert_eq!(rg.get("rank").and_then(Json::as_f64), cg.get("rank").and_then(Json::as_f64), "{what}: rank");
            let ints = |j: Option<&Json>| match j {
                Some(Json::Vector(v) | Json::List(v)) => v.iter().filter_map(Json::as_f64).collect::<Vec<_>>(),
                _ => Vec::new(),
            };
            assert_eq!(ints(rg.get("pivot")), ints(cg.get("pivot")), "{what}: pivot");
            if case == "pm" {
                // percent.mt 0 in every cell: lm()'s QR drops it (rank 1), as R's does
                assert_eq!(cg.get("rank").and_then(Json::as_f64), Some(1.0), "{what}: rank-deficient");
            }
        }
    }
}

/// `scrust sctransform` with the options, then the staged commands, write what `scrust cluster
/// --counts` writes with them.
#[test]
fn the_routes_with_the_regression() {
    let object = "soupx_tenx_test";
    let c = input(object).display().to_string();
    let i = m15b(object).join("inputs");
    let (mt, md) = (i.join("mt_genes.txt").display().to_string(), i.join("metadata.tsv").display().to_string());
    for (name, sct_opts, cluster_opts) in [
        ("pm_ncount", vec!["--vars-to-regress", "percent.mt,nCount_RNA", "--mt-genes", mt.as_str()], vec![]),
        ("score", vec!["--vars-to-regress", "score", "--percent-mt", "--metadata", md.as_str()], vec!["--metadata", md.as_str()]),
    ] {
        let d = tmp(&format!("routes_{name}"));
        let mut one: Vec<&str> = vec!["cluster", "--counts", &c, "--preset", "workflow", "--out", "one", "--keep-stages", "--threads", "2"];
        one.extend(&sct_opts);
        one.extend(&cluster_opts);
        if name == "score" {
            // --metadata once: for the regression and for metadata.tsv
            let k = one.iter().rposition(|a| *a == "--metadata").unwrap();
            one.drain(k..k + 2);
        }
        ok(&d, &one);
        let mut s: Vec<&str> = vec!["sctransform", "--counts", &c, "--out", "sct", "--seed", "12345", "--threads", "2"];
        s.extend(&sct_opts);
        ok(&d, &s);
        ok(&d, &["pca", "--sct", "sct", "--out", "pca", "--seed", "12345"]);
        ok(&d, &["neighbors", "--pca", "pca", "--out", "nn", "--dims", "1:30"]);
        for f in ["sct", "pca", "nn"] {
            let t = std::fs::read_to_string(d.join(f).join("cell_metadata.tsv")).unwrap();
            assert!(t.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tpercent.mt\n"), "{name}: {f}/cell_metadata.tsv carries percent.mt");
        }
        let mut a: Vec<&str> = vec!["cluster", "--sct", "sct", "--preset", "workflow", "--out", "from_sct", "--keep-stages"];
        a.extend(&cluster_opts);
        ok(&d, &a);
        let mut a: Vec<&str> = vec!["cluster", "--snn", "nn", "--pca", "pca", "--preset", "workflow", "--out", "from_snn"];
        a.extend(&cluster_opts);
        ok(&d, &a);
        let meta = std::fs::read_to_string(d.join("one/metadata.tsv")).unwrap();
        assert!(meta.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tpercent.mt\tSCT_snn_res.0.8"), "{name}: {}", &meta[..80]);
        let skip_neighbors = ["summary.json", "stages/", "sct/"];
        for (other, skip) in [("from_sct", &["summary.json", "stages/"][..]), ("from_snn", &skip_neighbors[..])] {
            let want: BTreeMap<String, PathBuf> = tree(&d.join("one")).into_iter().filter(|(k, _)| !skip.iter().any(|s| k == s || (s.ends_with('/') && k.starts_with(s)))).collect();
            let got = tree(&d.join(other));
            for (k, p) in &want {
                assert!(got.get(k).map(|q| std::fs::read(q).unwrap()) == Some(std::fs::read(p).unwrap()), "{name}: {other}: {k} differs or is missing");
            }
        }
        // SCTransform's result, its regression and percent.mt's record included, by every route
        let sct_of = |p: &str| read_json(&d.join(p).join("summary.json")).get("result").and_then(|r| r.get("sctransform")).cloned();
        let s1 = sct_of("one").unwrap();
        assert!(s1.get("regression").is_some() && s1.get("percent_mt").is_some(), "{name}");
        assert_eq!(read_json(&d.join("sct/summary.json")).get("result").cloned(), Some(s1.clone()), "{name}: scrust sctransform's result");
        assert_eq!(sct_of("from_sct"), Some(s1), "{name}: --sct copies it");
    }
}

/// The in-memory rule (A3) with the covariates handed over in memory: the entry point's result,
/// written by the library's writer, is `scrust cluster --counts`' output.
#[test]
fn the_entry_point_regresses_in_memory() {
    let object = "filtered_seurat_emptydrops_test";
    let inp = input(object);
    let i = m15b(object).join("inputs");
    let d = tmp("in_memory");
    let (c, mt, md) = (inp.display().to_string(), i.join("mt_genes.txt").display().to_string(), i.join("metadata.tsv").display().to_string());
    ok(&d, &["cluster", "--counts", &c, "--preset", "workflow", "--export", "--vars-to-regress", "percent.mt,score", "--mt-genes", &mt, "--metadata", &md, "--out", "command", "--threads", "2"]);
    let cmd_dir = d.join("command");
    let summary = read_json(&cmd_dir.join("summary.json"));
    let prov = summary.get("files").and_then(|f| f.get("provenance")).unwrap();
    // the command reproduces the values: --metadata stays when its columns are regressed
    let command = prov.get("command").and_then(Json::as_str).unwrap().to_string();
    assert!(command.contains("--vars-to-regress percent.mt,score") && command.contains("--metadata ") && command.contains("--mt-genes "), "{command}");
    let options: Vec<String> = match prov.get("inputs") {
        Some(Json::List(v)) => v.iter().filter_map(|x| x.get("option").and_then(Json::as_str).map(String::from)).collect(),
        _ => panic!("no inputs"),
    };
    assert_eq!(options, ["counts", "mt-genes", "metadata"], "each input once, digested");
    // in this process: the covariates in memory
    let params = cluster::chain_params::build(Some("workflow"), None, &[("vars-to-regress".into(), "percent.mt,score".into()), ("mt-genes".into(), mt.clone())]).unwrap();
    let counts = scio::Counts::read(&inp).unwrap();
    let sha = scio::sha256::file(Path::new(&mt)).unwrap().1;
    let cov = CellCovariates {
        mt_genes: Some(MtGenes::from_file_text(&std::fs::read_to_string(&mt).unwrap(), &mt, &sha)),
        metadata: Some(cluster::metadata::extra_columns(Path::new(&md), &counts.barcodes, &[]).unwrap()),
    };
    let opts = Options { correct: true, keep_reductions: true, keep_sct: true, ..Options::default() };
    let r = chain::cluster_counts_with(Cow::Owned(counts), &params, &cov, 2, opts, &mut chain::NoObserver).unwrap();
    assert!(r.notes.iter().any(|n| n.contains("not among the counts' features")), "{:?}", r.notes);
    let lib_dir = d.join("library");
    std::fs::create_dir_all(&lib_dir).unwrap();
    let mut w = Writer::new(&lib_dir, r.cells.clone(), 2);
    w.write_stages = false;
    w.export = true;
    w.command = command;
    let Some(Json::List(inputs)) = prov.get("inputs") else { panic!() };
    w.input_checksums = inputs
        .iter()
        .flat_map(|i| {
            let path = PathBuf::from(i.get("path").and_then(Json::as_str).unwrap());
            let files = match i.get("files") {
                Some(Json::List(f)) => f.clone(),
                _ => Vec::new(),
            };
            files.into_iter().map(move |f| (f.get("sha256").and_then(Json::as_str).unwrap().to_string(), path.join(f.get("name").and_then(Json::as_str).unwrap()).display().to_string()))
        })
        .collect();
    w.extra = Some(cluster::metadata::extra_columns(Path::new(&md), &r.cells, &["percent.mt"]).unwrap());
    w.write_clustered(&r).unwrap();
    w.finish_labels(&r.chain.clusters, r.chain.silhouette.as_ref(), r.chain.purity.as_ref()).unwrap();
    for (rel, p) in &tree(&lib_dir) {
        assert!(std::fs::read(p).unwrap() == std::fs::read(cmd_dir.join(rel)).unwrap(), "{rel} differs");
    }
    let as_text = |j: &Json| Json::parse(&scio::json::to_file_text(j)).unwrap();
    assert_eq!(&as_text(&cluster::chain_output::clustered_result_json(&r)), summary.get("result").unwrap(), "summary.json's result");
}

/// `scrust doubletfinder` with the options in its parameter files: clustering on the fly is
/// `scrust cluster --counts` with them, the re-clustering the two-command route's, DoubletFinder's
/// own files those of a run without them.
#[test]
fn doubletfinder_in_process_with_the_regression() {
    let object = "soupx_tenx_test";
    let inp = input(object).display().to_string();
    let clusters = input(object).join("seurat_clusters.tsv").display().to_string();
    let i = m15b(object).join("inputs");
    let d = tmp("doubletfinder");
    // parameter files next to the id list, which they name by a relative path
    std::fs::copy(i.join("mt_genes.txt"), d.join("mt_genes.txt")).unwrap();
    std::fs::copy(i.join("metadata.tsv"), d.join("metadata.tsv")).unwrap();
    std::fs::write(d.join("hook.toml"), "[sctransform]\nvars_to_regress = [\"percent.mt\", \"nCount_RNA\"]\nmt_genes = \"mt_genes.txt\"\n").unwrap();
    std::fs::write(d.join("re.toml"), "[sctransform]\nvars_to_regress = [\"percent.mt\", \"score\"]\nmt_genes = \"mt_genes.txt\"\n").unwrap();
    for sub in ["hook", "alone", "one", "two", "plain"] {
        std::fs::create_dir_all(d.join(sub)).unwrap();
    }
    let up = |f: &str| format!("../{f}");
    // clustering on the fly against scrust cluster --counts with the same file
    ok(&d.join("hook"), &["doubletfinder", "--counts", &inp, "--out", "run", "--cluster-params", &up("hook.toml"), "--cluster-stats", "silhouette,purity", "--threads", "2"]);
    ok(&d.join("alone"), &["cluster", "--counts", &inp, "--preset", "workflow", "--params", &up("hook.toml"), "--stats", "silhouette,purity", "--out", "run/clusters", "--threads", "2"]);
    let hook_meta = std::fs::read_to_string(d.join("hook/run/clusters/metadata.tsv")).unwrap();
    assert!(hook_meta.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tpercent.mt\t"), "{}", &hook_meta[..60]);
    assert_same_tree(&d.join("hook/run/clusters"), &d.join("alone/run/clusters"), &[], &["pca/", "neighbors/"], &["summary.json"]);
    let theirs = without_entries(&read_json(&d.join("alone/run/clusters/summary.json")), &["pca/", "neighbors/"]);
    assert_eq!(settled(&read_json(&d.join("hook/run/clusters/summary.json")), ""), settled(&theirs, ""), "clusters/summary.json");
    // DoubletFinder's own SCTransform does not regress: its files are a run given the hook's labels
    let labels = d.join("hook/run/clusters/seurat_clusters.tsv").display().to_string();
    ok(&d.join("plain"), &["doubletfinder", "--counts", &inp, "--clusters", &labels, "--out", "run", "--threads", "2"]);
    for f in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv", "singlets/matrix.mtx.gz"] {
        assert!(std::fs::read(d.join("hook/run").join(f)).unwrap() == std::fs::read(d.join("plain/run").join(f)).unwrap(), "{f}");
    }
    // the re-clustering in process, a --metadata column regressed, against the two commands
    ok(&d.join("one"), &["doubletfinder", "--counts", &inp, "--clusters", &clusters, "--out", "run", "--recluster", "--recluster-params", &up("re.toml"), "--recluster-metadata", &up("metadata.tsv"), "--threads", "2"]);
    ok(&d.join("two"), &["doubletfinder", "--counts", &inp, "--clusters", &clusters, "--out", "run", "--threads", "3"]);
    ok(&d.join("two"), &["cluster", "--counts", "run/singlets", "--preset", "workflow", "--params", &up("re.toml"), "--metadata", &up("metadata.tsv"), "--out", "run/singlets_clusters", "--threads", "2"]);
    assert_same_tree(&d.join("one/run/singlets_clusters"), &d.join("two/run/singlets_clusters"), &[], &[], &[]);
    assert_same_tree(&d.join("one/run"), &d.join("two/run"), &["singlets_clusters/"], &["singlets_clusters/"], &["summary.json"]);
    let s = read_json(&d.join("one/run/singlets_clusters/summary.json"));
    let cmd = s.get("files").and_then(|f| f.get("provenance")).and_then(|p| p.get("command")).and_then(Json::as_str).unwrap().to_string();
    assert_eq!(cmd, "scrust cluster --counts run/singlets --preset workflow --params ../re.toml --metadata ../metadata.tsv", "the command whose values it reproduces");
    let meta = std::fs::read_to_string(d.join("one/run/singlets_clusters/metadata.tsv")).unwrap();
    assert!(meta.starts_with("barcode\tnCount_RNA\tnFeature_RNA\tpercent.mt\t") && meta.lines().next().unwrap().ends_with("\tscore\tlabel"), "{}", meta.lines().next().unwrap());
    for dir in ["hook/run", "one/run"] {
        let v = scrust_in(&d, &["validate", "manifest", "--dir", dir]);
        assert!(v.status.success(), "{dir}: {}", String::from_utf8_lossy(&v.stderr));
    }
}

/// Without the new options nothing new is written; with --percent-mt alone only percent.mt is
/// added; and what is refused.
#[test]
fn the_defaults_and_the_refusals() {
    let object = "soupx_tenx_test";
    let c = input(object).display().to_string();
    let i = m15b(object).join("inputs");
    let d = tmp("defaults");
    ok(&d, &["cluster", "--counts", &c, "--preset", "workflow", "--out", "plain", "--threads", "2"]);
    ok(&d, &["cluster", "--counts", &c, "--preset", "workflow", "--out", "empty", "--vars-to-regress", "", "--threads", "2"]);
    let o = ok(&d, &["cluster", "--counts", &c, "--preset", "workflow", "--out", "pm", "--percent-mt", "--threads", "2"]);
    assert!(String::from_utf8_lossy(&o.stderr).contains("percent.mt is 0 in every cell; give the mitochondrial genes' ids with --mt-genes"));
    let plain = read_json(&d.join("plain/summary.json"));
    let text = std::fs::read_to_string(d.join("plain/summary.json")).unwrap();
    for key in ["\"regression\"", "\"percent_mt\"", "\"notes\"", "\"vars_to_regress\"", "\"mt_genes\""] {
        assert!(!text.contains(key), "a run without the options records {key}");
    }
    assert!(std::fs::read_to_string(d.join("plain/metadata.tsv")).unwrap().starts_with("barcode\tnCount_RNA\tnFeature_RNA\tSCT_snn_res.0.8\t"));
    // an empty list is no regression: the same files, and the same result
    assert_same_tree(&d.join("plain"), &d.join("empty"), &[], &[], &["summary.json"]);
    assert_eq!(plain.get("result"), read_json(&d.join("empty/summary.json")).get("result"));
    // percent.mt alone: the scale.data and every file but metadata.tsv (its column) are the same
    assert_same_tree(&d.join("plain"), &d.join("pm"), &["metadata.tsv", "summary.json"], &["metadata.tsv", "summary.json"], &[]);
    let pm = read_json(&d.join("pm/summary.json"));
    assert_eq!(plain.get("result").and_then(|r| r.get("scale_data_sha256")), pm.get("result").and_then(|r| r.get("scale_data_sha256")));
    let rec = pm.get("result").and_then(|r| r.get("sctransform")).and_then(|s| s.get("percent_mt")).unwrap();
    assert_eq!(rec.get("definition").and_then(Json::as_str), Some("pattern"));
    assert_eq!(rec.get("n_features").and_then(Json::as_f64), Some(0.0));
    let t = table(&std::fs::read_to_string(d.join("pm/metadata.tsv")).unwrap());
    assert!(column(&t, "percent.mt").unwrap().iter().all(|v| v == "0"), "zeros, as computed");
    // refusals
    let md = i.join("metadata.tsv").display().to_string();
    let mt = i.join("mt_genes.txt").display().to_string();
    let bad_md = d.join("bad.tsv");
    std::fs::write(&bad_md, "barcode\tgappy\n").unwrap();
    let feature = scio::tenx::read_features(&input(object)).unwrap().into_iter().map(|f| f.name).find(|n| cluster::covariates::is_syntactic(n)).unwrap();
    let nope = d.join("nope.txt");
    std::fs::write(&nope, "NOT_HERE\n").unwrap();
    let refused = |args: &[&str], needle: &str| {
        let o = scrust_in(&d, args);
        assert!(!o.status.success(), "{args:?} accepted");
        let e = String::from_utf8_lossy(&o.stderr);
        assert!(e.contains(needle), "{args:?}: {e}");
    };
    let cl = |extra: &[&str]| -> Vec<String> {
        let mut a: Vec<String> = ["cluster", "--counts", &c, "--preset", "workflow", "--out", "refused"].iter().map(|x| x.to_string()).collect();
        a.extend(extra.iter().map(|x| x.to_string()));
        a
    };
    fn as_str(v: &[String]) -> Vec<&str> {
        v.iter().map(String::as_str).collect()
    }
    refused(&as_str(&cl(&["--vars-to-regress", "percent.mt"])), "--mt-genes");
    refused(&as_str(&cl(&["--vars-to-regress", "nCount_RNA,nCount_RNA"])), "given twice");
    refused(&as_str(&cl(&["--vars-to-regress", "my-col"])), "not a syntactic name");
    refused(&as_str(&cl(&["--vars-to-regress", "GENE"])), "response");
    refused(&as_str(&cl(&["--vars-to-regress", &feature])), "no --metadata table");
    refused(&as_str(&cl(&["--vars-to-regress", "score"])), "no --metadata table is given");
    refused(&as_str(&cl(&["--vars-to-regress", "label", "--metadata", &md])), "non-numeric");
    refused(&as_str(&cl(&["--vars-to-regress", "nope", "--metadata", &md])), "not percent.mt, nCount_RNA, nFeature_RNA or a column of --metadata");
    refused(&as_str(&cl(&["--vars-to-regress", "gappy", "--metadata", bad_md.to_str().unwrap()])), "missing");
    refused(&as_str(&cl(&["--percent-mt", "--mt-genes", nope.to_str().unwrap()])), "none of its 1 id(s)");
    refused(&["sctransform", "--counts", &c, "--out", "refused_sct", "--metadata", &md], "reads it only for the columns --vars-to-regress names");
    refused(&["cluster", "--sct", "x", "--out", "refused_sct2", "--vars-to-regress", "nCount_RNA"], "act where SCTransform runs");
    refused(&["cluster", "--snn", "x", "--out", "refused_snn", "--mt-genes", &mt], "act where SCTransform runs");
    refused(&["doubletfinder", "--counts", &c, "--clusters", "x.tsv", "--out", "refused_df", "--cluster-metadata", &md], "clustering on the fly is for runs without --clusters");
    refused(&["doubletfinder", "--counts", &c, "--out", "refused_df2", "--recluster-metadata", &md], "goes with --recluster");
    // a column named like a feature, from the table: Seurat would add the feature's residuals too
    let named = d.join("named.tsv");
    let bc = scio::tenx::read_barcodes(&input(object)).unwrap();
    std::fs::write(&named, format!("barcode\t{feature}\n{}", bc.iter().enumerate().map(|(k, b)| format!("{b}\t{k}\n")).collect::<String>())).unwrap();
    refused(&as_str(&cl(&["--vars-to-regress", &feature, "--metadata", named.to_str().unwrap()])), "also a feature's name");
}
