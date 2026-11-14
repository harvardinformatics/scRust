//! `scrust neighbors`: Seurat 5.5.1's `FindNeighbors` on a PCA (docs/PLAN.md C3, C3a, C3b;
//! milestone M9): the Annoy forest, every cell's neighbours, `SCT_nn` and `SCT_snn`.

use crate::cli::{self, Args, Phases};
use neighbors::{parse_dims, parse_number, Metric, NeighborParams};
use scio::{fmt17, parse17, Json, Npy, Table};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: scrust neighbors --pca <dir> --out <dir> [--dims 1:10] [--k 20] [--trees 50]
                        [--metric euclidean|cosine|manhattan] [--prune 1/15] [--search-k -1]
                        [--annoy yes|no] [--threads N]

Seurat 5.5.1's FindNeighbors(object, dims, k.param, prune.SNN, n.trees, annoy.metric) on cell
embeddings, as R computes it, bit for bit: the embedding columns rounded to float, the Annoy
forest RcppAnnoy 0.0.23 builds (Annoy 1.17.3, Kiss64Random from its default seed), each cell's
k nearest neighbours (search.k = -1: k x n.trees candidates), the 0/1 kNN graph SCT_nn and the
shared-nearest-neighbour graph SCT_snn (ComputeSNN: s / (k + (k - s)) for s shared neighbours,
entries below prune.SNN removed). No BLAS is involved: the results do not depend on the machine
or the thread count. The workflow calls FindNeighbors(obj, dims = 1:30).

  --pca <dir>        cell embeddings, either `scrust pca`'s embeddings.tsv (barcode, PC_1 ..),
                     or embeddings.npy (cells x components, as the R reference runs write it)
                     with the barcodes in barcodes.tsv beside it or in ../metadata.tsv (either
                     may be gzipped)
  --out <dir>        output directory (must not hold a previous run)
  --dims D           the embedding columns used, 1-based: a:b, numbers, comma-separated lists
                     (default 1:10, Seurat's; the workflow uses 1:30)
  --k N              k.param (default 20); with fewer cells, the number of cells minus 1
  --trees N          n.trees (default 50)
  --metric M         annoy.metric: euclidean (default), cosine (Annoy's angular) or manhattan
  --prune X          prune.SNN (default 1/15; a fraction a/b is divided, as R evaluates it)
  --search-k N       AnnoySearch's search.k (default -1, the only value FindNeighbors uses)
  --annoy yes|no     also write annoy.ann, the index in Annoy's file format (default no)
  --threads N        threads for the search and the graphs (default: SLURM_CPUS_PER_TASK,
                     else 1); the build is sequential, as in RcppAnnoy

Outputs:
  knn_index.tsv      barcode, nn_1 .. nn_k: each cell's neighbours, 1-based, nearest first
  knn_distance.tsv   barcode, nn_1 .. nn_k: their distances as Annoy returns them (floats,
                     written with 17 significant digits, so exactly)
  nn.mtx.gz          SCT_nn (Matrix Market, cells x cells, entries in stored order)
  snn.mtx.gz         SCT_snn (Matrix Market, values with 17 significant digits)
  annoy.ann          with --annoy yes
  summary.json       parameters, sizes, the build's statistics, timings, peak memory
";

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// Cell embeddings: barcodes and the matrix (cells x components, column-major).
pub struct Embeddings {
    pub barcodes: Vec<String>,
    pub n: usize,
    pub ncol: usize,
    pub data: Vec<f64>,
    pub source: PathBuf,
}

/// Reads `scrust pca`'s `embeddings.tsv`, or an R run's `embeddings.npy` with its barcodes.
pub fn read_embeddings(dir: &Path) -> Result<Embeddings, String> {
    let tsv = dir.join("embeddings.tsv");
    if tsv.exists() {
        let t = Table::read(&tsv).map_err(|e| format!("{}: {e}", tsv.display()))?;
        let n = t.rows.len();
        let ncol = t.header.len().saturating_sub(1);
        if ncol == 0 {
            return Err(format!("{}: no components", tsv.display()));
        }
        let mut data = vec![0.0; n * ncol];
        for (i, row) in t.rows.iter().enumerate() {
            if row.len() != ncol + 1 {
                return Err(format!("{}: row {} has {} fields", tsv.display(), i + 1, row.len()));
            }
            for c in 0..ncol {
                data[c * n + i] = parse17(&row[c + 1]).ok_or_else(|| format!("{}: row {}: bad value '{}'", tsv.display(), i + 1, row[c + 1]))?;
            }
        }
        let barcodes = t.rows.iter().map(|r| r[0].clone()).collect();
        return Ok(Embeddings { barcodes, n, ncol, data, source: tsv });
    }
    let npy = dir.join("embeddings.npy");
    if npy.exists() {
        let a = Npy::read(&npy).map_err(|e| e.to_string())?;
        let (n, ncol) = (a.nrow(), a.ncol());
        let data = a.block_f64(0, n);
        let barcodes = if let Ok(p) = scio::tenx::find(dir, "barcodes.tsv") {
            scio::tenx::read_lines(&p).map_err(|e| e.to_string())?
        } else {
            let meta = scio::tenx::find(&dir.join(".."), "metadata.tsv").map_err(|e| format!("{e} (the barcodes of embeddings.npy: barcodes.tsv beside it, or ../metadata.tsv)"))?;
            let text = String::from_utf8(scio::tenx::read_maybe_gz(&meta).map_err(|e| e.to_string())?).map_err(|_| format!("{}: not text", meta.display()))?;
            let t = Table::parse(&text).ok_or_else(|| format!("{}: empty", meta.display()))?;
            let c = t.col("barcode").ok_or_else(|| format!("{}: no barcode column", meta.display()))?;
            t.rows.iter().map(|r| r[c].clone()).collect()
        };
        if barcodes.len() != n {
            return Err(format!("{}: {n} cells but {} barcodes", npy.display(), barcodes.len()));
        }
        return Ok(Embeddings { barcodes, n, ncol, data, source: npy });
    }
    Err(format!("{}: neither embeddings.tsv nor embeddings.npy", dir.display()))
}

/// Parses the FindNeighbors options shared by `scrust neighbors` (and, later, `scrust cluster`).
pub fn params_from_args(args: &Args) -> Result<NeighborParams, String> {
    let mut p = NeighborParams::default();
    if let Some(d) = args.get("dims") {
        p.dims = parse_dims(d)?;
    }
    p.k_param = args.number("k", p.k_param)?;
    p.n_trees = args.number("trees", p.n_trees)?;
    if let Some(m) = args.get("metric") {
        p.metric = Metric::parse(m)?;
    }
    if let Some(x) = args.get("prune") {
        p.prune_snn = parse_number(x).map_err(|e| format!("--prune: {e}"))?;
    }
    p.search_k = args.number("search-k", p.search_k)?;
    Ok(p)
}

fn knn_table(barcodes: &[String], k: usize, value: impl Fn(usize, usize) -> String) -> String {
    let mut t = String::from("barcode");
    for r in 1..=k {
        let _ = write!(t, "\tnn_{r}");
    }
    t.push('\n');
    for (i, b) in barcodes.iter().enumerate() {
        t.push_str(b);
        for r in 0..k {
            t.push('\t');
            t.push_str(&value(i, r));
        }
        t.push('\n');
    }
    t
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse(argv)?;
    args.only(&["pca", "out", "dims", "k", "trees", "metric", "prune", "search-k", "annoy", "threads"])?;
    let pca = args.path("pca")?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let params = params_from_args(&args)?;
    let write_annoy = match args.get("annoy").unwrap_or("no") {
        "yes" => true,
        "no" => false,
        v => return Err(format!("--annoy: expected yes or no, not '{v}'")),
    };
    if out.join("summary.json").exists() || out.join("snn.mtx.gz").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }
    let emb = read_embeddings(&pca)?;
    phases.mark("load");

    let r = neighbors::find_neighbors(&emb.data, emb.n, emb.ncol, &params, threads).map_err(|e| e.to_string())?;
    phases.mark("neighbors");

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let k = r.k;
    write(&out.join("knn_index.tsv"), &knn_table(&emb.barcodes, k, |i, j| (r.nn_idx[i * k + j] + 1).to_string()))?;
    write(&out.join("knn_distance.tsv"), &knn_table(&emb.barcodes, k, |i, j| fmt17(r.nn_dist[i * k + j] as f64)))?;
    let version = env!("CARGO_PKG_VERSION");
    let describe = |what: &str| {
        vec![
            format!(" scrust {version} neighbors: {what}"),
            format!(" rows and columns are the cells, in the order of knn_index.tsv ({} cells)", emb.n),
        ]
    };
    scio::mtx::write_gz(&out.join("nn.mtx.gz"), &r.nn, &describe("Seurat's SCT_nn graph (FindNeighbors): 1 at (cell, neighbour)"), threads).map_err(|e| e.to_string())?;
    scio::mtx::write_gz(&out.join("snn.mtx.gz"), &r.snn, &describe("Seurat's SCT_snn graph (FindNeighbors, ComputeSNN): s / (k + (k - s)) for s shared neighbours"), threads).map_err(|e| e.to_string())?;
    if write_annoy {
        r.index.save(&out.join("annoy.ann")).map_err(|e| format!("{}: {e}", out.join("annoy.ann").display()))?;
    }
    phases.mark("write");

    let self_first = (0..r.n).filter(|&i| k > 0 && r.nn_idx[i * k] == i as i32).count();
    let self_missing = (0..r.n).filter(|&i| !r.nn_idx[i * k..(i + 1) * k].contains(&(i as i32))).count();
    let mut distinct: Vec<u64> = r.snn.x.iter().map(|v| v.to_bits()).collect();
    distinct.sort_unstable();
    distinct.dedup();
    let st = r.index.stats();
    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let int = |v: usize| Json::Int(v as i64);
    let summary = Json::obj(vec![
        ("tool", Json::Str("scrust neighbors".into())),
        ("reference", Json::Str(neighbors::REFERENCE.into())),
        ("scrust_version", Json::Str(version.into())),
        (
            "parameters",
            Json::obj(vec![
                ("dims", Json::Vector(params.dims.iter().map(|&d| int(d)).collect())),
                ("k_param", int(params.k_param)),
                ("n_trees", Json::Int(params.n_trees as i64)),
                ("metric", Json::Str(params.metric.seurat_name().into())),
                ("annoy_index", Json::Str(format!("Annoy{}", match params.metric { Metric::Euclidean => "Euclidean", Metric::Angular => "Angular", Metric::Manhattan => "Manhattan" }))),
                ("prune_snn", Json::Num(params.prune_snn)),
                ("search_k", Json::Int(params.search_k as i64)),
                ("annoy_seed", Json::Str(neighbors::Kiss64Random::DEFAULT_SEED.to_string())),
            ]),
        ),
        (
            "result",
            Json::obj(vec![
                ("n_cells", int(r.n)),
                ("n_components", int(emb.ncol)),
                ("k", int(k)),
                ("annoy_f", int(r.index.f())),
                ("annoy_node_bytes", int(r.index.node_bytes())),
                ("annoy_n_nodes", int(r.index.n_nodes())),
                ("annoy_n_trees", int(r.index.n_trees())),
                ("annoy_file_bytes", int(r.index.n_nodes() * r.index.node_bytes())),
                (
                    "annoy_build",
                    Json::obj(vec![
                        ("split_nodes", int(st.split_nodes)),
                        ("leaves", int(st.leaves)),
                        ("accepted_at_attempt", Json::Vector(st.accepted_at_attempt.iter().map(|&v| int(v)).collect())),
                        ("kept_imbalanced", int(st.kept_imbalanced)),
                        ("random_fallbacks", int(st.random_fallbacks)),
                        ("random_rounds", int(st.random_rounds)),
                        ("zero_margins", int(st.zero_margins)),
                        ("draws", Json::Int(st.draws as i64)),
                        ("max_depth", int(st.max_depth)),
                    ]),
                ),
                ("knn_self_first", int(self_first)),
                ("knn_self_missing", int(self_missing)),
                ("nn_nnz", int(r.nn.nnz())),
                ("snn_nnz", int(r.snn.nnz())),
                ("snn_distinct_values", int(distinct.len())),
            ]),
        ),
        (
            "run",
            Json::obj(vec![
                ("pca", Json::Str(emb.source.display().to_string())),
                ("outdir", Json::Str(out.display().to_string())),
                ("threads", int(threads)),
                ("start_time", Json::Str(start_time)),
                ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
                ("build_s", Json::Num(r.timings.build)),
                ("search_s", Json::Num(r.timings.search)),
                ("graphs_s", Json::Num(r.timings.graphs)),
                ("elapsed_s", Json::Num(phases.elapsed())),
                ("phases", phases.to_json()),
                ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
                ("fingerprint", fingerprint),
            ]),
        ),
    ]);
    write(&out.join("summary.json"), &scio::json::to_file_text(&summary))?;
    eprintln!(
        "scrust neighbors: {} cells, {} dims, k = {k}, {} trees ({}): {} nodes; nn {} entries, snn {}; build {:.2} s, search {:.2} s, graphs {:.2} s; {threads} threads",
        r.n,
        params.dims.len(),
        params.n_trees,
        params.metric,
        r.index.n_nodes(),
        r.nn.nnz(),
        r.snn.nnz(),
        r.timings.build,
        r.timings.search,
        r.timings.graphs
    );
    Ok(())
}
