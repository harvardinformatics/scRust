//! `scrust cluster`: Seurat 5.5.1's `FindClusters` on an SNN graph (docs/PLAN.md C3, C3c, C3d;
//! milestone M10): the modularity optimizer, GroupSingletons, the columns and their levels.

use crate::cli::{self, Args, Phases};
use cluster::{parse_resolutions, ClusterParams};
use scio::Json;
use std::path::Path;

pub const USAGE: &str = "\
Usage: scrust cluster --snn <dir> --out <dir> [--resolution 0.8[,1.2,...]] [--algorithm 1|2|3]
                      [--starts 10] [--iterations 10] [--seed 0] [--modularity 1|2]
                      [--no-group-singletons] [--graph-name SCT_snn]
                      [--stats silhouette,purity|none] [--keep-stages] [--threads N]

Seurat 5.5.1's FindClusters(object, graph.name, modularity.fxn, resolution, algorithm, n.start,
n.iter, random.seed, group.singletons) on a shared-nearest-neighbour graph, as R computes it,
bit for bit: the network of the graph's lower triangle, the modularity optimizer
(RunModularityClusteringCpp: Seurat's C++ translation of ModularityOptimizer 1.3.0, with its
Java-style generator seeded from random.seed), the best of the random starts with its
clusters numbered by size, GroupSingletons (ties drawn by R's sample() after set.seed(1)),
and the columns and levels FindClusters.Seurat makes. No BLAS is involved: the results do not
depend on the machine or the thread count. The workflow calls FindClusters(obj,
random.seed = 12345) with every other default.

  --snn <dir>            the SNN graph: snn.mtx(.gz) as `scrust neighbors` writes it (with its
                         knn_index.tsv), or an R run's graphs/snn.mtx.gz (with ../metadata.tsv);
                         a directory holding graphs/ or neighbors/ is searched too
  --out <dir>            output directory (must not hold a previous run)
  --resolution R         resolution(s), comma-separated, each read as R reads it (default 0.8);
                         every resolution starts the optimizer afresh, as in Seurat
  --algorithm A          1: original Louvain (default); 2: Louvain with multilevel refinement;
                         3: SLM. 4 (Leiden) is not ported
  --starts N             n.start, random starts (default 10)
  --iterations N         n.iter, iterations per start (default 10)
  --seed N               random.seed (default 0, Seurat's; the workflow uses 12345)
  --modularity M         modularity.fxn: 1 standard (default), 2 alternative
  --no-group-singletons  group.singletons = FALSE: label single-cell clusters \"singleton\"
  --graph-name NAME      graph.name, the columns' prefix (default SCT_snn)
  --stats S              silhouette,purity (default, as in the workflow) or none. Silhouette
                         (milestone M11) and purity (M12) are not computed yet: when asked
                         for, the command says so and writes the labels only
  --keep-stages          also write stages/<column>_partitions.npy: every random start's own
                         partition (starts x cells), as the instrumented optimizer prints them
  --threads N            threads; resolutions run in parallel (default: SLURM_CPUS_PER_TASK,
                         else 1). Within a resolution the starts share one generator and run
                         in turn

Outputs:
  seurat_clusters.tsv         barcode, seurat_clusters (the last resolution's labels), the
                              layout `scrust doubletfinder --clusters` reads
  seurat_clusters_levels.txt  its levels, in order
  clusters.tsv                barcode and one column per resolution (SCT_snn_res.0.8, ...)
  <column>_levels.txt         each column's levels
  cluster_ids.txt, nclusters.txt   as the workflow's write_cluster_metadata writes them
  metadata.tsv                barcode, the columns and seurat_clusters
  optimizer_ids.tsv           each resolution's optimizer result before GroupSingletons
                              (0-based, numbered by size)
  summary.json                result: the network, every start's iterations and modularity, the
                              best start, what FindClusters printed, the singletons, the levels,
                              R's .Random.seed after GroupSingletons; run: parameters, timings,
                              peak memory, the environment
";

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// The clustering options shared by `scrust cluster`'s starting points.
pub fn params_from_args(args: &Args) -> Result<ClusterParams, String> {
    let mut p = ClusterParams::default();
    if let Some(r) = args.get("resolution") {
        p.resolutions = parse_resolutions(r)?;
    }
    let int = |name: &str, default: i32| -> Result<i32, String> {
        match args.get(name) {
            None => Ok(default),
            Some(v) => cluster::params::parse_int(&format!("--{name}"), v),
        }
    };
    p.algorithm = int("algorithm", p.algorithm)?;
    p.n_start = int("starts", p.n_start)?;
    p.n_iter = int("iterations", p.n_iter)?;
    p.random_seed = int("seed", p.random_seed)?;
    p.modularity_fxn = int("modularity", p.modularity_fxn)?;
    p.group_singletons = !args.flag("no-group-singletons");
    if let Some(g) = args.get("graph-name") {
        p.graph_name = g.to_string();
    }
    p.check()?;
    Ok(p)
}

/// `--stats`: which scores are asked for.
fn stats_from_args(args: &Args) -> Result<Vec<&'static str>, String> {
    let v = args.get("stats").unwrap_or("silhouette,purity");
    if v == "none" {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    for s in v.split(',') {
        match s.trim() {
            "silhouette" => out.push("silhouette"),
            "purity" => out.push("purity"),
            other => return Err(format!("--stats: '{other}' is not silhouette, purity or none")),
        }
    }
    Ok(out)
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse_with_flags(argv, &["no-group-singletons", "keep-stages"])?;
    for other in ["counts", "sct", "pca", "preset", "params", "export", "metadata", "purity-k", "purity-weighted"] {
        if args.get(other).is_some() {
            return Err(format!(
                "--{other}: this version clusters an SNN graph (--snn); starting from the counts or another stage, the preset, the TOML file and the scores' options come with milestones M11 to M13"
            ));
        }
    }
    args.only(&["snn", "out", "resolution", "algorithm", "starts", "iterations", "seed", "modularity", "no-group-singletons", "graph-name", "stats", "keep-stages", "threads"])?;
    let snn_dir = args.path("snn")?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let params = params_from_args(&args)?;
    let stats = stats_from_args(&args)?;
    let keep = args.flag("keep-stages");
    if out.join("summary.json").exists() || out.join("clusters.tsv").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }
    let input = cluster::input::read_snn_dir(&snn_dir)?;
    phases.mark("load");

    let r = cluster::find_clusters(&input.snn, &params, threads, keep).map_err(|e| e.to_string())?;
    phases.mark("cluster");
    // what FindClusters prints, as Seurat prints it
    for x in &r.resolutions {
        for line in x.printed() {
            eprintln!("{line}");
        }
        if let Some(m) = x.message() {
            eprintln!("{m}");
        }
    }

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let mut files = cluster::output::write_labels(&out, &input.barcodes, &r).map_err(|e| e.to_string())?;
    if keep {
        let stages = out.join("stages");
        std::fs::create_dir_all(&stages).map_err(|e| format!("{}: {e}", stages.display()))?;
        files.extend(cluster::output::write_partitions(&stages, &r).map_err(|e| e.to_string())?.into_iter().map(|f| format!("stages/{f}")));
    }
    phases.mark("write");

    let stats_note = if stats.is_empty() {
        "none asked for".to_string()
    } else {
        let msg = format!(
            "{} asked for, but not computed: this version writes the labels only (silhouette_width comes with milestone M11, neighborhood_purity with M12)",
            stats.join(" and ")
        );
        eprintln!("scrust cluster: --stats {}: {msg}", stats.join(","));
        msg
    };
    let version = env!("CARGO_PKG_VERSION");
    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let summary = Json::obj(vec![
        ("tool", Json::Str("scrust cluster".into())),
        ("reference", Json::Str(cluster::REFERENCE.into())),
        ("scrust_version", Json::Str(version.into())),
        ("result", cluster::output::result_json(&r)),
        (
            "run",
            Json::obj(vec![
                ("parameters", cluster::output::params_json(&params)),
                ("start", Json::Str("snn".into())),
                ("snn", Json::Str(input.source.display().to_string())),
                ("barcodes", Json::Str(input.barcodes_source.display().to_string())),
                ("outdir", Json::Str(out.display().to_string())),
                ("files", Json::strs(&files)),
                (
                    "stats",
                    Json::obj(vec![
                        ("requested", Json::strs(&stats)),
                        ("computed", Json::Vector(Vec::new())),
                        ("note", Json::Str(stats_note)),
                    ]),
                ),
                ("keep_stages", Json::Bool(keep)),
                ("threads", Json::Int(threads as i64)),
                ("start_time", Json::Str(start_time)),
                ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
                ("network_s", Json::Num(r.network_seconds)),
                ("resolutions_s", Json::Num(r.resolutions_seconds)),
                ("optimizer_s", Json::nums(&r.resolutions.iter().map(|x| x.optimizer.seconds).collect::<Vec<_>>())),
                ("groupsingletons_s", Json::nums(&r.resolutions.iter().map(|x| x.group_seconds).collect::<Vec<_>>())),
                ("printed_elapsed", Json::strs(&r.resolutions.iter().map(|x| x.printed().pop().unwrap_or_default()).collect::<Vec<_>>())),
                ("elapsed_s", Json::Num(phases.elapsed())),
                ("phases", phases.to_json()),
                ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
                ("fingerprint", fingerprint),
            ]),
        ),
    ]);
    write(&out.join("summary.json"), &scio::json::to_file_text(&summary))?;
    eprint!("scrust cluster: {} cells; {}", r.n_cells, cluster::output::describe(&r));
    eprintln!("scrust cluster: network {:.2} s, resolutions {:.2} s; {threads} threads", r.network_seconds, r.resolutions_seconds);
    Ok(())
}
