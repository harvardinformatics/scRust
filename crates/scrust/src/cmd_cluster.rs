//! `scrust cluster`: the workflow's clustering chain (docs/PLAN.md C1, C3; milestones M10 and
//! M11): SCTransform, RunPCA, FindNeighbors, FindClusters and the silhouette widths, from the
//! counts in one process, or from any stage's output (`scrust sctransform`, `scrust pca`,
//! `scrust neighbors`).

use crate::cli::{self, Args, Phases};
use cluster::chain::{self, Options};
use cluster::chain_output::Writer;
use cluster::chain_params::{self, ChainParams, PARAM_FLAGS, PARAM_SWITCHES};
use cluster::metadata::RnaMetadata;
use scio::stages::DumpLevel;
use scio::{tenx, Json, Table};
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: scrust cluster (--counts <10x dir> | --sct <dir> | --pca <dir> | --snn <dir> [--pca <dir>])
                      --out <dir> [--preset workflow] [--params <file.toml>] [parameters]
                      [--metadata <tsv>] [--stats silhouette,purity|silhouette|none]
                      [--keep-stages] [--dump-stages none|light|standard|full] [--annoy yes|no]
                      [--threads N]

The workflow's clustering chain as Seurat 5.5.1 runs it, bit for bit (under the reference's
OpenBLAS kernel, for the stages that use BLAS):

  SCTransform(obj, seed.use)              sctransform 0.4.3 v2, glmGamPoi 1.22.0
  require_min_cells_for_pca(obj)          the workflow's guard (on with --preset workflow)
  RunPCA(obj, npcs, seed.use)             irlba 2.4.1
  FindNeighbors(obj, dims, k.param, ...)  RcppAnnoy 0.0.23's Annoy, the SNN graph
  FindClusters(obj, resolution, ...)      the modularity optimizer, GroupSingletons
  add_silhouette_to_metadata(obj)         bluster 1.20.0's approxSilhouette (purity: M12)

RunUMAP, which the workflow runs after RunPCA, changes nothing downstream and is left out.

Starting points:
  --counts <dir>         a 10x Matrix Market directory: the whole chain in one process
  --sct <dir>            `scrust sctransform`'s output: RunPCA onwards
  --pca <dir>            `scrust pca`'s output (embeddings.tsv), or an R run's pca/: FindNeighbors
                         onwards
  --snn <dir>            `scrust neighbors`' output (snn.mtx.gz with knn_index.tsv), or an R run's
                         graphs/: FindClusters onwards; with --pca, the scores too
  Each stage's directory carries the counts' nCount_RNA and nFeature_RNA (cell_metadata.tsv)
  forward, so every route writes the same files.

Parameters (Seurat 5.5.1's defaults; --preset workflow sets the workflow's calls: seed 12345 in
SCTransform, RunPCA and FindClusters, dims 1:30, and the PCA guard). Applied in this order:
the defaults, the preset, the --params file (a TOML file with the sections [sctransform],
[pca], [neighbors], [cluster] and [stats], and optionally preset = \"workflow\"), the flags.
Numbers are read as R reads them; a fraction such as 1/15 is divided in double.
  SCTransform:   --sct-seed 1448145  --ncells 5000  --variable-features 3000
                 --clip-range lo,hi (default +-sqrt(N/30))
  RunPCA:        --npcs 50  --pca-seed 42  --weight-by-var true  --pca-guard none|true|N
  FindNeighbors: --dims 1:10  --k 20  --trees 50  --metric euclidean|cosine|manhattan
                 --prune 1/15  --search-k -1
  FindClusters:  --resolution 0.8[,1.2,...]  --algorithm 1|2|3  --starts 10  --iterations 10
                 --seed 0  --modularity 1|2  --no-group-singletons  --graph-name SCT_snn
  Scores:        --stats silhouette,purity (default) | silhouette | purity | none
                 --purity-k 50  --purity-weighted true
                 The scores use the neighbour search's dims. Purity comes with milestone M12:
                 asked for, the command says it is not computed yet.
  Other code paths are refused: --vars-to-regress, --vst-flavor other than v2, --sct-method
  other than glmGamPoi_offset, --do-scale true, --nn-method rann|hnsw, --l2-norm true,
  --metric hamming, --algorithm 4 (Leiden).

Outputs (--out, which must not hold a previous run):
  seurat_clusters.tsv, seurat_clusters_levels.txt   the labels, as scrust doubletfinder reads them
  clusters.tsv, <column>_levels.txt                 one column per resolution (SCT_snn_res.0.8)
  cluster_ids.txt, nclusters.txt                    as the workflow's write_cluster_metadata
  optimizer_ids.tsv                                 each resolution's optimizer result
  metadata.tsv      barcode, nCount_RNA, nFeature_RNA, each res. column, seurat_clusters,
                    silhouette_width, silhouette_other (the closest other cluster), then the
                    columns of --metadata (a TSV with a barcode column), by barcode
  pca/              embeddings.tsv, loadings.tsv, stdev.tsv, random_seed.txt (scrust pca's)
  neighbors/        knn_index.tsv, knn_distance.tsv, nn.mtx.gz, snn.mtx.gz (scrust neighbors'),
                    and annoy.ann with --annoy yes
  sct/              with --keep-stages, from the counts: scrust sctransform's layout
  stages/           with --keep-stages, every random start's partition; with --dump-stages,
                    the stage dumps of the instrumented R runs, for scrust validate cluster
  summary.json      result: FindClusters' (the network, every start, the best start, what
                    FindClusters printed, the singletons, the levels), each stage's, R's digests
                    of .Random.seed after SCTransform and RunPCA and of scale.data; run:
                    parameters, timings and peak memory per stage, the environment fingerprint
  --export (the corrected counts and the normalised data) comes with milestone M13.
";

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// The parameters from `--preset`, `--params` and the flags.
pub fn params_from_args(args: &Args) -> Result<ChainParams, String> {
    let toml = match args.get("params") {
        Some(p) => Some(std::fs::read_to_string(p).map_err(|e| format!("--params {p}: {e}"))?),
        None => None,
    };
    let mut flags: Vec<(String, String)> = Vec::new();
    for name in PARAM_FLAGS {
        if let Some(v) = args.get(name) {
            flags.push((name.to_string(), v.to_string()));
        }
    }
    for name in PARAM_SWITCHES {
        if args.flag(name) {
            flags.push((name.to_string(), String::new()));
        }
    }
    chain_params::build(args.get("preset"), toml.as_deref(), &flags)
}

/// The parameters as `summary.json` records them (`run.parameters`).
pub fn params_json(p: &ChainParams) -> Json {
    let mut sct = sctransform::output::params_json(&p.sct);
    if let Json::Object(m) = &mut sct {
        m.push(("clip_range".into(), p.sct.clip_range.map(|(a, b)| Json::nums(&[a, b])).unwrap_or(Json::Null)));
    }
    let mut pca = pca::output::params_json(&p.pca);
    if let Json::Object(m) = &mut pca {
        m.push(("guard_npcs".into(), p.pca_guard.map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)));
    }
    Json::obj(vec![
        ("preset", p.preset.clone().map(Json::Str).unwrap_or(Json::Null)),
        ("sctransform", sct),
        ("pca", pca),
        ("neighbors", neighbors::output::params_json(&p.neighbors)),
        ("cluster", cluster::output::params_json(&p.cluster)),
        (
            "stats",
            Json::obj(vec![
                ("requested", Json::strs(&p.stats.names())),
                ("purity_k", Json::Int(p.stats.purity_k as i64)),
                ("purity_weighted", Json::Bool(p.stats.purity_weighted)),
            ]),
        ),
    ])
}

/// The barcodes of a stage directory's cells.
fn table_barcodes(path: &Path) -> Result<Vec<String>, String> {
    let t = Table::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(t.rows.iter().map(|r| r[0].clone()).collect())
}

/// The counts' columns a stage directory carries, checked against its cells.
fn carried_rna(dirs: &[&Path], barcodes: &[String]) -> Result<Option<RnaMetadata>, String> {
    for d in dirs {
        let p = d.join(cluster::metadata::FILE);
        if p.is_file() {
            let (b, m) = RnaMetadata::read(&p)?;
            if b != barcodes {
                return Err(format!("{}: its cells are not the run's", p.display()));
            }
            return Ok(Some(m));
        }
    }
    Ok(None)
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let mut switches: Vec<&str> = PARAM_SWITCHES.to_vec();
    switches.extend(["keep-stages"]);
    let args = Args::parse_with_flags(argv, &switches)?;
    if args.get("export").is_some() || args.flag("export") {
        return Err("--export: the export of the corrected counts and the normalised data comes with milestone M13".into());
    }
    let mut known: Vec<&str> = vec!["counts", "sct", "pca", "snn", "out", "preset", "params", "metadata", "keep-stages", "dump-stages", "annoy", "threads"];
    known.extend(PARAM_FLAGS);
    known.extend(PARAM_SWITCHES);
    args.only(&known)?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let params = params_from_args(&args)?;
    let keep = args.flag("keep-stages");
    let dumps = DumpLevel::parse(args.get("dump-stages").unwrap_or("none"))?;
    let annoy = match args.get("annoy").unwrap_or("no") {
        "yes" => true,
        "no" => false,
        v => return Err(format!("--annoy: expected yes or no, not '{v}'")),
    };
    let start = match (args.get("counts"), args.get("sct"), args.get("pca"), args.get("snn")) {
        (Some(_), None, None, None) => "counts",
        (None, Some(_), None, None) => "sct",
        (None, None, Some(_), None) => "pca",
        (None, None, _, Some(_)) => "snn",
        _ => return Err("give one starting point: --counts, --sct, --pca, or --snn (with --pca for the scores)".into()),
    };
    if dumps.is_some() && start != "counts" {
        return Err("--dump-stages: the stage dumps describe the whole chain from the counts (--counts)".into());
    }
    if out.join("summary.json").exists() || out.join("clusters.tsv").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }
    let opts = Options { keep_partitions: keep, keep_residuals: dumps == Some(DumpLevel::Full) };
    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;

    let mut inputs: Vec<(&str, Json)> = Vec::new();
    let mut stats_note: Vec<String> = Vec::new();
    let result = match start {
        "counts" => {
            let dir = args.path("counts")?;
            let barcodes = tenx::read_barcodes(&dir).map_err(|e| e.to_string())?;
            let features = tenx::read_features(&dir).map_err(|e| e.to_string())?;
            let counts = scio::sparse::read_mtx(&dir).map_err(|e| e.to_string())?;
            if counts.nrow != features.len() || counts.ncol != barcodes.len() {
                return Err(format!("{}: the matrix is {} x {}, with {} features and {} barcodes", dir.display(), counts.nrow, counts.ncol, features.len(), barcodes.len()));
            }
            inputs.push(("counts", Json::Str(dir.display().to_string())));
            let mut w = Writer::new(&out, barcodes, threads);
            w.features = Some(features.iter().map(|f| (f.id.clone(), f.name.clone())).collect());
            w.rna = Some(RnaMetadata::of(&counts));
            w.keep_stages = keep;
            w.annoy = annoy;
            w.neighbor_params = Some(params.neighbors.clone());
            if let Some(level) = dumps {
                w = w.with_dumps(level);
            }
            phases.mark("load");
            let r = chain::from_counts(counts, &params, threads, opts, &mut w);
            (r, w)
        }
        "sct" => {
            let dir = args.path("sct")?;
            let barcodes = table_barcodes(&dir.join("cell_attr.tsv"))?;
            let vf = Table::read(&dir.join("variable_features.tsv")).map_err(|e| format!("{}: {e}", dir.join("variable_features.tsv").display()))?;
            let (rows, nf, nc) = scio::npy::read_f64_row_major(&dir.join("scale_data.npy")).map_err(|e| e.to_string())?;
            if nf != vf.rows.len() || nc != barcodes.len() {
                return Err(format!("{}: scale_data.npy is {nf} x {nc}, with {} variable features and {} cells", dir.display(), vf.rows.len(), barcodes.len()));
            }
            inputs.push(("sct", Json::Str(dir.display().to_string())));
            let rna = carried_rna(&[&dir], &barcodes)?;
            let mut w = Writer::new(&out, barcodes, threads);
            // the PCA's input rows are the variable features, in order
            w.features = Some(vf.rows.iter().map(|r| (r[0].clone(), r.get(1).cloned().unwrap_or_default())).collect());
            w.rna = rna;
            w.keep_stages = keep;
            w.annoy = annoy;
            w.neighbor_params = Some(params.neighbors.clone());
            w.sct_result = std::fs::read_to_string(dir.join("summary.json")).ok().and_then(|t| Json::parse(&t).ok()).and_then(|j| j.get("result").cloned());
            if keep {
                // scrust sctransform's files, as the one-process run writes them
                for f in ["cell_attr.tsv", "gene_attr.tsv", "model_pars.tsv", "model_pars_fit.tsv", "variable_features.tsv", "scale_data.npy"] {
                    std::fs::create_dir_all(out.join("sct")).map_err(|e| e.to_string())?;
                    std::fs::copy(dir.join(f), out.join("sct").join(f)).map_err(|e| format!("{}: {e}", dir.join(f).display()))?;
                    w.files.push(format!("sct/{f}"));
                }
            }
            phases.mark("load");
            let r = chain::from_scale_data(rows, (0..nf).collect(), nc, &params, threads, opts, &mut w);
            (r, w)
        }
        "pca" => {
            let dir = args.path("pca")?;
            let emb = crate::cmd_neighbors::read_embeddings(&dir)?;
            inputs.push(("pca", Json::Str(dir.display().to_string())));
            let rna = carried_rna(&[&dir], &emb.barcodes)?;
            let mut w = Writer::new(&out, emb.barcodes.clone(), threads);
            w.rna = rna;
            w.keep_stages = keep;
            w.annoy = annoy;
            w.neighbor_params = Some(params.neighbors.clone());
            w.pca_result = w.copy_stage(&dir, "pca", &pca::output::FILES).map_err(|e| e.to_string())?;
            phases.mark("load");
            let r = chain::from_embeddings(emb.data, emb.n, emb.ncol, &params, threads, opts, &mut w);
            (r, w)
        }
        _ => {
            let dir = args.path("snn")?;
            let input = cluster::input::read_snn_dir(&dir)?;
            inputs.push(("snn", Json::Str(input.source.display().to_string())));
            inputs.push(("barcodes", Json::Str(input.barcodes_source.display().to_string())));
            let emb = match args.get("pca") {
                Some(p) => {
                    let e = crate::cmd_neighbors::read_embeddings(Path::new(p))?;
                    if e.barcodes != input.barcodes {
                        return Err(format!("--pca {p}: its cells are not the SNN graph's"));
                    }
                    inputs.push(("pca", Json::Str(p.to_string())));
                    Some((PathBuf::from(p), e))
                }
                None => None,
            };
            let gdir = input.source.parent().map(Path::to_path_buf).unwrap_or_default();
            let mut dirs: Vec<&Path> = vec![gdir.as_path()];
            if let Some((p, _)) = &emb {
                dirs.push(p.as_path());
            }
            let rna = carried_rna(&dirs, &input.barcodes)?;
            let mut w = Writer::new(&out, input.barcodes.clone(), threads);
            w.rna = rna;
            w.keep_stages = keep;
            w.annoy = annoy;
            w.neighbor_params = Some(params.neighbors.clone());
            w.neighbors_result = w.copy_stage(&gdir, "neighbors", &["knn_index.tsv", "knn_distance.tsv", "nn.mtx.gz", "snn.mtx.gz", "annoy.ann"]).map_err(|e| e.to_string())?;
            if let Some((p, _)) = &emb {
                w.pca_result = w.copy_stage(p, "pca", &pca::output::FILES).map_err(|e| e.to_string())?;
            }
            if params.stats.silhouette && emb.is_none() {
                let m = "silhouette asked for, but not computed: the scores need the PCA (give --pca with --snn)".to_string();
                eprintln!("scrust cluster: {m}");
                stats_note.push(m);
            }
            phases.mark("load");
            let r = chain::from_snn(&input.snn, emb.as_ref().map(|(_, e)| (e.data.as_slice(), e.ncol)), &params, threads, opts, &mut w);
            (r, w)
        }
    };
    let (r, mut w) = result;
    let r = r.map_err(|e| e.to_string())?;
    for (name, t, rss) in &w.phases {
        phases.marks.push((name.clone(), *t + phases.marks.first().map_or(0.0, |m| m.1), *rss));
    }
    // what FindClusters prints, as Seurat prints it
    for x in &r.clusters.resolutions {
        for line in x.printed() {
            eprintln!("{line}");
        }
        if let Some(m) = x.message() {
            eprintln!("{m}");
        }
    }
    if let Some(m) = args.get("metadata") {
        let mut taken: Vec<&str> = vec!["barcode", "nCount_RNA", "nFeature_RNA", "seurat_clusters", "silhouette_width", "silhouette_other", "neighborhood_purity", "purity_maximum", "nCount_SCT", "nFeature_SCT"];
        let cols: Vec<String> = r.clusters.resolutions.iter().map(|x| x.column.name.clone()).collect();
        taken.extend(cols.iter().map(String::as_str));
        w.extra = Some(cluster::metadata::extra_columns(Path::new(m), &w.barcodes, &taken)?);
        inputs.push(("metadata", Json::Str(m.to_string())));
    }
    w.finish_labels(&r.clusters, r.silhouette.as_ref()).map_err(|e| e.to_string())?;
    w.finish_dumps().map_err(|e| e.to_string())?;
    phases.mark("write");

    let mut computed: Vec<&str> = Vec::new();
    if r.silhouette.is_some() {
        computed.push("silhouette");
    }
    if params.stats.purity {
        let m = "purity asked for, but not computed: neighborhood_purity comes with milestone M12".to_string();
        eprintln!("scrust cluster: {m}");
        stats_note.push(m);
    }
    if params.stats.names().is_empty() {
        stats_note.push("none asked for".into());
    }
    let version = env!("CARGO_PKG_VERSION");
    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let t = &r.timings;
    let summary = Json::obj(vec![
        ("tool", Json::Str("scrust cluster".into())),
        ("reference", Json::Str(format!("{}; {}; {}; {}; bluster 1.20.0 approxSilhouette", sctransform::REFERENCE, pca::REFERENCE, neighbors::REFERENCE, cluster::REFERENCE))),
        ("scrust_version", Json::Str(version.into())),
        ("result", w.result_json(&r.clusters)),
        (
            "run",
            Json::obj(vec![
                ("parameters", params_json(&params)),
                ("start", Json::Str(start.into())),
                ("inputs", Json::obj(inputs)),
                ("outdir", Json::Str(out.display().to_string())),
                ("files", Json::strs(&w.files)),
                (
                    "stats",
                    Json::obj(vec![("requested", Json::strs(&params.stats.names())), ("computed", Json::strs(&computed)), ("note", Json::Str(stats_note.join("; ")))]),
                ),
                ("keep_stages", Json::Bool(keep)),
                ("dump_stages", Json::Str(dumps.map_or("none", |d| d.name()).into())),
                ("threads", Json::Int(threads as i64)),
                ("start_time", Json::Str(start_time)),
                ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
                ("sctransform_s", Json::Num(t.sctransform)),
                ("pca_s", Json::Num(t.pca)),
                ("neighbors_s", Json::Num(t.neighbors)),
                ("cluster_s", Json::Num(t.cluster)),
                ("silhouette_s", Json::Num(t.silhouette)),
                ("network_s", Json::Num(r.clusters.network_seconds)),
                ("resolutions_s", Json::Num(r.clusters.resolutions_seconds)),
                ("optimizer_s", Json::nums(&r.clusters.resolutions.iter().map(|x| x.optimizer.seconds).collect::<Vec<_>>())),
                ("groupsingletons_s", Json::nums(&r.clusters.resolutions.iter().map(|x| x.group_seconds).collect::<Vec<_>>())),
                ("printed_elapsed", Json::strs(&r.clusters.resolutions.iter().map(|x| x.printed().pop().unwrap_or_default()).collect::<Vec<_>>())),
                ("elapsed_s", Json::Num(phases.elapsed())),
                ("phases", phases.to_json()),
                ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
                ("fingerprint", fingerprint),
            ]),
        ),
    ]);
    write(&out.join("summary.json"), &scio::json::to_file_text(&summary))?;
    eprint!("scrust cluster: {} cells; {}", r.n_cells, cluster::output::describe(&r.clusters));
    eprintln!(
        "scrust cluster: SCTransform {:.1} s, RunPCA {:.1} s, FindNeighbors {:.1} s, FindClusters {:.1} s, silhouette {:.2} s; {:.1} s in all, {threads} threads",
        t.sctransform,
        t.pca,
        t.neighbors,
        t.cluster,
        t.silhouette,
        phases.elapsed()
    );
    Ok(())
}
