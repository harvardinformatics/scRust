//! `scrust doubletfinder`: DoubletFinder as the workflow calls it (docs/PLAN.md B1 and B3).
//!
//! From a 10x counts directory and the cells' cluster labels, every stage is computed by Rust
//! (milestone M6): the sweep subsample and the artificial doublets, drawn from R's
//! random-number stream; SCTransform and RunPCA in each of the six sweep iterations and the final
//! run; neighbours, pANN, the sweep summary and pK, nExp, the calls. With `--from-stages <dir>`,
//! the sweep subsample and the PCA embeddings come from an instrumented R run instead (teacher
//! forcing, milestone M3). The outputs have the reference runs' layout (`cells.tsv`,
//! `sweep_stats.tsv`, `bcmvn.tsv`, `summary.json`, and with `--dump-stages` the stage arrays),
//! so runs compare file to file. `singlets/` holds the counts without the doublets, as a new
//! 10x directory (milestone M6b): the step's result at the level of the counts, which the next
//! tool reads.
//!
//! Milestone M15 (docs/PLAN.md C3h) adds two in-process uses of the clustering chain's library
//! entry point ([`cluster::chain::cluster_counts`]):
//!
//! - **clustering on the fly:** given no `--clusters`, the command clusters its own input first,
//!   as the workflow's upstream script did, and hands the labels with their levels to DoubletFinder
//!   in memory. Only the labels are kept (and the scores, with `--cluster-stats`); the counts read
//!   once serve both. The clustering is written to `clusters/`; at the hand-off the allocator is
//!   settled ([`scio::memory::settle`]), so that DoubletFinder's peak is its own;
//! - **re-clustering the singlets** (`--recluster`): after the calls, the singlets go from memory
//!   to the chain, as `doubletfinder.R` re-clusters them before saving, and `singlets_clusters/` gets
//!   what `scrust cluster --counts <out>/singlets` writes, through the same code.
//!
//! Both directories are nested outputs (docs/output_manifest.md): each is a `scrust cluster`
//! output with its own `summary.json` (the clustering's numbers) and manifest, which the parent's
//! manifest lists by that `summary.json`, with the role `nested_output`.

use crate::cli::{self, Args, Phases};
use crate::cmd_cluster::{self, Record, Settings};
use cluster::chain::{self, Clustered, Options};
use cluster::chain_output::Writer;
use cluster::chain_params::{self, ChainParams};
use doubletfinder::chain::{rng_digest, RustSource};
use doubletfinder::input::{read_counts, read_input};
use doubletfinder::pann::NeighbourDumps;
use doubletfinder::record::{DumpLevel, Recorder};
use doubletfinder::run::{run, Input, Params, PkChoice};
use doubletfinder::singlets::DIR as SINGLETS;
use doubletfinder::stages::{stage_dumps, StageSource};
use scio::manifest::{Entry, Ids, Input as ProvInput, InputFile, Manifest, Provenance};
use scio::Json;
use std::borrow::Cow;
use std::cell::RefCell;
use std::path::{Path, PathBuf};

/// Clustering on the fly's directory (C3h).
pub const CLUSTERS: &str = "clusters";
/// The singlets' re-clustering's directory (M15).
pub const SINGLETS_CLUSTERS: &str = "singlets_clusters";

pub const USAGE: &str = "\
Usage: scrust doubletfinder --counts <10x dir> [--clusters <tsv> [--cluster-levels <txt>]] --out <dir>
                            [--seed 12345] [--pcs 10] [--pn 0.25] [--nexp-rate 0.15]
                            [--pk auto|<value>] [--threads N]
                            [--cluster-preset workflow] [--cluster-params <file.toml>]
                            [--cluster-stats none|silhouette,purity|silhouette|purity]
                            [--recluster [--recluster-preset workflow] [--recluster-params <file.toml>]
                                         [--recluster-export]]
                            [--dump-stages none|light|standard|full] [--from-stages <dir>]

DoubletFinder 2.0.6 as the workflow calls it (docs/PLAN.md B1), every stage computed by Rust:
the sweep subsample and the artificial doublets drawn from R's random-number stream,
SCTransform and RunPCA in each of the six sweep iterations and in the final run, neighbours and
pANN, summarizeSweep, find.pK, nExp with the homotypic adjustment, and the calls. SCTransform's
and RunPCA's products are the OpenBLAS calls R makes, so the results are R's bit for bit under
the OpenBLAS kernel R ran with (OPENBLAS_CORETYPE; see env.sh).

  --counts <dir>          10x Matrix Market directory: matrix.mtx(.gz), barcodes.tsv(.gz),
                          features.tsv(.gz): the RNA counts, stored zeros kept
  --clusters <tsv>        cluster labels: a header line, then barcode and label per cell, in
                          the matrix's column order (the export's seurat_clusters.tsv, or
                          scrust cluster's). Without it, the input is clustered first (below)
  --cluster-levels <txt>  the factor levels in order, one per line (default:
                          seurat_clusters_levels.txt next to --clusters, else the order of
                          first appearance)
  --out <dir>             output directory (must not hold a previous run)
  --seed N                set.seed(N) before paramSweep and before doubletFinder (default
                          12345, the workflow's)
  --pcs N                 PCs = 1:N, RunPCA(npcs = N) (default 10)
  --pn X                  doubletFinder's pN (default 0.25)
  --nexp-rate X           nExp_poi = round(X * cells) (default 0.15)
  --pk auto|X             auto (default): find.pK's choice from the sweep; a value is read
                          with R's parser and skips the sweep
  --threads N             threads (default: SLURM_CPUS_PER_TASK, else 1); the results do not
                          depend on it
  --dump-stages <level>   none (default), light, standard or full: every stage array, RNG
                          checkpoint and scalar in the instrumented R runs' layout, in
                          <out>/stages/; standard adds scale.data per context, full the
                          residual matrices and the S4 order positions
  --from-stages <dir>     teacher forcing (milestone M3): the sweep subsample and the PCA
                          embeddings from an instrumented R run's stages/ (dumps: S4-S7);
                          needs --clusters

Clustering on the fly (without --clusters; docs/PLAN.md C3h): the input is clustered first, in
this process, by the chain scrust cluster runs (SCTransform, RunPCA, FindNeighbors, FindClusters),
as the workflow clustered DoubletFinder's inputs, and its labels and levels feed DoubletFinder's
homotypic adjustment. The counts are read once for both; only the labels are kept.
  --cluster-preset <name>     workflow (default; the workflow's calls, seed 12345) or seurat
                              (Seurat's defaults); applied first
  --cluster-params <file>     a parameter file as scrust cluster --params reads it (sections
                              [sctransform], [pca], [neighbors], [cluster], [stats])
  --cluster-stats <list>      the scores: none (default), silhouette,purity, silhouette or purity
                              (the file's [stats] compute, when given and this is not)
  Writes clusters/: seurat_clusters.tsv and seurat_clusters_levels.txt (the layout --clusters
  reads), the other label files and metadata.tsv (with the scores when asked), and its own
  summary.json (the clustering's numbers and manifest), as scrust cluster --counts <dir> writes
  them. No pca/, neighbors/ or sct/: only the labels are kept. At the hand-off to DoubletFinder,
  glibc's mmap threshold is fixed at 128 KiB and the free pages returned, so that the
  clustering's freed buffers do not raise DoubletFinder's peak (Linux with glibc).

Re-clustering the singlets (--recluster; M15): after the calls, the singlets are clustered in
this process, from memory, as the workflow's doubletfinder.R re-clusters them before saving
(SCTransform, require_min_cells_for_pca(context = \"doubletfinder\"), RunPCA, FindNeighbors,
FindClusters, add_silhouette_to_metadata; RunUMAP changes nothing downstream and is left out).
singlets_clusters/ holds what scrust cluster --counts <out>/singlets writes, byte for byte, apart
from the timings and the command line in its summary.json.
  --recluster-preset <name>   workflow (default) or seurat
  --recluster-params <file>   a parameter file, as scrust cluster --params reads it
  --recluster-export          also the export, as scrust cluster --export writes it (counts/: a
                              copy of singlets/; sct/: SCTransform's normalised data; README.txt)

Outputs: cells.tsv (barcode, pANN and classification of every cell), sweep_stats.tsv (pN, pK,
BCreal), bcmvn.tsv (find.pK's table), summary.json (\"result\": pK and its label as R writes it,
nExp_poi, the homotypic proportion, nExp, the doublets made, k, the column names R writes, the
digest of the final .Random.seed, and for clusters/ and singlets_clusters/ their directory,
cells, clusters and levels, the rest being in their own summary.json; \"run\": parameters,
versions, timings and peak memory per phase; \"files\": the manifest, which lists clusters/ and
singlets_clusters/ by their own summary.json, as nested outputs), and singlets/: the counts without the doublets, as a new
10x directory (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz): every singlet's column exactly
as stored, in input order, and the input's feature lines. The counts are read again from
--counts to write it, so that they need not stay in memory through the final run; with
--recluster the singlets then go to the chain from memory.
";

/// The switches (options without a value).
const SWITCHES: [&str; 2] = ["recluster", "recluster-export"];

/// The clustering's parameters from `--<prefix>-preset` (default: the file's `preset`, else
/// `workflow`) and `--<prefix>-params`, with `stats` (a `--stats` value) applied last. Returns the
/// parameters, the preset applied and the file.
fn chain_params_of(args: &Args, prefix: &str, stats: Option<&str>) -> Result<(ChainParams, String, Option<String>), String> {
    let file = args.get(&format!("{prefix}-params")).map(str::to_string);
    let toml = match &file {
        Some(p) => Some(std::fs::read_to_string(p).map_err(|e| format!("--{prefix}-params {p}: {e}"))?),
        None => None,
    };
    let file_preset = match &toml {
        Some(t) => chain_params::read_toml(t).map_err(|e| format!("--{prefix}-params: {e}"))?.0,
        None => None,
    };
    let preset = args.get(&format!("{prefix}-preset")).map(str::to_string).or(file_preset).unwrap_or_else(|| "workflow".into());
    let flags: Vec<(String, String)> = stats.map(|s| vec![("stats".to_string(), s.to_string())]).unwrap_or_default();
    let p = chain_params::build(Some(&preset), toml.as_deref(), &flags).map_err(|e| format!("--{prefix}-preset/--{prefix}-params: {e}"))?;
    Ok((p, preset, file))
}

/// Whether a parameter file asks for scores (`[stats] compute`).
fn file_sets_stats(args: &Args, prefix: &str) -> Result<bool, String> {
    let Some(p) = args.get(&format!("{prefix}-params")) else { return Ok(false) };
    let t = std::fs::read_to_string(p).map_err(|e| format!("--{prefix}-params {p}: {e}"))?;
    Ok(chain_params::read_toml(&t).map_err(|e| format!("--{prefix}-params: {e}"))?.1.iter().any(|(k, _)| k == "stats"))
}

/// The `scrust cluster` command whose values a nested clustering reproduces (its manifest's
/// `provenance.command`): `--counts <dir> --preset <p> [--params <file>] [--stats <s>] [--export]`.
fn equivalent_command(counts: &Path, preset: &str, params: Option<&str>, stats: Option<&str>, export: bool) -> String {
    let mut a: Vec<String> = vec!["--counts".into(), counts.display().to_string(), "--preset".into(), preset.into()];
    if let Some(p) = params {
        a.extend(["--params".into(), p.into()]);
    }
    if let Some(s) = stats {
        a.extend(["--stats".into(), s.into()]);
    }
    if export {
        a.push("--export".into());
    }
    let mut switches: Vec<&str> = chain_params::PARAM_SWITCHES.to_vec();
    switches.extend(["keep-stages", "export"]);
    scio::manifest::canonical_command("cluster", &a, &cmd_cluster::OUTPUT_ONLY, &switches)
}

/// Marks each stage of an in-process clustering on the command's clock.
struct Marks<'a> {
    phases: &'a RefCell<Phases>,
    prefix: &'static str,
}

impl chain::Observer for Marks<'_> {
    fn phase(&mut self, name: &str) {
        self.phases.borrow_mut().mark(&format!("{}{name}", self.prefix));
    }
}

/// The entry in the parent's manifest for a nested output's own `summary.json`.
fn nested_entry(dir: &str, what: &str) -> Entry {
    Entry::new(&format!("{dir}/summary.json"), scio::manifest::NESTED_OUTPUT, "none", "json", what)
}

/// What the parent's `result` says of a nested clustering: where it is, its cells, clusters and
/// levels. The clustering's numbers are in the nested output's own `summary.json`.
fn nested_result(dir: &str, c: &cluster::chain::Chain) -> Json {
    let levels = &c.clusters.seurat_clusters.levels;
    Json::obj(vec![
        ("dir", Json::Str(dir.into())),
        ("n_cells", Json::Int(c.n_cells as i64)),
        ("n_clusters", Json::Int(levels.len() as i64)),
        ("seurat_clusters_levels", Json::strs(levels)),
    ])
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let args = Args::parse_with_flags(argv, &SWITCHES)?;
    args.only(&[
        "counts", "clusters", "cluster-levels", "from-stages", "out", "threads", "pcs", "pn", "nexp-rate", "pk", "seed", "dump-stages",
        "cluster-preset", "cluster-params", "cluster-stats", "recluster", "recluster-preset", "recluster-params", "recluster-export",
    ])?;
    let phases = RefCell::new(Phases::new());
    let start_time = cli::utc_now();
    let counts = args.path("counts")?;
    let clusters_path = args.get("clusters").map(PathBuf::from);
    let out = args.path("out")?;
    let stages = args.get("from-stages").map(PathBuf::from);
    let threads = cli::threads(&args)?;
    let dump_arg = args.get("dump-stages").unwrap_or("none");
    let dump = DumpLevel::parse(dump_arg).map_err(|e| format!("--dump-stages: {e}"))?;
    let neighbour_dumps = match dump {
        None => NeighbourDumps::default(),
        Some(DumpLevel::Full) => NeighbourDumps { digests: true, positions: true },
        Some(_) => NeighbourDumps { digests: true, positions: false },
    };
    let pk = match args.get("pk").unwrap_or("auto") {
        "auto" => PkChoice::Auto,
        // R reads a pK given in code with its own parser
        v => PkChoice::Fixed(rcompat::as_numeric(v)),
    };
    if let PkChoice::Fixed(v) = pk {
        if !v.is_finite() || v <= 0.0 {
            return Err(format!("--pk: invalid value '{}'", args.get("pk").unwrap()));
        }
    }
    let params = Params {
        pcs: args.number("pcs", 10)?,
        pn: rcompat::as_numeric(args.get("pn").unwrap_or("0.25")),
        nexp_rate: rcompat::as_numeric(args.get("nexp-rate").unwrap_or("0.15")),
        pk,
        seed: args.number("seed", 12345)?,
        threads,
        keep_stages: dump.is_some(),
        neighbour_dumps,
    };
    if params.pcs < 2 {
        return Err("--pcs: at least 2 (RunPCA's embeddings are u %*% diag(d))".into());
    }
    if !(params.pn > 0.0 && params.pn < 1.0) {
        return Err(format!("--pn: expected a value in (0, 1), not {}", params.pn));
    }
    // clustering on the fly (C3h): only without --clusters
    let on_the_fly = clusters_path.is_none();
    if !on_the_fly {
        if let Some(o) = ["cluster-preset", "cluster-params", "cluster-stats"].into_iter().find(|o| args.get(o).is_some()) {
            return Err(format!("--{o}: clustering on the fly is for runs without --clusters"));
        }
    } else {
        if args.get("cluster-levels").is_some() {
            return Err("--cluster-levels goes with --clusters; clustering on the fly makes the levels itself".into());
        }
        if stages.is_some() {
            return Err("--from-stages replays an R run, whose cluster labels it needs: give --clusters".into());
        }
    }
    let recluster = args.flag("recluster");
    if !recluster {
        if let Some(o) = ["recluster-preset", "recluster-params", "recluster-export"].into_iter().find(|o| args.get(o).is_some()) {
            return Err(format!("--{o} goes with --recluster"));
        }
    }
    // the scores of clustering on the fly: none unless asked (C3h: only the labels are kept)
    let hook = if on_the_fly {
        let stats = match args.get("cluster-stats") {
            Some(s) => Some(s.to_string()),
            None if file_sets_stats(&args, "cluster")? => None,
            None => Some("none".to_string()),
        };
        let (p, preset, file) = chain_params_of(&args, "cluster", stats.as_deref())?;
        Some((p, preset, file, stats))
    } else {
        None
    };
    let re = if recluster {
        let (p, preset, file) = chain_params_of(&args, "recluster", None)?;
        Some((p, preset, file, args.flag("recluster-export")))
    } else {
        None
    };
    for f in ["cells.tsv", "summary.json", SINGLETS, CLUSTERS, SINGLETS_CLUSTERS] {
        if out.join(f).exists() {
            return Err(format!("{}: holds a previous run; choose another --out", out.display()));
        }
    }

    // the input: the counts, and the labels given or computed
    let mut hook_result: Option<Clustered> = None;
    let (input, levels_described, matrix) = match &clusters_path {
        Some(cp) => {
            let (input, levels) = read_input(&counts, cp, args.get("cluster-levels").map(Path::new)).map_err(|e| e.to_string())?;
            let m = if stages.is_none() { Some(read_counts(&counts, &input).map_err(|e| e.to_string())?) } else { None };
            (input, levels.describe(), m)
        }
        None => {
            let (hp, ..) = hook.as_ref().expect("clustering on the fly's parameters");
            let c = scio::Counts::read(&counts).map_err(|e| e.to_string())?;
            phases.borrow_mut().mark("load");
            // the chain borrows the counts, which DoubletFinder reads next; it keeps the labels
            // (and the scores when asked) and frees each stage once the next has what it needs
            let mut marks = Marks { phases: &phases, prefix: "clustering: " };
            let clustered = chain::cluster_counts(Cow::Borrowed(&c), hp, threads, Options::default(), &mut marks).map_err(|e| format!("clustering on the fly: {e}"))?;
            let f = clustered.labels();
            let scio::Counts { matrix, barcodes, feature_lines } = c;
            let input = Input { barcodes, clusters: f.label_strings(), cluster_levels: f.levels.clone(), n_genes: feature_lines.len() };
            drop(feature_lines);
            hook_result = Some(clustered);
            (input, format!("{CLUSTERS}/seurat_clusters_levels.txt (clustering on the fly)"), Some(matrix))
        }
    };
    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    // the input's files with their checksums, digested once for every manifest that names them
    let counts_digest = || -> Result<ProvInput, String> {
        ProvInput::digest("counts", &counts, &scio::tenx::input_files(&counts).map_err(|e| e.to_string())?.iter().map(|p| p.as_path()).collect::<Vec<_>>()).map_err(|e| e.to_string())
    };
    let mut counts_prov: Option<ProvInput> = None;
    // clusters/: the labels the run computed, as a nested scrust cluster output, whose own
    // summary.json holds the clustering's numbers; then the result goes, but for what the parent's
    // result says of it (result.clustering) and its number of clusters
    let mut hook_done: Option<(Json, usize)> = None;
    if let (Some(c), Some((hp, preset, file, stats))) = (hook_result.take(), &hook) {
        let dir = out.join(CLUSTERS);
        let s = Settings { params: hp.clone(), threads, keep: false, export: false, dumps: None, purity_lists: false, annoy: false };
        let rec = Record {
            command_line: cli::command_line("doubletfinder", argv),
            command: equivalent_command(&counts, preset, file.as_deref(), stats.as_deref(), false),
            start_time: start_time.clone(),
            inputs: vec![("counts", Json::Str(counts.display().to_string()))],
            provenance_inputs: vec![counts_prov.get_or_insert(counts_digest()?).clone()],
            counts_entries: Vec::new(),
            export_s: 0.0,
            metadata: None,
            stats_note: Vec::new(),
        };
        phases.borrow_mut().mark("clustering: digest");
        std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        let mut w = Writer::new(&dir, c.cells.clone(), threads);
        w.write_stages = false;
        w.write_clustered(&c).map_err(|e| e.to_string())?;
        let mut p = Phases::new();
        p.mark("load");
        let Clustered { chain: ch, .. } = c;
        let d = cmd_cluster::finish_run(&dir, "counts", &s, rec, w, ch, p)?;
        hook_done = Some((nested_result(CLUSTERS, &d.chain), d.chain.clusters.seurat_clusters.levels.len()));
        phases.borrow_mut().mark("clustering: write");
        // the hand-off to DoubletFinder: the clustering's freed memory goes back to the system and
        // glibc's mmap threshold is fixed, so that DoubletFinder's peak is its own (scio::memory)
        scio::memory::settle();
        phases.borrow_mut().mark("clustering: settle");
    }

    let mut rng_end = None;
    let mut recorder = None;
    let output = match &stages {
        Some(dir) => {
            phases.borrow_mut().mark("load");
            let mut source = StageSource::new(dir).map_err(|e| e.to_string())?;
            run(&input, &params, &mut source, &mut |p| phases.borrow_mut().mark(p)).map_err(|e| e.to_string())?
        }
        None => {
            let m = matrix.expect("the counts");
            if !on_the_fly {
                phases.borrow_mut().mark("load");
            }
            let mut source = RustSource::new(m, params.seed, params.pcs, threads).with_marks(|p| phases.borrow_mut().mark(p));
            if let Some(level) = dump {
                source = source.with_recorder(Recorder::new(&out.join("stages"), level));
            }
            let output = run(&input, &params, &mut source, &mut |p| phases.borrow_mut().mark(p)).map_err(|e| e.to_string())?;
            rng_end = Some(source.rng_state());
            recorder = source.take_recorder();
            output
        }
    };
    drop(input);

    output.write_tables(&out).map_err(|e| e.to_string())?;
    match recorder {
        Some(r) => r.finish(&output).map_err(|e| format!("{}: {e}", out.join("stages").display()))?,
        None if dump.is_some() => stage_dumps(&output).write(&out.join("stages")).map_err(|e| e.to_string())?,
        None => {}
    }
    let mut phases = phases.into_inner();
    phases.mark("write");
    // singlets/: the counts read again (they were freed for the final run), the singlets kept in
    // memory, written as a new 10x directory, each file digested as it is written
    let software = format!("scrust {} doubletfinder", env!("CARGO_PKG_VERSION"));
    let full = scio::Counts::read(&counts).map_err(|e| e.to_string())?;
    let singlets = doubletfinder::singlets::singlets(full, &output).map_err(|e| format!("{}: {e}", counts.display()))?;
    let written = doubletfinder::singlets::write_counts(&singlets, &out, &software, threads).map_err(|e| format!("{}: {e}", out.join(SINGLETS).display()))?;
    phases.mark("singlets");
    let result = output.result_json(rng_end.as_deref().map(rng_digest));
    let n_doublets = output.doublet.iter().filter(|&&d| d).count();
    let (n_cells, pk_label, nexp, best_row) = (output.n_cells, output.pk_label.clone(), output.nexp, output.sweep.as_ref().map(|s| s.best_row));

    // singlets_clusters/: the singlets, handed over in memory, through the chain's entry point, written as
    // scrust cluster --counts <out>/singlets writes them (doubletfinder.R's steps after the calls)
    let mut re_done: Option<cmd_cluster::Done> = None;
    if let Some((rp, preset, file, export)) = &re {
        // DoubletFinder's own stages are gone; only its result's tables remain
        drop(output);
        let sdir = out.join(SINGLETS);
        let rdir = out.join(SINGLETS_CLUSTERS);
        let (input, entries, export_s) = if *export {
            // the export's counts/ is singlets/ copied byte for byte, with its checksums
            cmd_cluster::counts_input(&sdir, &rdir, true)?
        } else {
            let files = written.files.iter().map(|f| InputFile { name: f.name.clone(), bytes: f.bytes, sha256: f.sha256.clone() }).collect();
            (ProvInput { option: "counts".into(), path: sdir.display().to_string(), files }, Vec::new(), 0.0)
        };
        let s = Settings { params: rp.clone(), threads, keep: false, export: *export, dumps: None, purity_lists: false, annoy: false };
        let rec = Record {
            command_line: cli::command_line("doubletfinder", argv),
            command: equivalent_command(&sdir, preset, file.as_deref(), None, *export),
            start_time: cli::utc_now(),
            inputs: vec![("counts", Json::Str(sdir.display().to_string()))],
            provenance_inputs: vec![input],
            counts_entries: entries,
            export_s,
            metadata: None,
            stats_note: Vec::new(),
        };
        let mut p = Phases::new();
        p.mark("load");
        let d = cmd_cluster::run_counts(singlets, &rdir, &s, rec, p, Some("doubletfinder")).map_err(|e| format!("re-clustering the singlets: {e}"))?;
        // the re-clustering's own phases, on this command's clock
        let offset = d.started.saturating_duration_since(phases.started()).as_secs_f64();
        for (name, t, rss) in &d.phases {
            phases.marks.push((format!("recluster: {name}"), t + offset, *rss));
        }
        re_done = Some(d);
        phases.mark("recluster");
    } else {
        drop(singlets);
        drop(output);
    }

    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let mode = match &stages {
        Some(_) => "teacher-forced: R's sweep subsample and PCA embeddings (stages S1, S3)",
        None => "rust: every stage (S1-S7) computed by Rust from the counts, the cluster labels and the seed",
    };
    let pk_text = match params.pk {
        PkChoice::Auto => Json::Str("auto".into()),
        PkChoice::Fixed(v) => Json::Num(v),
    };
    let mut parameters = vec![
        ("seed", Json::Int(params.seed as i64)),
        ("pcs", Json::Int(params.pcs as i64)),
        ("pN", Json::Num(params.pn)),
        ("nexp_rate", Json::Num(params.nexp_rate)),
        ("pK", pk_text),
    ];
    if let Some((hp, preset, file, stats)) = &hook {
        parameters.push((
            "clustering",
            Json::obj(vec![
                ("preset", Json::Str(preset.clone())),
                ("params_file", file.clone().map(Json::Str).unwrap_or(Json::Null)),
                ("stats", stats.clone().map(Json::Str).unwrap_or(Json::Null)),
                ("parameters", cmd_cluster::params_json(hp)),
            ]),
        ));
    }
    if let Some((rp, preset, file, export)) = &re {
        parameters.push((
            "recluster",
            Json::obj(vec![
                ("preset", Json::Str(preset.clone())),
                ("params_file", file.clone().map(Json::Str).unwrap_or(Json::Null)),
                ("export", Json::Bool(*export)),
                ("parameters", cmd_cluster::params_json(rp)),
            ]),
        ));
    }
    let parameters = Json::obj(parameters);
    let mut reproduces = vec![doubletfinder::REFERENCE, sctransform::REFERENCE, pca::REFERENCE];
    if hook.is_some() || re.is_some() {
        reproduces.extend([neighbors::REFERENCE, cluster::REFERENCE]);
    }
    let clusters_json = match &clusters_path {
        Some(cp) => Json::Str(cp.display().to_string()),
        None => Json::Str(format!("{CLUSTERS}/seurat_clusters.tsv (clustering on the fly)")),
    };
    let mut run_items = vec![
        ("kind", Json::Str("scrust doubletfinder".into())),
        ("mode", Json::Str(mode.into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        ("reproduces", Json::strs(&reproduces)),
        ("counts", Json::Str(counts.display().to_string())),
        ("clusters", clusters_json),
        ("cluster_levels", Json::Str(levels_described)),
        ("from_stages", stages.as_ref().map(|d| Json::Str(d.display().to_string())).unwrap_or(Json::Null)),
        ("outdir", Json::Str(out.display().to_string())),
        ("parameters", parameters.clone()),
        ("threads", Json::Int(threads as i64)),
        ("dump_stages", Json::Str(dump_arg.into())),
        (
            "singlets",
            Json::obj(vec![
                ("dir", Json::Str(SINGLETS.into())),
                ("features", Json::Int(written.nrow as i64)),
                ("cells", Json::Int(written.ncol as i64)),
                ("entries", Json::Int(written.nnz as i64)),
                ("field", Json::Str(written.field.into())),
            ]),
        ),
    ];
    if hook.is_some() {
        run_items.push(("clustering", Json::obj(vec![("dir", Json::Str(CLUSTERS.into())), ("on_the_fly", Json::Bool(true))])));
    }
    if re.is_some() {
        run_items.push(("recluster", Json::obj(vec![("dir", Json::Str(SINGLETS_CLUSTERS.into())), ("singlets_from", Json::Str("memory".into()))])));
    }
    run_items.extend([
        ("start_time", Json::Str(start_time)),
        ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
        ("elapsed_s", Json::Num(phases.elapsed())),
        ("phases", phases.to_json()),
        ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
        ("fingerprint", fingerprint),
    ]);
    let run_info = Json::obj(run_items);
    // the manifest (docs/PLAN.md C3f, M13): every file, its role, and how the output was made
    let mut inputs = vec![match counts_prov.take() {
        Some(i) => i,
        None => counts_digest()?,
    }];
    if let Some(cp) = &clusters_path {
        inputs.push(ProvInput::digest("clusters", cp.parent().unwrap_or(Path::new("")), &[cp.as_path()]).map_err(|e| e.to_string())?);
        let levels_file = match args.get("cluster-levels") {
            Some(l) => Some(PathBuf::from(l)),
            None => cp.parent().map(|d| d.join("seurat_clusters_levels.txt")).filter(|p| p.is_file()),
        };
        if let Some(l) = &levels_file {
            inputs.push(ProvInput::digest("cluster-levels", l.parent().unwrap_or(Path::new("")), &[l.as_path()]).map_err(|e| e.to_string())?);
        }
    }
    for (option, file) in [("cluster-params", args.get("cluster-params")), ("recluster-params", args.get("recluster-params"))] {
        if let Some(f) = file {
            let p = Path::new(f);
            inputs.push(ProvInput::digest(option, p.parent().unwrap_or(Path::new("")), &[p]).map_err(|e| e.to_string())?);
        }
    }
    if let Some(d) = &stages {
        inputs.push(ProvInput { option: "from-stages".into(), path: d.display().to_string(), files: Vec::new() });
    }
    let mut seeds = vec![
        ("set.seed before paramSweep and doubletFinder", Json::Int(params.seed as i64)),
        ("SCTransform seed.use", Json::Int(1448145)),
        ("RunPCA seed.use", Json::Int(42)),
    ];
    if let Some((hp, ..)) = &hook {
        seeds.push(("clustering on the fly", cmd_cluster::seeds_json(hp)));
    }
    if let Some((rp, ..)) = &re {
        seeds.push(("re-clustering the singlets", cmd_cluster::seeds_json(rp)));
    }
    let mut m = Manifest::new(
        "doubletfinder",
        Provenance {
            scrust_version: env!("CARGO_PKG_VERSION").into(),
            command_line: crate::cli::command_line("doubletfinder", argv),
            command: scio::manifest::canonical_command("doubletfinder", argv, &["out", "threads", "dump-stages"], &SWITCHES),
            inputs,
            reproduces: reproduces.iter().map(|s| s.to_string()).collect(),
            parameters,
            seeds: Json::obj(seeds),
        },
    );
    m.push(Manifest::summary_entry());
    let cells = |f: &str| Ids::named("cells", f, "barcode");
    m.push(Entry::new("cells.tsv", "doublet_calls", "none", "tsv", "every input cell: barcode, pANN, classification (Singlet or Doublet), as doubletFinder computes them").header(true).rows(cells("cells.tsv")));
    for (f, role, what) in [("sweep_stats.tsv", "doublet_sweep_stats", "summarizeSweep's table: pN, pK, BCreal"), ("bcmvn.tsv", "doublet_bcmvn", "find.pK's table: pK and its BCmetric")] {
        if out.join(f).is_file() {
            m.push(Entry::new(f, role, "none", "tsv", what).header(true));
        }
    }
    let (mf, bf, ff) = (format!("{SINGLETS}/matrix.mtx.gz"), format!("{SINGLETS}/barcodes.tsv.gz"), format!("{SINGLETS}/features.tsv.gz"));
    let digest = |name: &str, e: Entry| -> Entry {
        let mut e = e;
        if let Some(f) = written.files.iter().find(|f| f.name == name) {
            e.bytes = Some(f.bytes);
            e.sha256 = Some(f.sha256.clone());
        }
        e
    };
    m.push(digest("matrix.mtx.gz", Entry::new(&mf, "raw_counts", "raw_counts", "matrix_market", "the singlets' counts: the input counts without the doublets, every kept column exactly as stored (a new 10x directory, Cell Ranger's layout)").shape(&[written.nrow, written.ncol]).nnz(written.nnz).rows(Ids::numbered("features", &ff, 2)).columns(Ids::numbered("cells", &bf, 1))));
    m.push(digest("barcodes.tsv.gz", Entry::new(&bf, "raw_counts_cells", "none", "tsv", "the singlets' barcodes, one per line, in input order").header(false)));
    m.push(digest("features.tsv.gz", Entry::new(&ff, "raw_counts_features", "none", "tsv", "the input's feature lines: id, name, type").header(false)));
    m.entries.extend(scio::manifest::stage_entries(&out, "stages").map_err(|e| e.to_string())?);
    if hook.is_some() {
        m.push(nested_entry(CLUSTERS, "clustering on the fly: the input's clustering, whose labels and levels DoubletFinder used (seurat_clusters.tsv, seurat_clusters_levels.txt), as scrust cluster writes its labels and metadata.tsv; a nested scrust output, whose own manifest (this file) lists its files"));
    }
    if re.is_some() {
        m.push(nested_entry(SINGLETS_CLUSTERS, "the singlets re-clustered, as doubletfinder.R does before saving: what scrust cluster --counts singlets/ writes; a nested scrust output, whose own manifest (this file) lists its files"));
    }
    m.fill_shapes(&out).map_err(|e| e.to_string())?;
    let files = crate::cli::finish_manifest(&out, &mut m, threads)?;
    let mut result = result;
    if let (Json::Object(items), Some((r, _))) = (&mut result, &hook_done) {
        items.push(("clustering".into(), r.clone()));
    }
    if let (Json::Object(items), Some(d)) = (&mut result, &re_done) {
        items.push(("recluster".into(), nested_result(SINGLETS_CLUSTERS, &d.chain)));
    }
    let summary = Json::obj(vec![("result", result), ("run", run_info), ("files", files)]);
    std::fs::write(out.join("summary.json"), scio::json::to_file_text(&summary)).map_err(|e| e.to_string())?;
    eprintln!(
        "scrust doubletfinder: pK {} (BCmetric row {}), nExp {}, {} doublets of {} cells; {}/: {} singlets, {} entries;{}{} {:.1} s, {} threads, peak RSS {:.2} GiB",
        pk_label,
        best_row.map_or("-".into(), |r| (r + 1).to_string()),
        nexp,
        n_doublets,
        n_cells,
        SINGLETS,
        written.ncol,
        written.nnz,
        hook_done.as_ref().map(|(_, n)| format!(" {CLUSTERS}/: {n} clusters (on the fly);")).unwrap_or_default(),
        re_done.as_ref().map(|d| format!(" {SINGLETS_CLUSTERS}/: {} clusters;", d.chain.clusters.seurat_clusters.levels.len())).unwrap_or_default(),
        phases.elapsed(),
        threads,
        cli::peak_rss_kib().unwrap_or(0) as f64 / 1048576.0
    );
    Ok(())
}
