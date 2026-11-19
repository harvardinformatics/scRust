//! `scrust cluster`: the workflow's clustering chain (docs/PLAN.md C1, C3; milestones M10, M11
//! and M12): SCTransform, RunPCA, FindNeighbors, FindClusters, the silhouette widths and the
//! neighbourhood purities, from the counts in one process, or from any stage's output
//! (`scrust sctransform`, `scrust pca`, `scrust neighbors`). With `--export` (C3f, M13), also
//! SCTransform's corrected counts and the export: the input counts, the normalised data and the
//! per-cell metadata, as plain files that cannot be taken for raw counts.

use crate::cli::{self, Args, Phases};
use cluster::chain::{self, Chain, Options};
use cluster::chain_output::Writer;
use cluster::chain_params::{self, ChainParams, PARAM_FLAGS, PARAM_SWITCHES};
use cluster::covariates::{CellCovariates, MtGenes, FROM_COUNTS};
use cluster::metadata::RnaMetadata;
use scio::manifest::{Entry, Ids, Input, InputFile, Manifest, Provenance};
use scio::stages::DumpLevel;
use scio::{tenx, Json, Table};
use std::borrow::Cow;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: scrust cluster (--counts <10x dir> | --sct <dir> | --pca <dir> | --snn <dir> [--pca <dir>])
                      --out <dir> [--preset workflow] [--params <file.toml>] [parameters]
                      [--metadata <tsv>] [--stats silhouette,purity|silhouette|none] [--export]
                      [--keep-stages] [--dump-stages none|light|standard|full]
                      [--purity-lists auto|yes|no] [--annoy yes|no] [--threads N]

The workflow's clustering chain as Seurat 5.5.1 runs it, bit for bit (under the reference's
OpenBLAS kernel, for the stages that use BLAS):

  SCTransform(obj, seed.use)              sctransform 0.4.3 v2, glmGamPoi 1.22.0
  require_min_cells_for_pca(obj)          the workflow's guard (on with --preset workflow)
  RunPCA(obj, npcs, seed.use)             irlba 2.4.1
  FindNeighbors(obj, dims, k.param, ...)  RcppAnnoy 0.0.23's Annoy, the SNN graph
  FindClusters(obj, resolution, ...)      the modularity optimizer, GroupSingletons
  add_silhouette_to_metadata(obj)         bluster 1.20.0's approxSilhouette and neighborPurity
                                          (its search exhaustive, with knncolle's arithmetic, and
                                          a certificate that BiocNeighbors' KMKNN finds the same)

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
                 --vars-to-regress <name>[,<name>...]  vars.to.regress (default: none): the
                     per-cell covariates SCTransform's second regression removes from scale.data
                     (ScaleData's RegressOutMatrix: lm(GENE ~ vars)'s QR, qr.resid per variable
                     feature, then the centring): percent.mt, nCount_RNA, nFeature_RNA, or any
                     numeric column of --metadata (read as R reads numbers). An unknown or
                     repeated name, a name R cannot use in a formula, or a missing, non-numeric or
                     non-finite value is an error
                 --percent-mt  compute percent.mt without regressing it, as the workflow's
                     objects hold it (metadata.tsv; implied by --mt-genes and by regressing it)
                 --mt-genes <file>  percent.mt's mitochondrial features by id: one id per line
                     (column 1 of the counts' features.tsv; blank lines and lines starting with #
                     are skipped). Ids that no feature has are reported and skipped; a list that
                     matches none is an error. Without it: the workflow's \"(?i)^mt-\" on Seurat's
                     row names (features.tsv's column 2 made unique, _ and | replaced by -).
                     percent.mt is PercentageFeatureSet(): their counts over nCount_RNA, times
                     100. Regressing percent.mt when no feature matches the pattern is an error
                     that names --mt-genes; otherwise the column is written as computed (zeros
                     included, with a warning)
  RunPCA:        --npcs 50  --pca-seed 42  --weight-by-var true  --pca-guard none|true|N
  FindNeighbors: --dims 1:10  --k 20  --trees 50  --metric euclidean|cosine|manhattan
                 --prune 1/15  --search-k -1
  FindClusters:  --resolution 0.8[,1.2,...]  --algorithm 1|2|3  --starts 10  --iterations 10
                 --seed 0  --modularity 1|2  --no-group-singletons  --graph-name SCT_snn
  Scores:        --stats silhouette,purity (default) | silhouette | purity | none
                 --purity-k 50  --purity-weighted true
                 The scores use the neighbour search's dims (1:30 with the preset).
  The TOML file's [sctransform] takes vars_to_regress = [\"percent.mt\", ...], percent_mt = true
  and mt_genes = \"<file>\" (a relative path is read from the TOML file's directory).
  Other code paths are refused: --vst-flavor other than v2, --sct-method other than
  glmGamPoi_offset, --do-scale true, --nn-method rann|hnsw, --l2-norm true, --metric hamming,
  --algorithm 4 (Leiden). --vars-to-regress, --mt-genes and --percent-mt act where SCTransform
  runs: --counts here, or scrust sctransform (whose cell_metadata.tsv carries percent.mt to the
  staged starts).

Outputs (--out, which must not hold a previous run):
  seurat_clusters.tsv, seurat_clusters_levels.txt   the labels, as scrust doubletfinder reads them
  clusters.tsv, <column>_levels.txt                 one column per resolution (SCT_snn_res.0.8)
  cluster_ids.txt, nclusters.txt                    as the workflow's write_cluster_metadata
  optimizer_ids.tsv                                 each resolution's optimizer result
  metadata.tsv      barcode, nCount_RNA, nFeature_RNA, percent.mt (when asked for), each res.
                    column, seurat_clusters, silhouette_width, silhouette_other (the closest
                    other cluster), neighborhood_purity, purity_maximum (the cluster with the
                    largest share), then the columns of --metadata (a TSV with a barcode column),
                    by barcode
  pca/              embeddings.tsv, loadings.tsv, stdev.tsv, random_seed.txt (scrust pca's)
  neighbors/        knn_index.tsv, knn_distance.tsv, nn.mtx.gz, snn.mtx.gz (scrust neighbors'),
                    and annoy.ann with --annoy yes
  sct/              with --keep-stages, from the counts: scrust sctransform's layout
  stages/           with --keep-stages, every random start's partition; with --dump-stages,
                    the stage dumps of the instrumented R runs, for scrust validate cluster
                    (--purity-lists: purity's neighbour lists among them; auto: at level full)
  summary.json      result: FindClusters' (the network, every start, the best start, what
                    FindClusters printed, the singletons, the levels), each stage's (SCTransform's
                    with the regression and percent.mt's definition, features and ids not found,
                    when asked for), R's digests of .Random.seed after SCTransform and RunPCA
                    and of scale.data, the scores
                    (purity: the radius and the certificate, every flagged pair); run:
                    parameters, timings and peak memory per stage, the environment fingerprint;
                    files: the manifest, every file's role, scale, shape, identifiers, size and
                    SHA-256, and how the output was made (docs/output_manifest.md)

--export (from --counts only) writes the export (docs/PLAN.md C3f): SCTransform's normalised data
with the input counts and the per-cell metadata, as plain files that no 10x reader takes for raw
counts. In addition to the files above:
  counts/           the input counts (matrix, barcodes, features), copied byte for byte with their
                    checksums: the export's only 10x directory, and the only raw counts
  sct/              scrust sctransform's layout (the model, scale_data.npy), and
                      corrected_counts.mtx.gz        Seurat's SCT counts layer: sctransform's
                                                     correct() counts (integers), NOT raw counts
                      log1p_corrected_counts.mtx.gz  Seurat's SCT data layer: their log1p
                      sct_features.tsv.gz, sct_cells.tsv.gz   their rows and columns
                    each matrix saying in its comment lines what it is, scrust's version and
                    command, and the input's checksums
  metadata.tsv      gains nCount_SCT and nFeature_SCT (the corrected counts' column sums and
                    nonzero counts), after nFeature_RNA
  README.txt        what each file holds and how to load it
The other routes refuse --export: it copies the input counts and corrects every model gene from
vst's residuals, which only a run from the counts has.
";

/// Options that change no value of any file the normalised matrices' comment lines describe:
/// left out of the command those lines give, so that the files do not depend on them.
pub const OUTPUT_ONLY: [&str; 7] = ["out", "threads", "keep-stages", "dump-stages", "purity-lists", "annoy", "metadata"];

/// The files of a staged start's directories that the run reads, for the provenance.
fn digest_existing(option: &str, dir: &Path, names: &[&str]) -> Result<Input, String> {
    let files: Vec<PathBuf> = names.iter().map(|n| dir.join(n)).filter(|p| p.is_file()).collect();
    Input::digest(option, dir, &files.iter().map(|p| p.as_path()).collect::<Vec<_>>()).map_err(|e| e.to_string())
}

/// The seeds the chain reproduces, as the manifest records them.
pub fn seeds_json(p: &ChainParams) -> Json {
    let opt = |v: Option<i32>| v.map(|x| Json::Int(x as i64)).unwrap_or(Json::Null);
    Json::obj(vec![
        ("SCTransform seed.use", opt(p.sct.seed_use)),
        ("RunPCA seed.use", opt(p.pca.seed_use)),
        ("FindClusters random.seed", Json::Int(p.cluster.random_seed as i64)),
        ("GroupSingletons", Json::Str("set.seed(1) before each draw".into())),
        ("Annoy", Json::Str("Kiss64Random's default seed, 1234567890987654321 (Seurat never calls setSeed)".into())),
    ])
}

/// The export's README.txt: what it is, then every file.
fn export_intro() -> Vec<String> {
    vec![
        "This directory is an export of scrust cluster: the input counts, SCTransform's normalised data and the clustering, with one row per cell in metadata.tsv. It reproduces the workflow's chain as Seurat 5.5.1 runs it (SCTransform, RunPCA, FindNeighbors, FindClusters, and bluster's silhouette widths and neighbourhood purities).".into(),
        "THE ONLY RAW COUNTS ARE IN counts/, the input copied byte for byte (its checksums are in summary.json): the export's only 10x directory. sct/ holds SCTransform's normalised data, which are NOT raw counts: the corrected counts (integers, corrected to the median sequencing depth by sctransform's correct()), their log1p and scale.data. No file outside counts/ has a name a 10x reader looks for (matrix.mtx, barcodes.tsv, features.tsv, genes.tsv), so Seurat's Read10X and scanpy's read_10x_mtx fail on sct/ and on this top directory, and each normalised matrix says in its comment lines what it holds and how it was made.".into(),
        "summary.json's files section is the manifest: every file's role, its scale (raw_counts, corrected_counts, log1p_corrected_counts, residuals, ...), its shape, where its row and column names are, its size and SHA-256, and how the export was made (scrust's version, the command, the input's checksums, the R functions reproduced, the parameters and the seeds). The schema is described in docs/output_manifest.md of the scRust repository.".into(),
        "scRust's loaders (loaders/ in the scRust repository) build objects from an export, each matrix in a slot meant for its scale, through the manifest: in R, scrustio::load_seurat(\"<this directory>\") (a Seurat object: the raw counts in the RNA assay, the SCT assay with its SCTModel, the metadata, the PCA, the graphs and the identities) and scrustio::load_sce() (a SingleCellExperiment); in Python, scrustio.load_anndata(\"<this directory>\") (AnnData). As one line: Rscript -e 'scrustio::save_rds(\"<this directory>\", \"object.rds\")', or python -m scrustio h5ad <this directory> object.h5ad. Each file below also loads with standard readers, as shown; counts/ loads as a 10x directory: Seurat::Read10X(\"counts\"), scanpy.read_10x_mtx(\"counts\").".into(),
    ]
}

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// The parameters from `--preset`, `--params` and the flags.
pub fn params_from_args(args: &Args) -> Result<ChainParams, String> {
    let toml = match args.get("params") {
        Some(p) => Some(std::fs::read_to_string(p).map_err(|e| format!("--params {p}: {e}"))?),
        None => None,
    };
    let toml_dir = args.get("params").map(|p| Path::new(p).parent().unwrap_or(Path::new("")).to_path_buf());
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
    chain_params::build_in(args.get("preset"), toml.as_deref(), toml_dir.as_deref(), &flags)
}

/// Whether a covariate of `vars.to.regress` is a column of `--metadata` (M15b): any name but
/// percent.mt, nCount_RNA and nFeature_RNA.
pub fn metadata_regressed(params: &ChainParams) -> bool {
    params.sct.vars_to_regress.iter().any(|v| !FROM_COUNTS.contains(&v.as_str()))
}

/// What the chain needs besides the counts for `percent.mt` and `vars.to.regress` (M15b), read
/// here so that the chain reads no file: `--mt-genes`' ids, with the file's digest, and, when a
/// covariate is one of `--metadata`'s columns, the table's columns for `barcodes` (each cell's
/// text). Returns them and the inputs read, for the provenance.
pub fn read_covariates(params: &ChainParams, metadata: Option<&str>, barcodes: &[String]) -> Result<(CellCovariates, Vec<Input>), String> {
    let mut cov = CellCovariates::default();
    let mut inputs = Vec::new();
    if let Some(f) = &params.mt_genes {
        let p = Path::new(f);
        let text = std::fs::read_to_string(p).map_err(|e| format!("--mt-genes {f}: {e}"))?;
        let input = Input::digest("mt-genes", p.parent().unwrap_or(Path::new("")), &[p]).map_err(|e| format!("--mt-genes {f}: {e}"))?;
        cov.mt_genes = Some(MtGenes::from_file_text(&text, f, &input.files[0].sha256));
        inputs.push(input);
    }
    if metadata_regressed(params) {
        let m = metadata.ok_or_else(|| {
            let cols: Vec<&str> = params.sct.vars_to_regress.iter().map(String::as_str).filter(|v| !FROM_COUNTS.contains(v)).collect();
            format!("--vars-to-regress {}: not percent.mt, nCount_RNA or nFeature_RNA, so a column of --metadata, but no --metadata table is given", cols.join(","))
        })?;
        cov.metadata = Some(cluster::metadata::extra_columns(Path::new(m), barcodes, &[])?);
        inputs.push(Input::digest("metadata", Path::new(m).parent().unwrap_or(Path::new("")), &[Path::new(m)]).map_err(|e| e.to_string())?);
    }
    Ok((cov, inputs))
}

/// The options left out of the normalised matrices' comment lines' command (M15b): `--metadata`
/// changes values when its columns are regressed, so it stays in the command then.
pub fn output_only(params: &ChainParams) -> Vec<&'static str> {
    OUTPUT_ONLY.iter().copied().filter(|&o| !(o == "metadata" && metadata_regressed(params))).collect()
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
    // percent.mt (M15b): only when asked for, so that a run without it records what earlier runs did
    if let Json::Object(m) = &mut sct {
        if p.percent_mt {
            m.push(("percent_mt".into(), Json::Bool(true)));
        }
        if let Some(f) = &p.mt_genes {
            m.push(("mt_genes".into(), Json::Str(f.clone())));
        }
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

/// How a run of the chain writes its output: the options that shape the files (the same for
/// every route).
#[derive(Debug, Clone)]
pub struct Settings {
    pub params: ChainParams,
    pub threads: usize,
    /// `--keep-stages`: `sct/` and every random start's partition.
    pub keep: bool,
    /// `--export` (from the counts only).
    pub export: bool,
    pub dumps: Option<DumpLevel>,
    pub purity_lists: bool,
    /// `--annoy yes`: `neighbors/annoy.ann`.
    pub annoy: bool,
}

impl Settings {
    /// What the chain computes and keeps for these files ([`chain::cluster_counts`]'s result keeps
    /// every stage the files hold).
    pub fn options(&self) -> Options {
        Options {
            keep_partitions: self.keep,
            keep_residuals: self.dumps == Some(DumpLevel::Full),
            keep_purity_lists: self.purity_lists,
            correct: self.export,
            keep_reductions: true,
            keep_sct: self.keep || self.export,
            guard_context: None,
        }
    }
}

/// What a run's `summary.json` records of how it was asked for.
#[derive(Debug, Clone)]
pub struct Record {
    /// The command line as given (`provenance.command_line`).
    pub command_line: Vec<String>,
    /// The command as the normalised matrices' comment lines give it (`provenance.command`).
    pub command: String,
    pub start_time: String,
    /// `run.inputs`: each input option and its path.
    pub inputs: Vec<(&'static str, Json)>,
    /// `provenance.inputs`: each input's files, with their sizes and checksums.
    pub provenance_inputs: Vec<Input>,
    /// With `--export`: `counts/`'s entries in the manifest.
    pub counts_entries: Vec<Entry>,
    /// Seconds spent copying `counts/` (`run.export_counts_copy_s`).
    pub export_s: f64,
    /// `--metadata`: a table whose columns `metadata.tsv` carries over, by barcode.
    pub metadata: Option<String>,
    /// Notes on the scores (`run.stats.note`), such as scores asked for that a route cannot make.
    pub stats_note: Vec<String>,
    /// Warnings about `percent.mt` and the regression (`run.notes`, written only when there are
    /// some; M15b).
    pub notes: Vec<String>,
}

/// The input counts of a run from a 10x directory, as the manifest records them: with `export`,
/// `counts/` (the directory's three files copied into `out` byte for byte, with their checksums)
/// and its entries; otherwise the files' checksums. Returns the provenance input, `counts/`'s
/// entries and the seconds the copy took.
pub fn counts_input(dir: &Path, out: &Path, export: bool) -> Result<(Input, Vec<Entry>, f64), String> {
    let files = tenx::input_files(dir).map_err(|e| e.to_string())?;
    if !export {
        let input = Input::digest("counts", dir, &files.iter().map(|p| p.as_path()).collect::<Vec<_>>()).map_err(|e| e.to_string())?;
        return Ok((input, Vec::new(), 0.0));
    }
    let t = std::time::Instant::now();
    let copied = scio::copy::copy_files(&files.iter().map(|p| p.as_path()).collect::<Vec<_>>(), &out.join("counts")).map_err(|e| e.to_string())?;
    let export_s = t.elapsed().as_secs_f64();
    let names: Vec<String> = copied.iter().map(|c| c.name.clone()).collect();
    let input = Input { option: "counts".into(), path: dir.display().to_string(), files: copied.into_iter().map(|c| InputFile { name: c.name, bytes: c.bytes, sha256: c.sha256 }).collect() };
    let (m, b, f) = (format!("counts/{}", names[0]), format!("counts/{}", names[1]), format!("counts/{}", names[2]));
    let entries = vec![
        Entry::new(&m, "raw_counts", "raw_counts", "matrix_market", "the input counts, copied byte for byte: genes x cells, Cell Ranger's layout (the export's only 10x directory, and its only raw counts)").rows(Ids::numbered("features", &f, 2)).columns(Ids::numbered("cells", &b, 1)),
        Entry::new(&b, "raw_counts_cells", "none", "tsv", "the input counts' barcodes, one per line (their columns)").header(false),
        Entry::new(&f, "raw_counts_features", "none", "tsv", "the input counts' features: id, name, type (their rows); Seurat's rows are the names").header(false),
    ];
    Ok((input, entries, export_s))
}

/// What a run of the chain from the counts leaves to its caller, besides its files.
#[derive(Debug, Clone)]
pub struct Done {
    /// FindClusters' result and the scores (the stages are freed once written).
    pub chain: Chain,
    /// The run's phases (`run.phases`: name, seconds since `started`, peak RSS), for a caller
    /// that records them on its own clock.
    pub phases: Vec<(String, f64, Option<u64>)>,
    pub started: std::time::Instant,
}

/// The chain from counts in memory, written into `out` as `scrust cluster --counts` writes it: the
/// entry point ([`chain::cluster_counts`]; docs/PLAN.md A3, C3h), then its result written by the
/// same code as every route, the manifest and `summary.json`. `scrust cluster --counts` reads the
/// counts and calls this; `scrust doubletfinder --recluster` hands over its singlets from memory
/// (M15). `phases` holds the caller's marks so far (`load`); `guard_context` is the PCA guard's
/// (`doubletfinder.R` calls it with "doubletfinder").
pub fn run_counts(counts: scio::Counts, out: &Path, s: &Settings, mut rec: Record, phases: Phases, guard_context: Option<&'static str>, cov: &CellCovariates) -> Result<Done, String> {
    if let Some(e) = rec.counts_entries.first_mut() {
        *e = e.clone().shape(&[counts.n_features(), counts.n_cells()]).nnz(counts.matrix.nnz());
    }
    std::fs::create_dir_all(out).map_err(|e| format!("{}: {e}", out.display()))?;
    let mut w = Writer::new(out, counts.barcodes.clone(), s.threads);
    // the files come from the result; the writer observes only the stage dumps and the timings
    w.write_stages = false;
    w.keep_stages = s.keep;
    w.export = s.export;
    w.command = rec.command.clone();
    w.input_checksums = cli::input_checksums(&rec.provenance_inputs);
    w.sct_bin_size = s.params.sct.bin_size;
    w.annoy = s.annoy;
    w.purity_lists = s.purity_lists;
    w.neighbor_params = Some(s.params.neighbors.clone());
    if let Some(level) = s.dumps {
        w = w.with_dumps(level);
    }
    let opts = Options { guard_context, ..s.options() };
    let c = chain::cluster_counts_with(Cow::Owned(counts), &s.params, cov, s.threads, opts, &mut w).map_err(|e| e.to_string())?;
    for n in &c.notes {
        eprintln!("scrust cluster: {n}");
    }
    rec.notes.extend(c.notes.iter().cloned());
    w.write_clustered(&c).map_err(|e| e.to_string())?;
    // the stages are written: only the clustering and the scores go on
    let r = c.chain;
    finish_run(out, "counts", s, rec, w, r, phases)
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let mut switches: Vec<&str> = PARAM_SWITCHES.to_vec();
    switches.extend(["keep-stages", "export"]);
    let args = Args::parse_with_flags(argv, &switches)?;
    let mut known: Vec<&str> = vec!["counts", "sct", "pca", "snn", "out", "preset", "params", "metadata", "keep-stages", "export", "dump-stages", "purity-lists", "annoy", "threads"];
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
    // neighborPurity's neighbour lists in the dumps, as R's instrumented runs choose them
    let purity_lists = match args.get("purity-lists").unwrap_or("auto") {
        "auto" => dumps == Some(DumpLevel::Full),
        "yes" => dumps.is_some(),
        "no" => false,
        v => return Err(format!("--purity-lists: expected auto, yes or no, not '{v}'")),
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
    let export = args.flag("export");
    if export && start != "counts" {
        return Err("--export needs --counts: the export copies the input counts (counts/) and makes the corrected counts from vst's residuals of every model gene, which only a run from the counts has (scrust sctransform keeps only the variable features' clipped residuals); run scrust cluster --counts <dir> --export, or scrust sctransform --corrected-counts for the SCT part alone".into());
    }
    if start != "counts" && params.needs_covariates() {
        return Err("--vars-to-regress, --mt-genes and --percent-mt act where SCTransform runs: give them to scrust cluster --counts, or to scrust sctransform, whose cell_metadata.tsv carries percent.mt on to scrust pca, scrust neighbors and the staged starts (the regression is in its scale_data.npy)".into());
    }
    let command = scio::manifest::canonical_command("cluster", argv, &output_only(&params), &switches);
    if out.join("summary.json").exists() || out.join("clusters.tsv").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }
    let settings = Settings { params, threads, keep, export, dumps, purity_lists, annoy };
    let opts = settings.options();
    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let mut rec = Record {
        command_line: cli::command_line("cluster", argv),
        command,
        start_time,
        inputs: Vec::new(),
        provenance_inputs: Vec::new(),
        counts_entries: Vec::new(),
        export_s: 0.0,
        metadata: args.get("metadata").map(str::to_string),
        stats_note: Vec::new(),
        notes: Vec::new(),
    };
    let params = &settings.params;

    let (r, w) = match start {
        "counts" => {
            // read, call the entry point, write (docs/PLAN.md A3): run_counts
            let dir = args.path("counts")?;
            let (input, entries, export_s) = counts_input(&dir, &out, export)?;
            rec.provenance_inputs.push(input);
            rec.counts_entries = entries;
            rec.export_s = export_s;
            rec.inputs.push(("counts", Json::Str(dir.display().to_string())));
            let counts = scio::Counts::read(&dir).map_err(|e| e.to_string())?;
            // percent.mt's ids and --metadata's regressed columns (M15b), before the chain
            let (cov, more) = read_covariates(params, args.get("metadata"), &counts.barcodes)?;
            rec.provenance_inputs.extend(more);
            phases.mark("load");
            run_counts(counts, &out, &settings, rec, phases, None, &cov)?;
            return Ok(());
        }
        "sct" => {
            let dir = args.path("sct")?;
            let barcodes = table_barcodes(&dir.join("cell_attr.tsv"))?;
            let vf = Table::read(&dir.join("variable_features.tsv")).map_err(|e| format!("{}: {e}", dir.join("variable_features.tsv").display()))?;
            let (rows, nf, nc) = scio::npy::read_f64_row_major(&dir.join("scale_data.npy")).map_err(|e| e.to_string())?;
            if nf != vf.rows.len() || nc != barcodes.len() {
                return Err(format!("{}: scale_data.npy is {nf} x {nc}, with {} variable features and {} cells", dir.display(), vf.rows.len(), barcodes.len()));
            }
            rec.inputs.push(("sct", Json::Str(dir.display().to_string())));
            rec.provenance_inputs.push(digest_existing("sct", &dir, &["cell_attr.tsv", "variable_features.tsv", "scale_data.npy", "summary.json", cluster::metadata::FILE])?);
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
            let r = chain::from_scale_data(rows, (0..nf).collect(), nc, params, threads, opts, &mut w);
            (r, w)
        }
        "pca" => {
            let dir = args.path("pca")?;
            let emb = crate::cmd_neighbors::read_embeddings(&dir)?;
            rec.inputs.push(("pca", Json::Str(dir.display().to_string())));
            rec.provenance_inputs.push(digest_existing("pca", &dir, &["embeddings.tsv", "embeddings.npy", cluster::metadata::FILE])?);
            let rna = carried_rna(&[&dir], &emb.barcodes)?;
            let mut w = Writer::new(&out, emb.barcodes.clone(), threads);
            w.rna = rna;
            w.keep_stages = keep;
            w.annoy = annoy;
            w.neighbor_params = Some(params.neighbors.clone());
            w.pca_result = w.copy_stage(&dir, "pca", &pca::output::FILES).map_err(|e| e.to_string())?;
            phases.mark("load");
            let r = chain::from_embeddings(emb.data, emb.n, emb.ncol, params, threads, opts, &mut w);
            (r, w)
        }
        _ => {
            let dir = args.path("snn")?;
            let input = cluster::input::read_snn_dir(&dir)?;
            rec.inputs.push(("snn", Json::Str(input.source.display().to_string())));
            rec.inputs.push(("barcodes", Json::Str(input.barcodes_source.display().to_string())));
            let mut read: Vec<&Path> = vec![input.source.as_path()];
            if input.barcodes_source != input.source {
                read.push(input.barcodes_source.as_path());
            }
            rec.provenance_inputs.push(Input::digest("snn", &dir, &read).map_err(|e| e.to_string())?);
            let emb = match args.get("pca") {
                Some(p) => {
                    let e = crate::cmd_neighbors::read_embeddings(Path::new(p))?;
                    if e.barcodes != input.barcodes {
                        return Err(format!("--pca {p}: its cells are not the SNN graph's"));
                    }
                    rec.inputs.push(("pca", Json::Str(p.to_string())));
                    rec.provenance_inputs.push(digest_existing("pca", Path::new(p), &["embeddings.tsv", "embeddings.npy", cluster::metadata::FILE])?);
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
            for (asked, name) in [(params.stats.silhouette, "silhouette"), (params.stats.purity, "purity")] {
                if asked && emb.is_none() {
                    let m = format!("{name} asked for, but not computed: the scores need the PCA (give --pca with --snn)");
                    eprintln!("scrust cluster: {m}");
                    rec.stats_note.push(m);
                }
            }
            phases.mark("load");
            let r = chain::from_snn(&input.snn, emb.as_ref().map(|(_, e)| (e.data.as_slice(), e.ncol)), params, threads, opts, &mut w);
            (r, w)
        }
    };
    let r = r.map_err(|e| e.to_string())?;
    finish_run(&out, start, &settings, rec, w, r, phases)?;
    Ok(())
}

/// The tail of every route: what FindClusters printed, the labels' files and `metadata.tsv` (with
/// `--metadata`'s columns), the stage dumps, the manifest (with `--export`, `counts/` and the
/// README) and `summary.json`. `w` has written or copied the stages; `r` is the clustering.
pub fn finish_run(out: &Path, start: &str, s: &Settings, rec: Record, mut w: Writer, r: Chain, mut phases: Phases) -> Result<Done, String> {
    let params = &s.params;
    let threads = s.threads;
    let mut inputs = rec.inputs;
    let mut provenance_inputs = rec.provenance_inputs;
    let mut stats_note = rec.stats_note;
    // the chain's stages, timed from the writer's making, on the command's clock
    let offset = w.started().saturating_duration_since(phases.started()).as_secs_f64();
    for (name, t, rss) in &w.phases {
        phases.marks.push((name.clone(), *t + offset, *rss));
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
    if let Some(m) = &rec.metadata {
        let mut taken: Vec<&str> = vec!["barcode", "nCount_RNA", "nFeature_RNA", "seurat_clusters", "silhouette_width", "silhouette_other", "neighborhood_purity", "purity_maximum", "nCount_SCT", "nFeature_SCT"];
        if w.rna.as_ref().is_some_and(|x| x.percent_mt.is_some()) {
            taken.push(cluster::covariates::PERCENT_MT);
        }
        let cols: Vec<String> = r.clusters.resolutions.iter().map(|x| x.column.name.clone()).collect();
        taken.extend(cols.iter().map(String::as_str));
        w.extra = Some(cluster::metadata::extra_columns(Path::new(m), &w.barcodes, &taken)?);
        inputs.push(("metadata", Json::Str(m.to_string())));
        // read before the chain when its columns are regressed (M15b): digested once
        if !provenance_inputs.iter().any(|i| i.option == "metadata") {
            provenance_inputs.push(Input::digest("metadata", Path::new(m).parent().unwrap_or(Path::new("")), &[Path::new(m)]).map_err(|e| e.to_string())?);
        }
    }
    w.finish_labels(&r.clusters, r.silhouette.as_ref(), r.purity.as_ref()).map_err(|e| e.to_string())?;
    w.finish_dumps().map_err(|e| e.to_string())?;
    phases.mark("write");

    // the manifest (C3f): every file, its role and scale, and how the output was made; with
    // --export, counts/ and README.txt too
    let t_manifest = std::time::Instant::now();
    let mut manifest = Manifest::new(
        "cluster",
        Provenance {
            scrust_version: env!("CARGO_PKG_VERSION").into(),
            command_line: rec.command_line,
            command: rec.command,
            inputs: provenance_inputs,
            reproduces: vec![
                sctransform::REFERENCE.into(),
                pca::REFERENCE.into(),
                neighbors::REFERENCE.into(),
                cluster::REFERENCE.into(),
                "bluster 1.20.0 approxSilhouette and neighborPurity (BiocNeighbors 2.4.0, knncolle's KMKNN: searched exhaustively with its arithmetic), as the workflow's add_silhouette_to_metadata calls them; the workflow's write_cluster_metadata".into(),
            ],
            parameters: params_json(params),
            seeds: seeds_json(params),
        },
    );
    manifest.push(Manifest::summary_entry());
    manifest.entries.extend(w.manifest_entries().map_err(|e| e.to_string())?);
    manifest.entries.extend(rec.counts_entries);
    if s.export {
        manifest.push(Entry::new("README.txt", "readme", "none", "text", "this file: what the export holds and how to load each file"));
    }
    manifest.fill_shapes(out).map_err(|e| e.to_string())?;
    if s.export {
        let text = scio::manifest::readme(&manifest, "scrust cluster --export", &export_intro());
        write(&out.join("README.txt"), &text)?;
    }
    let files_json = cli::finish_manifest(out, &mut manifest, threads)?;
    let manifest_s = t_manifest.elapsed().as_secs_f64();
    phases.mark("manifest");

    let mut computed: Vec<&str> = Vec::new();
    if r.silhouette.is_some() {
        computed.push("silhouette");
    }
    if let Some(p) = &r.purity {
        computed.push("purity");
        if p.k < p.k_requested {
            // BiocNeighbors' warning
            let m = format!("neighborPurity: 'k' capped at the number of observations minus 1 ({})", p.k);
            eprintln!("scrust cluster: {m}");
            stats_note.push(m);
        }
        let c = &p.certificate;
        if c.n_flagged() > 0 {
            let m = format!(
                "neighborPurity's certificate flags {} pair(s) ({} the median pair, on the radius by construction; {} other), reported in summary.json (result.purity.certificate): KMKNN's floating-point bounds could treat them differently from this exhaustive search, whose result is kept",
                c.n_flagged(),
                c.n_flagged_kind(knn::purity::Kind::MedianPair),
                c.n_flagged_kind(knn::purity::Kind::Other)
            );
            eprintln!("scrust cluster: {m}");
            stats_note.push(m);
        }
    }
    if params.stats.names().is_empty() {
        stats_note.push("none asked for".into());
    }
    let version = env!("CARGO_PKG_VERSION");
    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let t = &r.timings;
    let result = w.result_json(&r.clusters);
    let notes = rec.notes;
    let mut summary = Json::obj(vec![
        ("tool", Json::Str("scrust cluster".into())),
        ("reference", Json::Str(format!("{}; {}; {}; {}; bluster 1.20.0 approxSilhouette and neighborPurity (BiocNeighbors 2.4.0, knncolle's KMKNN: searched exhaustively with its arithmetic, with a certificate)", sctransform::REFERENCE, pca::REFERENCE, neighbors::REFERENCE, cluster::REFERENCE))),
        ("scrust_version", Json::Str(version.into())),
        ("result", result),
        (
            "run",
            Json::obj(vec![
                ("parameters", params_json(params)),
                ("start", Json::Str(start.into())),
                ("inputs", Json::obj(inputs)),
                ("outdir", Json::Str(out.display().to_string())),
                ("files", Json::strs(&w.files)),
                (
                    "stats",
                    Json::obj(vec![("requested", Json::strs(&params.stats.names())), ("computed", Json::strs(&computed)), ("note", Json::Str(stats_note.join("; ")))]),
                ),
                ("keep_stages", Json::Bool(s.keep)),
                ("export", Json::Bool(s.export)),
                ("export_counts_copy_s", Json::Num(rec.export_s)),
                ("manifest_s", Json::Num(manifest_s)),
                ("dump_stages", Json::Str(s.dumps.map_or("none", |d| d.name()).into())),
                ("purity_lists", Json::Bool(s.purity_lists)),
                ("threads", Json::Int(threads as i64)),
                ("start_time", Json::Str(rec.start_time)),
                ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
                ("sctransform_s", Json::Num(t.sctransform)),
                ("pca_s", Json::Num(t.pca)),
                ("neighbors_s", Json::Num(t.neighbors)),
                ("cluster_s", Json::Num(t.cluster)),
                ("silhouette_s", Json::Num(t.silhouette)),
                ("purity_s", Json::Num(t.purity)),
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
        ("files", files_json),
    ]);
    // M15b: the warnings about percent.mt and the regression, only when there are some
    if !notes.is_empty() {
        if let Json::Object(m) = &mut summary {
            if let Some((_, Json::Object(run))) = m.iter_mut().find(|(k, _)| k == "run") {
                run.push(("notes".into(), Json::strs(&notes)));
            }
        }
    }
    write(&out.join("summary.json"), &scio::json::to_file_text(&summary))?;
    eprint!("scrust cluster: {} cells; {}", r.n_cells, cluster::output::describe(&r.clusters));
    if let Some(p) = r.purity.as_ref().filter(|p| p.computed) {
        let c = &p.certificate;
        eprintln!(
            "scrust cluster: neighborPurity: radius {}, {:.2} neighbours per cell; certificate: {}",
            scio::fmt17(p.radius),
            p.mean_neighbors(),
            if c.n_flagged() == 0 {
                "no pair flagged".to_string()
            } else {
                format!("{} pair(s) flagged and reported ({} the median pair, {} other)", c.n_flagged(), c.n_flagged_kind(knn::purity::Kind::MedianPair), c.n_flagged_kind(knn::purity::Kind::Other))
            }
        );
    }
    eprintln!(
        "scrust cluster: SCTransform {:.1} s, RunPCA {:.1} s, FindNeighbors {:.1} s, FindClusters {:.1} s, silhouette {:.2} s, purity {:.2} s; {:.1} s in all, {threads} threads",
        t.sctransform,
        t.pca,
        t.neighbors,
        t.cluster,
        t.silhouette,
        t.purity,
        phases.elapsed()
    );
    Ok(Done { chain: r, started: phases.started(), phases: phases.marks })
}
