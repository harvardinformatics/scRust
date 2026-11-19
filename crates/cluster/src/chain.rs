//! The workflow's clustering chain as one library call (docs/PLAN.md C1, C3 and C3h; milestone
//! M11; purity: M12): SCTransform, the workflow's PCA guard, RunPCA, FindNeighbors, FindClusters,
//! the silhouette widths and the neighbourhood purities, from the counts, or from any stage's
//! output (`scrust cluster --sct`, `--pca`, `--snn`).
//!
//! ```text
//! obj <- SCTransform(obj, seed.use = 12345)            sctransform::start / finish
//! require_min_cells_for_pca(obj)                         the guard: stop at 50 cells or fewer
//! obj <- RunPCA(obj, seed.use = 12345)                   pca::run_pca (npcs 50)
//! obj <- FindNeighbors(obj, dims = 1:30)                 neighbors::find_neighbors
//! obj <- FindClusters(obj, random.seed = 12345)          crate::find_clusters
//! obj <- add_silhouette_to_metadata(obj)                 crate::silhouette, crate::purity
//! ```
//!
//! RunUMAP, which the workflow runs between RunPCA and FindNeighbors, is left out: it changes
//! nothing downstream (C1; M8 checked it on all 34 reference objects).
//!
//! **The entry point** (C3h; A3's in-memory hand-offs; milestone M15) is [`cluster_counts`]: from
//! counts in memory ([`scio::Counts`]: the matrix with its barcodes and features) to a complete
//! typed result, [`Clustered`]. It opens no file. The result holds the labels with their levels,
//! the scores, every stage's summary for `summary.json`, and, as [`Options`] asks, RunPCA's and
//! FindNeighbors' results and SCTransform's model, `scale.data` and corrected counts: everything
//! `scrust cluster --counts`' files hold, which [`crate::chain_output::Writer::write_clustered`]
//! writes. `scrust cluster --counts` and `scrust doubletfinder` (clustering on the fly, and the
//! singlets' re-clustering) call it; the decontX and SoupX ports will too. The staged starts
//! ([`from_scale_data`], [`from_embeddings`], [`from_snn`]) and [`from_counts`] keep their
//! observer-only form.
//!
//! **Memory.** Each stage is freed once the next has what it needs: the counts once
//! SCTransform's setup has kept the kept genes' rows (their per-cell totals are taken first),
//! SCTransform's working data once `scale.data` is made, `scale.data` once the PCA is done (unless
//! the result keeps SCTransform's model), FindNeighbors' index and lists once the SNN graph is made
//! (unless the result keeps them). [`Observer`]s see each stage's result by reference as it is
//! made (to record its stage dumps, or to write it on the staged starts) and keep nothing.
//!
//! **Random numbers.** R's stream (`&mut RRng`) is reseeded by SCTransform and by RunPCA, as in
//! R; Annoy and the modularity optimizer have their own generators; GroupSingletons reseeds R's
//! with `set.seed(1)`. The thread count changes no result.
//!
//! **`vars.to.regress` and `percent.mt`** (milestone M15b; [`crate::covariates`]):
//! [`cluster_counts_with`] takes what the parameters' covariates need in memory
//! ([`CellCovariates`]: `--mt-genes`' ids, `--metadata`'s columns), computes `percent.mt` from
//! the counts when asked (before SCTransform, as the workflow's scripts do), takes the QR of the
//! covariates once and hands it to SCTransform, which regresses `scale.data`'s clipped rows on
//! them before centring them. Without them the code path is the one before M15b, unchanged.

use crate::chain_params::ChainParams;
use crate::covariates::{covariate_values, CellCovariates, MtGenes, PercentMt};
use crate::find_clusters::{find_clusters, ClusterError, Factor, FindClusters};
use crate::metadata::RnaMetadata;
use crate::purity::{neighbor_purity, Purity};
use crate::silhouette::{approx_silhouette, Silhouette};
use neighbors::{find_neighbors, Neighbors, NeighborsError};
use pca::{Pca, PcaError};
use rcompat::RRng;
use scio::{CscMatrix, Json};
use sctransform::{Regression, Sct, SctError};
use std::borrow::Cow;
use std::fmt;
use std::io;
use std::time::Instant;

/// Why the chain stops.
#[derive(Debug)]
pub enum ChainError {
    Params(String),
    Sct(SctError),
    /// The workflow's `require_min_cells_for_pca()`, with its message.
    Guard(String),
    Pca(PcaError),
    Neighbors(NeighborsError),
    Cluster(ClusterError),
    Silhouette(String),
    Purity(String),
    Io(io::Error),
}

impl fmt::Display for ChainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChainError::Params(m) | ChainError::Guard(m) => f.write_str(m),
            ChainError::Sct(e) => write!(f, "SCTransform: {e}"),
            ChainError::Pca(e) => write!(f, "RunPCA: {e}"),
            ChainError::Neighbors(e) => write!(f, "FindNeighbors: {e}"),
            ChainError::Cluster(e) => write!(f, "FindClusters: {e}"),
            ChainError::Silhouette(m) | ChainError::Purity(m) => write!(f, "add_silhouette_to_metadata: {m}"),
            ChainError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for ChainError {}

impl From<io::Error> for ChainError {
    fn from(e: io::Error) -> Self {
        ChainError::Io(e)
    }
}

/// The workflow's `require_min_cells_for_pca(seurat_obj, context = "", npcs = 50L)`: R stops
/// with this message when there are `npcs` cells or fewer.
pub fn pca_guard(n_cells: usize, npcs: usize) -> Result<(), ChainError> {
    pca_guard_in(n_cells, npcs, None)
}

/// [`pca_guard`] with the call's `context`, which R puts before the message (`doubletfinder.R`
/// calls it with `context = "doubletfinder"`).
pub fn pca_guard_in(n_cells: usize, npcs: usize, context: Option<&str>) -> Result<(), ChainError> {
    if n_cells <= npcs {
        let prefix = context.filter(|c| !c.is_empty()).map(|c| format!("{c}: ")).unwrap_or_default();
        return Err(ChainError::Guard(format!(
            "{prefix}[LOW_QUALITY_SAMPLE] only {n_cells} cell(s) remain - too few to compute {npcs} principal components (RunPCA and downstream UMAP/clustering require more cells than PCs). This usually means upstream QC/filtering removed nearly all cells for this sample; consider dataset-specific thresholds or excluding this sample."
        )));
    }
    Ok(())
}

/// Sees each stage's result as the chain makes it. Every method does nothing by default.
pub trait Observer {
    /// SCTransform's result (its RNG states at `vst`'s entry and SCTransform's exit included).
    fn sct(&mut self, _sct: &Sct) -> io::Result<()> {
        Ok(())
    }
    /// RunPCA's result: the input's rows (rows of the counts, in `VariableFeatures` order), the
    /// input (row-major, rows x cells), the generator's states at RunPCA's entry (`None` when the
    /// chain started from SCTransform's output, where R's state is not known) and exit.
    fn pca(&mut self, _input_rows: &[usize], _input: &[f64], _p: &Pca, _entry: Option<&[i32]>, _exit: &[i32]) -> io::Result<()> {
        Ok(())
    }
    /// FindNeighbors' result, on embeddings of `ncol` columns.
    fn neighbors(&mut self, _r: &Neighbors, _ncol: usize) -> io::Result<()> {
        Ok(())
    }
    fn clusters(&mut self, _r: &FindClusters) -> io::Result<()> {
        Ok(())
    }
    fn silhouette(&mut self, _s: &Silhouette) -> io::Result<()> {
        Ok(())
    }
    fn purity(&mut self, _p: &Purity) -> io::Result<()> {
        Ok(())
    }
    /// A stage ended (`sctransform`, `pca`, `neighbors`, `cluster`, `silhouette`, `purity`).
    fn phase(&mut self, _name: &str) {}
}

/// An observer that keeps nothing.
pub struct NoObserver;

impl Observer for NoObserver {}

/// What a run computes and keeps beyond the labels and the scores: every random start's partition
/// (FindClusters), the full residual matrix (for the `full` stage dumps; small inputs only),
/// purity's neighbour lists (for the `standard` and `full` stage dumps), SCTransform's corrected
/// counts (`correct()`, for the export; M13), which the observer sees with SCTransform's result;
/// and what [`cluster_counts`]'s result keeps (M15).
#[derive(Debug, Clone, Copy, Default)]
pub struct Options {
    pub keep_partitions: bool,
    pub keep_residuals: bool,
    pub keep_purity_lists: bool,
    pub correct: bool,
    /// Keep RunPCA's and FindNeighbors' results in [`Clustered`] (what `pca/` and `neighbors/`
    /// hold). Without it they are freed once the next stage has what it needs, and the result
    /// holds the labels, the scores and every stage's summary: clustering on the fly (C3h).
    pub keep_reductions: bool,
    /// Keep SCTransform's model and `scale.data` in [`Clustered`] (`sct/`: `--keep-stages`,
    /// `--export`), with the corrected counts when `correct` makes them. The kept genes' counts,
    /// the residual matrix and the fits' working data are freed all the same.
    pub keep_sct: bool,
    /// The context the workflow's PCA guard names in its message ([`pca_guard_in`]).
    pub guard_context: Option<&'static str>,
}

/// Wall-clock seconds of each stage the run made (0 for stages it started after).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Timings {
    pub sctransform: f64,
    pub pca: f64,
    pub neighbors: f64,
    pub cluster: f64,
    pub silhouette: f64,
    pub purity: f64,
}

/// The chain's result: the clustering and the scores. The stages before it were seen by the
/// observer and are gone.
#[derive(Debug, Clone)]
pub struct Chain {
    pub n_cells: usize,
    pub clusters: FindClusters,
    /// `None` when not asked for, or when the chain started from the SNN graph without
    /// embeddings.
    pub silhouette: Option<Silhouette>,
    /// Likewise.
    pub purity: Option<Purity>,
    /// The embeddings' columns (the PCA's components).
    pub n_components: Option<usize>,
    pub timings: Timings,
}

/// SCTransform's part of [`Clustered`]: its summary always, its model when kept.
#[derive(Debug, Clone)]
pub struct SctResult {
    /// `summary.json`'s `result.sctransform` ([`sctransform::output::result_json`]).
    pub summary: Json,
    /// R's `.Random.seed` when SCTransform returns.
    pub rng_exit: Vec<i32>,
    /// R's digest of `scale.data` (`result.scale_data_sha256`).
    pub scale_data_sha256: String,
    /// `nCount_SCT` and `nFeature_SCT`, when the corrected counts were made (`correct`).
    pub sct_columns: Option<(Vec<f64>, Vec<u32>)>,
    /// The model, `scale.data` and the corrected counts (`sct/`), when [`Options::keep_sct`].
    /// The kept genes' counts (`rows`) and the full residual matrix are freed.
    pub kept: Option<Sct>,
}

/// RunPCA's part of [`Clustered`].
#[derive(Debug, Clone)]
pub struct PcaResult {
    /// `result.pca` ([`pca::output::result_json`]).
    pub summary: Json,
    /// R's `.Random.seed` when RunPCA returns (`pca/random_seed.txt`).
    pub rng_exit: Vec<i32>,
    /// The PCA's input rows (rows of the counts, in `VariableFeatures` order: the loadings' rows)
    /// and the result (`pca/`), when [`Options::keep_reductions`].
    pub kept: Option<(Vec<usize>, Pca)>,
}

/// FindNeighbors' part of [`Clustered`].
#[derive(Debug, Clone)]
pub struct NeighborsResult {
    /// `result.neighbors` ([`neighbors::output::result_json`]).
    pub summary: Json,
    /// The stored entries of the nn and SNN graphs.
    pub nnz: (usize, usize),
    /// The result (`neighbors/`), when [`Options::keep_reductions`].
    pub kept: Option<Neighbors>,
}

/// The chain's complete result from the counts ([`cluster_counts`]; docs/PLAN.md A3's in-memory
/// hand-offs, C3h): everything `scrust cluster --counts`' files hold, as typed values, as far as
/// [`Options`] keeps the stages.
#[derive(Debug, Clone)]
pub struct Clustered {
    /// The cells (the counts' barcodes), in column order.
    pub cells: Vec<String>,
    /// The counts' features (id, name): the rows of `sct/` and of the loadings.
    pub features: Vec<(String, String)>,
    /// `nCount_RNA` and `nFeature_RNA`, as `CreateSeuratObject()` computes them, and `percent.mt`
    /// when asked for (M15b).
    pub rna: RnaMetadata,
    /// `percent.mt`'s definition and record, when computed (M15b; its values are in `rna`).
    pub percent_mt: Option<PercentMt>,
    /// Warnings about `percent.mt` and the regression (ids not found, no feature matched, a
    /// rank-deficient model), for the caller to print and record.
    pub notes: Vec<String>,
    pub sct: SctResult,
    pub pca: PcaResult,
    pub neighbors: NeighborsResult,
    /// FindClusters' result and the scores.
    pub chain: Chain,
}

impl Clustered {
    /// `seurat_clusters` with its levels: what a tool given no labels uses (C3h).
    pub fn labels(&self) -> &Factor {
        &self.chain.clusters.seurat_clusters
    }

    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }
}

/// The chain from counts in memory (C3h's entry point; A3's in-memory hand-offs): SCTransform, the
/// workflow's PCA guard (with the preset), RunPCA, FindNeighbors, FindClusters, and the silhouette
/// widths and purities when asked, into a complete result. It opens no file.
///
/// `counts` is borrowed or taken: taken (`Cow::Owned`), the matrix is freed as soon as
/// SCTransform's setup has kept the rows it needs; borrowed, the caller keeps it (clustering on the
/// fly, where DoubletFinder reads the same counts next). Each later stage is freed once the next
/// has what it needs, unless `opts` keeps it in the result. `obs` sees each stage as it is made
/// (stage dumps, timings); [`NoObserver`] sees nothing.
pub fn cluster_counts(counts: Cow<'_, scio::Counts>, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Clustered, ChainError> {
    cluster_counts_with(counts, params, &CellCovariates::default(), threads, opts, obs)
}

/// [`cluster_counts`] with what `percent.mt` and `vars.to.regress` need besides the counts
/// (M15b): `--mt-genes`' ids and `--metadata`'s columns, in memory (`cov`; the pattern when no
/// ids are given). When the parameters ask for neither, this is [`cluster_counts`].
pub fn cluster_counts_with(counts: Cow<'_, scio::Counts>, params: &ChainParams, cov: &CellCovariates, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Clustered, ChainError> {
    params.check().map_err(ChainError::Params)?;
    let mut rna = RnaMetadata::of(&counts.matrix);
    let (percent_mt, regression, notes) = covariates_for(&counts, params, cov, &mut rna)?;
    let features: Vec<(String, String)> = counts.features().into_iter().map(|f| (f.id, f.name)).collect();
    let (matrix, cells) = match counts {
        Cow::Owned(c) => {
            let scio::Counts { matrix, barcodes, .. } = c;
            (Cow::Owned(matrix), barcodes)
        }
        Cow::Borrowed(c) => (Cow::Borrowed(&c.matrix), c.barcodes.clone()),
    };
    let mut s = stages_from_counts(matrix, params, regression, threads, opts, obs)?;
    if let (Some(p), Json::Object(m)) = (&percent_mt, &mut s.sct.summary) {
        m.push(("percent_mt".into(), p.record_json()));
    }
    Ok(Clustered { cells, features, rna, percent_mt, notes, sct: s.sct, pca: s.pca, neighbors: s.neighbors, chain: s.chain })
}

/// `percent.mt` (when the parameters ask for it) into `rna`, and `vars.to.regress`' regression
/// (when they name covariates), from the counts and `cov` (M15b). Returns `percent.mt`'s record,
/// the regression, and the warnings for the caller.
pub fn covariates_for(counts: &scio::Counts, params: &ChainParams, cov: &CellCovariates, rna: &mut RnaMetadata) -> Result<(Option<PercentMt>, Option<Regression>, Vec<String>), ChainError> {
    let mut notes = Vec::new();
    let vars = &params.sct.vars_to_regress;
    if params.mt_genes.is_some() != cov.mt_genes.is_some() {
        return Err(ChainError::Params("--mt-genes: the parameters name a file, but its ids were not handed to the chain (or the other way round)".into()));
    }
    let percent_mt = if params.wants_percent_mt() {
        let mt = cov.mt_genes.clone().unwrap_or(MtGenes::Pattern);
        let p = PercentMt::compute(counts, &mt).map_err(ChainError::Params)?;
        if p.rows.is_empty() && vars.iter().any(|v| v == crate::covariates::PERCENT_MT) {
            return Err(ChainError::Params(format!(
                "--vars-to-regress percent.mt: no feature's name matches the workflow's pattern \"{}\", so percent.mt would be 0 in every cell; give the mitochondrial genes' ids (column 1 of features.tsv) with --mt-genes <file>",
                scio::features::MT_PATTERN
            )));
        }
        notes.extend(p.notes());
        rna.percent_mt = Some(p.values.clone());
        Some(p)
    } else {
        None
    };
    let regression = if vars.is_empty() {
        None
    } else {
        let names = counts.seurat_names();
        let values = covariate_values(vars, rna, cov.metadata.as_ref(), &names).map_err(ChainError::Params)?;
        let r = Regression::new(vars.clone(), &values, counts.n_cells()).map_err(ChainError::Params)?;
        let (rank, p) = r.rank();
        if rank < p {
            notes.push(format!(
                "vars.to.regress: the model matrix (the intercept and {}) has rank {rank} of {p}: lm()'s QR moved {} column(s) to the end as negligible (pivot {}), and the regression uses the first {rank}, as R's does",
                vars.join(", "),
                p - rank,
                r.qr.pivot.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",")
            ));
        }
        Some(r)
    };
    Ok((percent_mt, regression, notes))
}

/// From the counts (genes x cells, every stored entry as read): the whole chain, keeping only the
/// clustering and the scores (and what `opts` keeps, which this form drops); the observer sees the
/// stages. [`cluster_counts`] returns them.
pub fn from_counts(counts: CscMatrix, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.check().map_err(ChainError::Params)?;
    if params.needs_covariates() {
        return Err(ChainError::Params("percent.mt and vars.to.regress need the counts' features and the covariates: use cluster_counts_with".into()));
    }
    Ok(stages_from_counts(Cow::Owned(counts), params, None, threads, opts, obs)?.chain)
}

struct Stages {
    sct: SctResult,
    pca: PcaResult,
    neighbors: NeighborsResult,
    chain: Chain,
}

/// The chain from the counts, stage by stage, freeing what no later stage or the result needs.
fn stages_from_counts(matrix: Cow<'_, CscMatrix>, params: &ChainParams, regression: Option<Regression>, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Stages, ChainError> {
    let t0 = Instant::now();
    let n_features = matrix.nrow;
    // SCTransform.Seurat calls set.seed(seed.use) first: R's state before it does not matter
    let mut rng = RRng::new(params.sct.seed_use.expect("checked"));
    let started = sctransform::start(&matrix, &params.sct, &mut rng).map_err(ChainError::Sct)?;
    // taken counts are freed here; borrowed ones stay with the caller
    drop(matrix);
    let mut sct = sctransform::finish_full(started, &params.sct, &mut rng, opts.keep_residuals, opts.correct, regression, threads).map_err(ChainError::Sct)?;
    let mut sct_result = SctResult {
        summary: sctransform::output::result_json(&sct, n_features),
        rng_exit: sct.rng_exit.clone(),
        scale_data_sha256: sctransform::output::scale_data_sha256_r_order(&sct.residuals),
        sct_columns: sct.corrected.as_ref().map(|c| (c.n_count(), c.n_feature())),
        kept: None,
    };
    obs.sct(&sct)?;
    let sct_s = t0.elapsed().as_secs_f64();
    obs.phase("sctransform");
    let rows: Vec<usize> = sct.residuals.top_features.iter().map(|&k| sct.setup.genes[k]).collect();
    let n = sct.residuals.ncells;
    // only scale.data goes on (and the model, when kept): the kept genes' counts, the full
    // residual matrix and the fits' working data are freed here
    let scale = if opts.keep_sct {
        sct.rows = sctransform::matrix::GeneRows { ncol: 0, p: vec![0], j: Vec::new(), x: Vec::new() };
        sct.residuals.y = None;
        if !opts.correct {
            sct.corrected = None;
        }
        sct_result.kept = Some(sct);
        None
    } else {
        let s = std::mem::take(&mut sct.residuals.scale_data);
        drop(sct);
        Some(s)
    };
    if let Some(npcs) = params.pca_guard {
        pca_guard_in(n, npcs, opts.guard_context)?;
    }
    let t1 = Instant::now();
    let entry = rng.random_seed();
    let scale_ref: &[f64] = match (&scale, &sct_result.kept) {
        (Some(s), _) => s,
        (None, Some(k)) => &k.residuals.scale_data,
        (None, None) => unreachable!("scale.data is kept or held"),
    };
    // RunPCA.default calls set.seed(seed.use) itself
    let p = pca::run_pca(scale_ref, rows.len(), n, &params.pca, &mut rng, threads).map_err(ChainError::Pca)?;
    let pca_exit = rng.random_seed();
    obs.pca(&rows, scale_ref, &p, Some(&entry), &pca_exit)?;
    drop(scale);
    let pca_summary = pca::output::result_json(&p, rows.len());
    let pca_s = t1.elapsed().as_secs_f64();
    obs.phase("pca");
    let ncol = p.npcs;
    let t2 = Instant::now();
    let nb = find_neighbors(&p.embeddings, n, ncol, &params.neighbors, threads).map_err(ChainError::Neighbors)?;
    obs.neighbors(&nb, ncol)?;
    let nb_summary = neighbors::output::result_json(&nb, ncol);
    let nnz = (nb.nn.nnz(), nb.snn.nnz());
    let neighbors_s = t2.elapsed().as_secs_f64();
    obs.phase("neighbors");
    // FindClusters reads the SNN graph; the index and the lists go unless kept
    let (kept_nb, snn) = if opts.keep_reductions { (Some(nb), None) } else { (None, Some(nb.snn)) };
    let snn_ref = match (&kept_nb, &snn) {
        (Some(k), _) => &k.snn,
        (None, Some(s)) => s,
        (None, None) => unreachable!("the SNN graph is kept or held"),
    };
    let mut chain = from_snn(snn_ref, Some((&p.embeddings, ncol)), params, threads, opts, obs)?;
    drop(snn);
    chain.timings.sctransform = sct_s;
    chain.timings.pca = pca_s;
    chain.timings.neighbors = neighbors_s;
    let pca_kept = if opts.keep_reductions { Some((rows, p)) } else { None };
    Ok(Stages {
        sct: sct_result,
        pca: PcaResult { summary: pca_summary, rng_exit: pca_exit, kept: pca_kept },
        neighbors: NeighborsResult { summary: nb_summary, nnz, kept: kept_nb },
        chain,
    })
}

/// From SCTransform's `scale.data` (`scrust cluster --sct`): `scale` holds the variable features'
/// rows, row-major (`rows.len() x n`), `rows` their rows in the counts.
pub fn from_scale_data(scale: Vec<f64>, rows: Vec<usize>, n: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.check().map_err(ChainError::Params)?;
    from_scale_data_with(scale, rows, n, None, params, threads, opts, obs)
}

#[allow(clippy::too_many_arguments)]
fn from_scale_data_with(scale: Vec<f64>, rows: Vec<usize>, n: usize, rng: Option<(RRng, Vec<i32>)>, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    if let Some(npcs) = params.pca_guard {
        pca_guard_in(n, npcs, opts.guard_context)?;
    }
    let t0 = Instant::now();
    // RunPCA.default calls set.seed(seed.use) itself
    let (mut rng, entry) = match rng {
        Some((g, e)) => (g, Some(e)),
        None => (RRng::new(params.pca.seed_use.expect("checked")), None),
    };
    let nf = rows.len();
    let p = pca::run_pca(&scale, nf, n, &params.pca, &mut rng, threads).map_err(ChainError::Pca)?;
    obs.pca(&rows, &scale, &p, entry.as_deref(), &rng.random_seed())?;
    drop(scale);
    let pca_s = t0.elapsed().as_secs_f64();
    obs.phase("pca");
    let ncol = p.npcs;
    let emb = p.embeddings;
    let mut chain = from_embeddings_inner(emb, n, ncol, params, threads, opts, obs)?;
    chain.timings.pca = pca_s;
    Ok(chain)
}

/// From the PCA's embeddings (`scrust cluster --pca`): `n x ncol`, column-major.
pub fn from_embeddings(emb: Vec<f64>, n: usize, ncol: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.neighbors.check(ncol).map_err(|e| ChainError::Params(format!("FindNeighbors: {e}")))?;
    params.cluster.check().map_err(|e| ChainError::Params(format!("FindClusters: {e}")))?;
    from_embeddings_inner(emb, n, ncol, params, threads, opts, obs)
}

fn from_embeddings_inner(emb: Vec<f64>, n: usize, ncol: usize, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    let t0 = Instant::now();
    let nb = find_neighbors(&emb, n, ncol, &params.neighbors, threads).map_err(ChainError::Neighbors)?;
    obs.neighbors(&nb, ncol)?;
    let neighbors_s = t0.elapsed().as_secs_f64();
    obs.phase("neighbors");
    let snn = nb.snn;
    let mut chain = from_snn(&snn, Some((&emb, ncol)), params, threads, opts, obs)?;
    chain.timings.neighbors = neighbors_s;
    Ok(chain)
}

/// From the SNN graph (`scrust cluster --snn`), with the embeddings (`n x ncol`, column-major)
/// when the scores are to be computed.
pub fn from_snn(snn: &CscMatrix, emb: Option<(&[f64], usize)>, params: &ChainParams, threads: usize, opts: Options, obs: &mut dyn Observer) -> Result<Chain, ChainError> {
    params.cluster.check().map_err(|e| ChainError::Params(format!("FindClusters: {e}")))?;
    let t0 = Instant::now();
    let clusters = find_clusters(snn, &params.cluster, threads, opts.keep_partitions).map_err(ChainError::Cluster)?;
    obs.clusters(&clusters)?;
    let cluster_s = t0.elapsed().as_secs_f64();
    obs.phase("cluster");
    let n = clusters.n_cells;
    let t1 = Instant::now();
    let silhouette = match (params.stats.silhouette, emb) {
        (true, Some((e, ncol))) => {
            if e.len() != n * ncol {
                return Err(ChainError::Silhouette(format!("embeddings of {} values for {n} cells x {ncol} columns", e.len())));
            }
            let s = approx_silhouette(e, n, ncol, &params.neighbors.dims, &clusters.seurat_clusters).map_err(ChainError::Silhouette)?;
            obs.silhouette(&s)?;
            obs.phase("silhouette");
            Some(s)
        }
        _ => None,
    };
    let silhouette_s = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let purity = match (params.stats.purity, emb) {
        (true, Some((e, ncol))) => {
            if e.len() != n * ncol {
                return Err(ChainError::Purity(format!("embeddings of {} values for {n} cells x {ncol} columns", e.len())));
            }
            let st = &params.stats;
            let p = neighbor_purity(e, n, ncol, &params.neighbors.dims, &clusters.seurat_clusters, st.purity_k, st.purity_weighted, opts.keep_purity_lists, threads).map_err(ChainError::Purity)?;
            obs.purity(&p)?;
            obs.phase("purity");
            Some(p)
        }
        _ => None,
    };
    let purity_s = t2.elapsed().as_secs_f64();
    Ok(Chain {
        n_cells: n,
        clusters,
        silhouette,
        purity,
        n_components: emb.map(|(_, c)| c),
        timings: Timings { cluster: cluster_s, silhouette: silhouette_s, purity: purity_s, ..Timings::default() },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_guard_stops_at_50_cells_or_fewer() {
        assert!(pca_guard(51, 50).is_ok());
        let e = pca_guard(50, 50).unwrap_err().to_string();
        assert!(e.starts_with("[LOW_QUALITY_SAMPLE] only 50 cell(s) remain - too few to compute 50 principal components"), "{e}");
    }
}
