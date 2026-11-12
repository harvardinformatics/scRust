//! The whole of DoubletFinder in Rust (milestone M6): the embeddings [`crate::run::run`] ranks
//! neighbours on, made from the counts alone, with R's random-number stream reproduced from
//! the seed (docs/PLAN.md B1 and B2).
//!
//! [`RustSource`] is called in R's order and draws exactly what R draws, from one generator:
//!
//! 1. `set.seed(seed)`, then `paramSweep`: above 10,000 cells, `sample(1:n, 10000)` (stage S1);
//! 2. for each pN in turn: two `sample(real.cells, n_doublets, replace = TRUE)` draws (S1), the
//!    real cells and their doublets as DoubletFinder builds them ([`crate::sct::sct_input`]),
//!    `SCTransform()` (S2: three `set.seed(1448145)`, the cell and gene sampling, two draws per
//!    fitting bin) and `RunPCA(npcs = 10)` (S3: `set.seed(42)`, irlba's start vector). The next
//!    pN's doublets continue from the state this PCA leaves;
//! 3. `set.seed(seed)`, then `doubletFinder`: the same on all cells.
//!
//! Stages S4 to S7 are [`crate::run`]'s. The products inside SCTransform and PCA are OpenBLAS
//! calls (rcompat's binding), so the results are R's bit for bit under R's OpenBLAS kernel.
//!
//! **Memory.** Nothing of a context outlives it but its embeddings (cells x 10). Within one,
//! the doublet-augmented counts are freed once SCTransform has kept the rows of its genes
//! ([`sctransform::start`]); those rows once the residuals and `scale.data` are made; and
//! `scale.data`, the one copy of the PCA input (`RunPCA` reads it in place), once the PCA is
//! done. The counts themselves are freed when the final run's input is built, the last time
//! they are read. No N x N distance matrix exists: neighbours are ranked by the [`knn`] crate.

use crate::pann::Embeddings;
use crate::record::{DumpLevel, Recorder};
use crate::run::EmbeddingSource;
use crate::sct::sct_input;
use crate::DfError;
use rcompat::{as_character, sample_int, RRng};
use scio::CscMatrix;
use sctransform::{Sct, SctParams};

/// R's `.Random.seed` as `digest::digest(paste(.Random.seed, collapse = ","), algo = "sha256",
/// serialize = FALSE)`: the `rng_end_sha256` field of the reference runs' `summary.json`.
pub fn rng_digest(state: &[i32]) -> String {
    let text = state.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
    scio::sha256::hex(text.as_bytes())
}

/// paramSweep's subsample size: `sample(1:n, 10000)` when there are more than 10,000 cells.
pub const SWEEP_CELLS: usize = 10_000;

/// Embeddings computed by Rust from the counts: doublets, SCTransform and PCA, drawing from
/// R's random-number stream as R does. See the module documentation.
pub struct RustSource<'a> {
    counts: CscMatrix,
    seed: i32,
    rng: RRng,
    sct: SctParams,
    pca: pca::PcaParams,
    threads: usize,
    recorder: Option<Recorder>,
    mark: Box<dyn FnMut(&str) + 'a>,
}

impl<'a> RustSource<'a> {
    /// DoubletFinder's settings: SCTransform with Seurat's defaults, `RunPCA(npcs = pcs)`, and
    /// the generator in `set.seed(seed)`'s state (the workflow calls `set.seed(12345)` before
    /// `paramSweep`).
    pub fn new(counts: CscMatrix, seed: i32, pcs: usize, threads: usize) -> RustSource<'a> {
        RustSource {
            counts,
            seed,
            rng: RRng::new(seed),
            sct: SctParams::default(),
            pca: pca::PcaParams { npcs: pcs, ..pca::PcaParams::default() },
            threads: threads.max(1),
            recorder: None,
            mark: Box::new(|_| {}),
        }
    }

    /// Records the stage dumps (S1 to S3 and the RNG checkpoints here; S4 to S7 when the
    /// recorder is finished with the run's output).
    pub fn with_recorder(mut self, recorder: Recorder) -> Self {
        self.recorder = Some(recorder);
        self
    }

    /// `mark(phase)` is called after each SCTransform and PCA (`<ctx>: SCTransform`, `<ctx>:
    /// RunPCA`), for timings and memory.
    pub fn with_marks(mut self, mark: impl FnMut(&str) + 'a) -> Self {
        self.mark = Box::new(mark);
        self
    }

    /// The generator's state (R's `.Random.seed`).
    pub fn rng_state(&self) -> Vec<i32> {
        self.rng.random_seed()
    }

    /// The recorder, once the run is over (to add S4-S7 and write).
    pub fn take_recorder(&mut self) -> Option<Recorder> {
        self.recorder.take()
    }

    fn checkpoint(&mut self, ctx: &str, stage: &str, label: &str) {
        if let Some(r) = self.recorder.as_mut() {
            r.rng(ctx, stage, label, self.rng.random_seed());
        }
    }

    /// `real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)`, then `real.cells2`: pairs
    /// of positions in the real cells.
    fn doublets(&mut self, ctx: &str, n_real: usize, n_doublets: usize) -> Result<Vec<(usize, usize)>, DfError> {
        self.checkpoint(ctx, "S1", "before_doublets");
        let draw = |rng: &mut RRng| sample_int(rng, n_real, n_doublets, true).map_err(|e| DfError::Input(format!("{ctx}: doublets: {e}")));
        let a = draw(&mut self.rng)?;
        let b = draw(&mut self.rng)?;
        let pairs: Vec<(usize, usize)> = a.into_iter().zip(b).collect();
        if let Some(r) = self.recorder.as_mut() {
            r.context(ctx);
            r.s1_doublets(ctx, &pairs);
        }
        // SCTransform.Seurat's entry: CreateSeuratObject draws nothing
        self.checkpoint(ctx, "S1", "after_doublets");
        Ok(pairs)
    }

    /// `SCTransform()` then `RunPCA(npcs)` on one context's input (real cells, then doublets);
    /// returns the cell embeddings. The input is consumed: it is freed as soon as SCTransform's
    /// setup has kept the rows it needs.
    fn embed(&mut self, ctx: &str, input: CscMatrix) -> Result<Embeddings, DfError> {
        let ncells = input.ncol;
        let started = sctransform::start(&input, &self.sct, &mut self.rng).map_err(|e| DfError::Input(format!("{ctx}: SCTransform: {e}")))?;
        drop(input);
        let keep_y = self.recorder.as_ref().is_some_and(|r| r.level() == DumpLevel::Full);
        let sct = sctransform::finish(started, &self.sct, &mut self.rng, keep_y, self.threads).map_err(|e| DfError::Input(format!("{ctx}: SCTransform: {e}")))?;
        let Sct { setup, rows, fits, mean_cell_sum, reg, residuals, rng_entry, rng_exit } = sct;
        drop(rows);
        (self.mark)(&format!("{ctx}: SCTransform"));
        if let Some(r) = self.recorder.as_mut() {
            r.sct(ctx, &setup, &fits, mean_cell_sum, &reg, &residuals, &rng_entry, &rng_exit)?;
        }
        drop((fits, reg));
        // RunPCA(npcs = length(PCs)) on scale.data's rows, VariableFeatures order, read in place
        let rows: Vec<usize> = residuals.top_features.iter().map(|&k| setup.genes[k]).collect();
        drop(setup);
        let entry = self.rng.random_seed();
        let p = pca::run_pca(&residuals.scale_data, rows.len(), residuals.ncells, &self.pca, &mut self.rng, self.threads)
            .map_err(|e| DfError::Input(format!("{ctx}: RunPCA: {e}")))?;
        let exit = self.rng.random_seed();
        if let Some(r) = self.recorder.as_mut() {
            r.pca(ctx, &rows, &residuals, &p, &entry, &exit)?;
        }
        drop(residuals);
        (self.mark)(&format!("{ctx}: RunPCA"));
        let pca::Pca { embeddings, npcs, .. } = p;
        Embeddings::new(embeddings, ncells, npcs)
    }
}

impl EmbeddingSource for RustSource<'_> {
    fn sweep_cells(&mut self, n: usize) -> Result<Vec<usize>, DfError> {
        if n != self.counts.ncol {
            return Err(DfError::Input(format!("{n} cells, the counts have {}", self.counts.ncol)));
        }
        self.checkpoint("sweep", "S1", "paramSweep_entry");
        if n <= SWEEP_CELLS {
            return Ok((0..n).collect());
        }
        // real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace = FALSE)]
        let cells = sample_int(&mut self.rng, n, SWEEP_CELLS, false).map_err(|e| DfError::Input(format!("sweep subsample: {e}")))?;
        if let Some(r) = self.recorder.as_mut() {
            r.s1_subsample(&cells);
        }
        Ok(cells)
    }

    fn sweep_embeddings(&mut self, _pn_index: usize, pn: f64, real_cells: &[usize], n_doublets: usize) -> Result<Embeddings, DfError> {
        let ctx = format!("sweep_pN{}", as_character(pn));
        let pairs = self.doublets(&ctx, real_cells.len(), n_doublets)?;
        let input = sct_input(&self.counts, real_cells, &pairs);
        drop(pairs);
        let emb = self.embed(&ctx, input)?;
        // parallel_paramSweep's exit: the distances, order() and pANN draw nothing
        self.checkpoint(&ctx, "S1", "sweep_iteration_exit");
        Ok(emb)
    }

    fn final_embeddings(&mut self, _pn: f64, n_doublets: usize) -> Result<Embeddings, DfError> {
        // the workflow's set.seed(12345) before doubletFinder()
        self.rng.set_seed(self.seed);
        let n = self.counts.ncol;
        let pairs = self.doublets("final", n, n_doublets)?;
        let all: Vec<usize> = (0..n).collect();
        let input = sct_input(&self.counts, &all, &pairs);
        // the last read of the counts: free them for the final run
        self.counts = CscMatrix::zeros(0, 0);
        drop((all, pairs));
        let emb = self.embed("final", input)?;
        self.checkpoint("final", "S7", "doubletFinder_exit");
        Ok(emb)
    }
}
