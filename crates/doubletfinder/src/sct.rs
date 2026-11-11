//! SCTransform inside DoubletFinder (stage S2, milestone M4): the matrix DoubletFinder passes to
//! `SCTransform`, and the comparison of Rust's SCTransform with an instrumented R run's stage
//! dumps, sub-stage by sub-stage (teacher forcing, docs/PLAN.md A5b) and chained.
//!
//! **The SCT input.** `paramSweep` and `doubletFinder` build it as
//! `cbind(data, (data[, real.cells1] + data[, real.cells2]) / 2)`, where `data` is the counts
//! matrix restricted to the real cells (the 10,000-cell sweep subsample, or all cells) and the
//! doublet columns are named `X1..Xn`. Matrix 1.7-6 computes the sum of two `dgCMatrix` by
//! concatenating their stored entries as triplets and summing duplicates (`.Arith.Csparse`,
//! `tsort`), so a doublet column stores the union of its two parents' stored entries, explicit
//! zeros included; `/ 2` divides the stored values in place, and `cbind` and
//! `CreateSeuratObject` keep every stored entry (checked on R's own construction by
//! `validation/R/sctransform_vectors.R`).

use crate::stages::StageDir;
use crate::validate::{Check, Status};
use crate::DfError;
use rcompat::RRng;
use scio::{CscMatrix, Npy, Table};

/// DoubletFinder's SCT input: the real cells (columns `real_cells` of `counts`, in that order)
/// followed by one doublet per pair `(a, b)` of positions in the real cells.
pub fn sct_input(counts: &CscMatrix, real_cells: &[usize], pairs: &[(usize, usize)]) -> CscMatrix {
    let data = counts.select_cols(real_cells);
    let mut doublets = CscMatrix::zeros(data.nrow, pairs.len());
    doublets.p.clear();
    doublets.p.push(0);
    for &(a, b) in pairs {
        let (ra, va) = data.col(a);
        let (rb, vb) = data.col(b);
        let (mut i, mut k) = (0, 0);
        while i < ra.len() || k < rb.len() {
            let take_a = k >= rb.len() || (i < ra.len() && ra[i] < rb[k]);
            let take_b = i >= ra.len() || (k < rb.len() && rb[k] < ra[i]);
            if take_a {
                doublets.i.push(ra[i]);
                doublets.x.push(va[i] / 2.0);
                i += 1;
            } else if take_b {
                doublets.i.push(rb[k]);
                doublets.x.push(vb[k] / 2.0);
                k += 1;
            } else {
                doublets.i.push(ra[i]);
                doublets.x.push((va[i] + vb[k]) / 2.0);
                i += 1;
                k += 1;
            }
        }
        doublets.p.push(doublets.i.len());
    }
    data.cbind(&doublets)
}

/// `.Random.seed` checkpoints of an instrumented run (`rng.npy`, `rng.tsv`).
pub struct RngCheckpoints {
    states: Npy,
    labels: Table,
}

impl RngCheckpoints {
    pub fn open(stages: &StageDir) -> Result<RngCheckpoints, DfError> {
        let states = stages.read("rng.npy")?;
        let labels = Table::read(&stages.dir.join("rng.tsv")).map_err(DfError::Io)?;
        Ok(RngCheckpoints { states, labels })
    }

    /// The state labelled `label` in context `ctx` (the first one, if several).
    pub fn state(&self, ctx: &str, label: &str) -> Option<Vec<i32>> {
        let row: usize = self.labels.rows.iter().find(|r| r[1] == ctx && r[3] == label)?[0].parse().ok()?;
        Some((0..self.states.ncol()).map(|c| self.states.get_i32(row, c)).collect())
    }

    /// Every label of context `ctx`, in order.
    pub fn labels(&self, ctx: &str) -> Vec<String> {
        self.labels.rows.iter().filter(|r| r[1] == ctx).map(|r| r[3].clone()).collect()
    }

    pub fn rng(&self, ctx: &str, label: &str) -> Result<RRng, DfError> {
        let s = self.state(ctx, label).ok_or_else(|| DfError::Input(format!("no RNG checkpoint {label} in context {ctx}")))?;
        RRng::from_random_seed(&s).map_err(|e| DfError::Input(format!("RNG checkpoint {ctx}/{label}: {e}")))
    }
}

/// One context's block of a stacked array, column by column (`ncol` columns).
pub fn block_f64(stages: &StageDir, file: &str, ctx: &str) -> Result<Option<(Vec<Vec<f64>>, usize)>, DfError> {
    let Some(b) = stages.block(file, ctx).cloned() else { return Ok(None) };
    let a = stages.read(file)?;
    let cols = (0..b.ncol).map(|c| a.col_f64(c, b.row_start, b.nrow)).collect();
    Ok(Some((cols, b.nrow)))
}

pub fn block_i32(stages: &StageDir, file: &str, ctx: &str) -> Result<Option<(Vec<Vec<i32>>, usize)>, DfError> {
    let Some(b) = stages.block(file, ctx).cloned() else { return Ok(None) };
    let a = stages.read(file)?;
    let cols = (0..b.ncol).map(|c| a.col_i32(c, b.row_start, b.nrow)).collect();
    Ok(Some((cols, b.nrow)))
}

/// The contexts of a run that called SCTransform, in run order (`sweep_pN<pN>` ..., `final`).
pub fn sct_contexts(stages: &StageDir) -> Vec<String> {
    stages.blocks_of("S2a_cells.npy").iter().map(|b| b.ctx.clone()).collect()
}

/// The real cells of a context: the sweep subsample (or all cells) for sweep contexts, all cells
/// for the final run.
pub fn real_cells(stages: &StageDir, ctx: &str, n_cells: usize) -> Result<Vec<usize>, DfError> {
    if ctx != "final" {
        if let Some((c, _)) = block_i32(stages, "S1_subsample.npy", "sweep")? {
            return Ok(c[0].iter().map(|&i| i as usize).collect());
        }
    }
    Ok((0..n_cells).collect())
}

/// The doublet pairs (positions in the real cells) of a context, from `S1_doublets`.
pub fn doublet_pairs(stages: &StageDir, ctx: &str) -> Result<Vec<(usize, usize)>, DfError> {
    let (c, _) = block_i32(stages, "S1_doublets.npy", ctx)?.ok_or_else(|| DfError::Input(format!("no S1_doublets for {ctx}")))?;
    Ok(c[0].iter().zip(&c[1]).map(|(&a, &b)| (a as usize, b as usize)).collect())
}

// --- comparisons ---------------------------------------------------------------------------

fn bits_equal(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

/// Bitwise comparison of two double vectors.
pub fn cmp_f64(stage: &str, name: &str, got: &[f64], want: &[f64]) -> Check {
    if got.len() != want.len() {
        return Check { stage: stage.into(), name: name.into(), status: Status::Differ, n: want.len(), n_differ: want.len(), detail: format!("length {} against R's {}", got.len(), want.len()) };
    }
    let diff: Vec<usize> = (0..want.len()).filter(|&k| !bits_equal(got[k], want[k])).collect();
    let mut detail = String::new();
    if let Some(&k) = diff.first() {
        let max_rel = diff.iter().map(|&k| ((got[k] - want[k]) / want[k]).abs()).fold(0.0f64, |m, v| if v > m || v.is_nan() { v } else { m });
        detail = format!("first at {k}: {:e} against R's {:e}; max relative difference {max_rel:e}", got[k], want[k]);
    }
    Check { stage: stage.into(), name: name.into(), status: if diff.is_empty() { Status::Identical } else { Status::Differ }, n: want.len(), n_differ: diff.len(), detail }
}

/// Equality of two integer vectors.
pub fn cmp_idx(stage: &str, name: &str, got: &[usize], want: &[i32]) -> Check {
    if got.len() != want.len() {
        return Check { stage: stage.into(), name: name.into(), status: Status::Differ, n: want.len(), n_differ: want.len(), detail: format!("length {} against R's {}", got.len(), want.len()) };
    }
    let diff: Vec<usize> = (0..want.len()).filter(|&k| got[k] as i64 != want[k] as i64).collect();
    let detail = diff.first().map(|&k| format!("first at {k}: {} against R's {}", got[k], want[k])).unwrap_or_default();
    Check { stage: stage.into(), name: name.into(), status: if diff.is_empty() { Status::Identical } else { Status::Differ }, n: want.len(), n_differ: diff.len(), detail }
}

/// Equality of an RNG state with R's checkpoint.
pub fn cmp_rng(stage: &str, name: &str, got: &[i32], want: Option<Vec<i32>>) -> Check {
    match want {
        None => Check { stage: stage.into(), name: name.into(), status: Status::Missing, n: 0, n_differ: 0, detail: "no checkpoint".into() },
        Some(w) => {
            let same = got == w.as_slice();
            Check { stage: stage.into(), name: name.into(), status: if same { Status::Identical } else { Status::Differ }, n: 1, n_differ: (!same) as usize, detail: String::new() }
        }
    }
}

fn missing(stage: &str, name: &str) -> Check {
    Check { stage: stage.into(), name: name.into(), status: Status::Missing, n: 0, n_differ: 0, detail: "not dumped".into() }
}

/// The SCT input of one context and R's RNG state at `SCTransform`'s entry.
pub struct Context {
    pub ctx: String,
    pub input: CscMatrix,
    pub rng_entry: RRng,
}

/// Rebuilds a context's SCT input from the counts and R's S1 draws.
pub fn context(stages: &StageDir, rng: &RngCheckpoints, counts: &CscMatrix, ctx: &str) -> Result<Context, DfError> {
    let real = real_cells(stages, ctx, counts.ncol)?;
    let pairs = doublet_pairs(stages, ctx)?;
    let input = sct_input(counts, &real, &pairs);
    let rng_entry = rng.rng(ctx, "after_doublets")?;
    Ok(Context { ctx: ctx.to_string(), input, rng_entry })
}

/// `SCTransform`'s three `set.seed(seed.use)` calls before `vst` (the last one counts).
pub fn reseed(rng: &mut RRng, params: &sctransform::SctParams) {
    if let Some(s) = params.seed_use {
        rng.set_seed(s);
    }
}

/// S2a on one context, from R's state at SCTransform's entry; compared with every S2a array of
/// the context and the RNG checkpoints at vst's entry, after the cell sampling and at the fits'
/// entry.
pub fn check_s2a(stages: &StageDir, rng: &RngCheckpoints, c: &Context, params: &sctransform::SctParams) -> Result<(sctransform::Setup, sctransform::GeneRows, Vec<Check>), DfError> {
    let ctx = &c.ctx;
    let mut checks = Vec::new();
    let mut g = c.rng_entry.clone();
    reseed(&mut g, params);
    checks.push(cmp_rng("S2a", &format!("{ctx}: RNG at vst entry"), &g.random_seed(), rng.state(ctx, "vst_entry")));
    let (setup, rows) = sctransform::setup(&c.input, params, &mut g).map_err(|e| DfError::Input(format!("{ctx}: {e}")))?;
    checks.push(cmp_rng("S2a", &format!("{ctx}: RNG after cell sampling"), &setup.rng_after_cells, rng.state(ctx, "after_cell_sampling")));
    checks.push(cmp_rng("S2a", &format!("{ctx}: RNG at the fits' entry"), &setup.rng_after_genes, rng.state(ctx, "fits_entry")));
    match block_f64(stages, "S2a_cells.npy", ctx)? {
        Some((cols, _)) => {
            checks.push(cmp_f64("S2a", &format!("{ctx}: umi"), &setup.cell_attr.umi, &cols[0]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: log_umi"), &setup.cell_attr.log_umi, &cols[1]));
        }
        None => checks.push(missing("S2a", &format!("{ctx}: S2a_cells"))),
    }
    match block_i32(stages, "S2a_genes.npy", ctx)? {
        Some((cols, _)) => checks.push(cmp_idx("S2a", &format!("{ctx}: genes"), &setup.genes, &cols[0])),
        None => checks.push(missing("S2a", &format!("{ctx}: S2a_genes"))),
    }
    match block_f64(stages, "S2a_gene_stats.npy", ctx)? {
        Some((cols, _)) => {
            checks.push(cmp_f64("S2a", &format!("{ctx}: log_gmean"), &setup.log_gmean, &cols[0]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: amean"), &setup.amean, &cols[1]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: var"), &setup.var, &cols[2]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: detection_rate"), &setup.detection_rate(), &cols[3]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: gmean"), &setup.gmean(), &cols[4]));
        }
        None => checks.push(missing("S2a", &format!("{ctx}: S2a_gene_stats"))),
    }
    let idx_checks = [
        ("S2a_cells_step1.npy", "cells_step1", &setup.cells_step1),
        ("S2a_genes_step1_pool.npy", "genes_step1 pool", &setup.genes_step1_pool),
        ("S2a_genes_step1.npy", "genes_step1", &setup.genes_step1),
    ];
    for (file, name, got) in idx_checks {
        match block_i32(stages, file, ctx)? {
            Some((cols, _)) => checks.push(cmp_idx("S2a", &format!("{ctx}: {name}"), got, &cols[0])),
            None => checks.push(missing("S2a", &format!("{ctx}: {name}"))),
        }
    }
    match block_f64(stages, "S2a_log_gmean_step1.npy", ctx)? {
        Some((cols, _)) => checks.push(cmp_f64("S2a", &format!("{ctx}: log_gmean_step1"), &setup.log_gmean_step1, &cols[0])),
        None => checks.push(missing("S2a", &format!("{ctx}: log_gmean_step1"))),
    }
    let r_prob = block_f64(stages, "S2a_sampling_prob.npy", ctx)?;
    let r_dens = block_f64(stages, "S2a_density.npy", ctx)?;
    match (&setup.sampling, r_prob, r_dens) {
        (Some(s), Some((p, _)), Some((d, _))) => {
            checks.push(cmp_f64("S2a", &format!("{ctx}: density x"), &s.density_x, &d[0]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: density y"), &s.density_y, &d[1]));
            checks.push(cmp_f64("S2a", &format!("{ctx}: sampling_prob"), &s.sampling_prob, &p[0]));
        }
        (None, None, None) => {}
        (s, p, _) => checks.push(Check {
            stage: "S2a".into(),
            name: format!("{ctx}: gene sampling"),
            status: Status::Differ,
            n: 1,
            n_differ: 1,
            detail: format!("Rust sampled genes: {}; R: {}", s.is_some(), p.is_some()),
        }),
    }
    Ok((setup, rows, checks))
}

/// Which SCTransform sub-stages to check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Which {
    pub s2a: bool,
}

impl Which {
    pub const ALL: Which = Which { s2a: true };
}

/// Stage S2 of an instrumented R run against Rust's SCTransform, in every context (or the
/// contexts listed): the SCT input rebuilt from the counts in `input` (a 10x directory) and R's
/// S1 draws, each sub-stage teacher-forced on R's own inputs.
pub fn validate(input: &std::path::Path, stages_dir: &std::path::Path, params: &sctransform::SctParams, which: Which, only: Option<&[String]>) -> Result<crate::validate::Report, DfError> {
    let counts = scio::sparse::read_mtx(input).map_err(DfError::Io)?;
    let stages = StageDir::open(stages_dir)?;
    let rng = RngCheckpoints::open(&stages)?;
    let mut checks = Vec::new();
    for ctx in sct_contexts(&stages) {
        if only.is_some_and(|o| !o.iter().any(|c| c == &ctx)) {
            continue;
        }
        let c = context(&stages, &rng, &counts, &ctx)?;
        if which.s2a {
            let (_, _, ch) = check_s2a(&stages, &rng, &c, params)?;
            checks.extend(ch);
        }
    }
    Ok(crate::validate::Report { reference: stages_dir.to_path_buf(), candidate: input.to_path_buf(), checks })
}
