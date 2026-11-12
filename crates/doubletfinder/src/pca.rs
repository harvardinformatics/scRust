//! RunPCA inside DoubletFinder (stage S3, milestone M5): Seurat's `RunPCA(npcs = 10)` on the SCT
//! assay of each context, and the comparison of Rust's PCA with an instrumented R run's dumps.
//!
//! - **Teacher-forced:** the PCA input is R's own (`S2d_scale_data.<ctx>.npy`, the rows
//!   `RunPCA.default` receives, dumped by the full and standard runs); for the light runs, which
//!   did not dump it, Rust's S2d on R's S2c parameters stands in (M4 showed it identical to R's).
//!   The generator starts from R's state at `RunPCA`'s entry.
//! - **Chained:** Rust's whole SCTransform from the context's SCT input (the counts and R's S1
//!   draws) and R's state at SCTransform's entry, then Rust's PCA on its `scale.data` with the
//!   generator SCTransform left: what the full Rust pipeline (M6) runs.
//!
//! Compared bit for bit: irlba's `d`, the standard deviations, the embeddings and the loadings
//! (S3), the PCA features (S2d_pca_features: PrepDR's filter), the scalars `pca_n_features`,
//! `irlba_iter`, `irlba_mprod` and `total_variance`, and the RNG checkpoints at RunPCA's entry,
//! at irlba's entry (after `set.seed(42)`) and exit, and at RunPCA's exit. The products are
//! OpenBLAS calls, so the run must use the reference run's kernel.

use crate::sct::{block_f64, block_i32, cmp_f64, cmp_idx, cmp_rng, cmp_scalar, context, sct_contexts, Context, RngCheckpoints, Scalars, RS2a};
use crate::stages::StageDir;
use crate::validate::{Check, Report, Status};
use crate::DfError;
use rcompat::RRng;

/// DoubletFinder's `RunPCA(seu, npcs = length(PCs))` with `PCs = 1:10`: Seurat's other defaults.
pub fn params() -> pca::PcaParams {
    pca::PcaParams { npcs: 10, ..pca::PcaParams::default() }
}

/// R's PCA input of a context, as `pca::run_pca` takes it (row-major features x cells), with the
/// features' input rows (`S2d_pca_features`); `None` when the run did not dump it.
pub fn r_pca_input(stages: &StageDir, ctx: &str) -> Result<Option<(Vec<f64>, usize, usize, Vec<usize>)>, DfError> {
    let file = format!("S2d_scale_data.{ctx}.npy");
    let Some(b) = stages.block(&file, ctx).cloned() else { return Ok(None) };
    if !stages.dir.join(&file).is_file() {
        return Ok(None);
    }
    let a = stages.read(&file)?;
    let (nf, nc) = (b.nrow, b.ncol);
    let mut rows = vec![0.0; nf * nc];
    for c in 0..nc {
        for f in 0..nf {
            rows[f * nc + c] = a.get(b.row_start + f, c);
        }
    }
    let feats = block_i32(stages, "S2d_pca_features.npy", ctx)?.ok_or_else(|| DfError::Input(format!("{ctx}: no S2d_pca_features")))?;
    let genes = feats.0[0].iter().map(|&g| g as usize).collect();
    Ok(Some((rows, nf, nc, genes)))
}

fn missing(stage: &str, name: &str) -> Check {
    Check { stage: stage.into(), name: name.into(), status: Status::Missing, n: 0, n_differ: 0, detail: "not dumped".into() }
}

/// Rust's PCA of one context against R's dumps, under `stage`. `genes` are the input rows of the
/// PCA's input rows (before PrepDR's filter); `rng_entry` and `rng_exit` are the generator's
/// states when RunPCA starts and ends.
pub fn compare_s3(stages: &StageDir, rng: &RngCheckpoints, scalars: &Scalars, ctx: &str, stage: &str, genes: &[usize], pca: &pca::Pca, rng_entry: &[i32], rng_exit: &[i32]) -> Result<Vec<Check>, DfError> {
    let mut checks = Vec::new();
    checks.push(cmp_rng(stage, &format!("{ctx}: RNG at RunPCA's entry"), rng_entry, rng.state(ctx, "pca_entry")));
    checks.push(cmp_rng(stage, &format!("{ctx}: RNG at irlba's entry"), &pca.rng_irlba_entry, rng.state(ctx, "irlba_entry")));
    checks.push(cmp_rng(stage, &format!("{ctx}: RNG at irlba's exit"), rng_exit, rng.state(ctx, "irlba_exit")));
    checks.push(cmp_rng(stage, &format!("{ctx}: RNG at RunPCA's exit"), rng_exit, rng.state(ctx, "pca_exit")));
    let feats: Vec<usize> = pca.features.iter().map(|&k| genes[k]).collect();
    match block_i32(stages, "S2d_pca_features.npy", ctx)? {
        Some((cols, _)) => checks.push(cmp_idx(stage, &format!("{ctx}: PCA features"), &feats, &cols[0])),
        None => checks.push(missing(stage, &format!("{ctx}: S2d_pca_features"))),
    }
    let ir = &pca.irlba;
    for (name, v) in [("pca_n_features", feats.len() as f64), ("irlba_iter", ir.iter as f64), ("irlba_mprod", ir.mprod as f64), ("total_variance", pca.total_variance)] {
        checks.push(cmp_scalar(stage, &format!("{ctx}: {name}"), v, scalars.get(ctx, name)));
    }
    let npcs = pca.npcs;
    for (file, name, got, nrow) in [
        ("S3_irlba_d.npy", "irlba d", &ir.d, 1),
        ("S3_sdev.npy", "sdev", &pca.sdev, 1),
        ("S3_embeddings.npy", "embeddings", &pca.embeddings, pca.ncells),
        ("S3_loadings.npy", "loadings", &pca.loadings, feats.len()),
    ] {
        let label = format!("{ctx}: {name}");
        match block_f64(stages, file, ctx)? {
            Some((cols, rows)) if cols.len() == npcs && rows == nrow => {
                // R's block is nrow x npcs; Rust's arrays are column-major nrow x npcs too
                let want: Vec<f64> = cols.into_iter().flatten().collect();
                checks.push(cmp_f64(stage, &label, got, &want));
            }
            Some((cols, rows)) => checks.push(Check {
                stage: stage.into(),
                name: label,
                status: Status::Differ,
                n: rows * cols.len(),
                n_differ: rows * cols.len(),
                detail: format!("{nrow} x {npcs} against R's {rows} x {}", cols.len()),
            }),
            None => checks.push(missing(stage, &label)),
        }
    }
    Ok(checks)
}

/// The teacher-forced PCA of one context: R's PCA input (or, for the light runs, Rust's S2d on
/// R's S2c parameters) and R's RNG state at RunPCA's entry. Returns the checks and where the
/// input came from (`"R"` or `"Rust S2d"`).
pub fn check_teacher(stages: &StageDir, rng: &RngCheckpoints, scalars: &Scalars, c: &Context, sct_params: &sctransform::SctParams, threads: usize) -> Result<(Vec<Check>, &'static str), DfError> {
    let ctx = &c.ctx;
    let (rows, nf, nc, genes, source) = match r_pca_input(stages, ctx)? {
        Some((rows, nf, nc, genes)) => (rows, nf, nc, genes, "R"),
        None => {
            let r = RS2a::read(stages, ctx)?;
            let (res, _) = crate::sct::check_s2d(stages, scalars, c, &r, sct_params, threads)?;
            let genes: Vec<usize> = res.top_features.iter().map(|&k| r.genes[k]).collect();
            let (nf, nc) = (genes.len(), res.ncells);
            (res.scale_data, nf, nc, genes, "Rust S2d")
        }
    };
    let mut g = rng.rng(ctx, "pca_entry")?;
    let entry = g.random_seed();
    let p = pca::run_pca(&rows, nf, nc, &params(), &mut g, threads).map_err(|e| DfError::Input(format!("{ctx}: {e}")))?;
    let checks = compare_s3(stages, rng, scalars, ctx, "S3", &genes, &p, &entry, &g.random_seed())?;
    Ok((checks, source))
}

/// The chained run of one context: Rust's SCTransform from the SCT input and R's state at
/// SCTransform's entry, then Rust's PCA with the generator SCTransform left (stage "S3-chain").
pub fn check_chained(stages: &StageDir, rng: &RngCheckpoints, scalars: &Scalars, c: &Context, sct_params: &sctransform::SctParams, threads: usize) -> Result<Vec<Check>, DfError> {
    let ctx = &c.ctx;
    let mut g: RRng = c.rng_entry.clone();
    let sct = sctransform::sctransform(&c.input, sct_params, &mut g, false, threads).map_err(|e| DfError::Input(format!("{ctx}: {e}")))?;
    let res = &sct.residuals;
    let genes: Vec<usize> = res.top_features.iter().map(|&k| sct.setup.genes[k]).collect();
    let entry = g.random_seed();
    let p = pca::run_pca(&res.scale_data, genes.len(), res.ncells, &params(), &mut g, threads).map_err(|e| DfError::Input(format!("{ctx}: {e}")))?;
    compare_s3(stages, rng, scalars, ctx, "S3-chain", &genes, &p, &entry, &g.random_seed())
}

/// Which PCA checks to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Which {
    pub teacher: bool,
    pub chain: bool,
}

impl Which {
    pub const ALL: Which = Which { teacher: true, chain: true };
    pub const TEACHER: Which = Which { teacher: true, chain: false };
    pub const CHAIN: Which = Which { teacher: false, chain: true };
}

/// Stage S3 of an instrumented R run against Rust's PCA, in every context (or those listed), and
/// per context where the teacher-forced input came from.
pub fn validate(input: &std::path::Path, stages_dir: &std::path::Path, which: Which, only: Option<&[String]>, threads: usize) -> Result<(Report, Vec<(String, &'static str)>), DfError> {
    let counts = scio::sparse::read_mtx(input).map_err(DfError::Io)?;
    let stages = StageDir::open(stages_dir)?;
    let rng = RngCheckpoints::open(&stages)?;
    let scalars = Scalars::open(&stages)?;
    let sct_params = sctransform::SctParams::default();
    let mut checks = Vec::new();
    let mut sources = Vec::new();
    for ctx in sct_contexts(&stages) {
        if only.is_some_and(|o| !o.iter().any(|c| c == &ctx)) {
            continue;
        }
        let c = context(&stages, &rng, &counts, &ctx)?;
        if which.teacher {
            let (ch, source) = check_teacher(&stages, &rng, &scalars, &c, &sct_params, threads)?;
            checks.extend(ch);
            sources.push((ctx.clone(), source));
        }
        if which.chain {
            checks.extend(check_chained(&stages, &rng, &scalars, &c, &sct_params, threads)?);
        }
    }
    Ok((Report { reference: stages_dir.to_path_buf(), candidate: input.to_path_buf(), checks }, sources))
}
