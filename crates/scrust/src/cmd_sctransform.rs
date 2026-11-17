//! `scrust sctransform`: SCTransform v2 as Seurat 5.5.1's `SCTransform()` runs it with its
//! defaults (docs/PLAN.md B3a, milestone M4), on a 10x counts directory.

use crate::cli::{self, Args, Phases};
use rcompat::RRng;
use scio::{tenx, Json};
use std::path::Path;

pub const USAGE: &str = "\
Usage: scrust sctransform --counts <10x dir> --out <dir> [--threads N] [--seed 1448145]
                          [--variable-features 3000] [--ncells 5000]
                          [--from-stages <dir> --context <ctx>]

SCTransform v2 as Seurat 5.5.1's SCTransform() runs it with its defaults: sctransform 0.4.3's
vst(vst.flavor = \"v2\") with glmGamPoi 1.22.0's fits, Pearson residuals, the variable features
and scale.data, computed as R computes them (bit for bit under the same OpenBLAS kernel).

  --counts <dir>             10x Matrix Market directory (matrix.mtx(.gz), features.tsv(.gz),
                             barcodes.tsv(.gz)): the genes x cells counts SCTransform reads
  --out <dir>                output directory (must not hold a previous run)
  --threads N                threads (default: SLURM_CPUS_PER_TASK, else 1); results do not
                             depend on it
  --seed N                   seed.use, set before vst (default 1448145, Seurat's)
  --variable-features N      variable.features.n (default 3000)
  --ncells N                 cells used to fit the per-gene models (default 5000)
  --from-stages <dir>        with --context: run on the SCT input of one DoubletFinder context
  --context <ctx>            (sweep_pN<pN> or final) of an instrumented R run, built from the
                             counts and the run's draws (S1): the real cells, then the doublets
                             X1, X2, ... (the benchmark's like-for-like input)

Outputs:
  cell_attr.tsv              barcode, umi, log_umi
  gene_attr.tsv              the genes vst keeps: id, name, detection_rate, gmean, amean,
                             variance, residual_mean, residual_variance
  model_pars.tsv             the fitting genes' parameters (theta after the 1e-7 floor,
                             (Intercept), log_umi), in sampling order
  model_pars_fit.tsv         the regularised parameters of every kept gene
  variable_features.tsv      id and name of the variable features, highest residual variance first
  scale_data.npy             scale.data: variable features x cells, float64, column-major (R's
                             layout; NumPy reads it with its shape)
  summary.json               parameters, sizes, min_variance, bandwidth, timings, peak memory
Doubles are written with 17 significant digits.
";

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse(argv)?;
    args.only(&["counts", "out", "threads", "seed", "variable-features", "ncells", "from-stages", "context"])?;
    let counts_dir = args.path("counts")?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let mut params = sctransform::SctParams::default();
    params.seed_use = Some(args.number("seed", 1448145i32)?);
    params.variable_features_n = Some(args.number("variable-features", 3000usize)?);
    params.ncells = args.number("ncells", 5000usize)?;
    params.check()?;
    if out.join("summary.json").exists() || out.join("gene_attr.tsv").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }

    let mut barcodes = tenx::read_barcodes(&counts_dir).map_err(|e| e.to_string())?;
    let features = tenx::read_features(&counts_dir).map_err(|e| e.to_string())?;
    let mut counts = scio::sparse::read_mtx(&counts_dir).map_err(|e| e.to_string())?;
    if counts.nrow != features.len() || counts.ncol != barcodes.len() {
        return Err(format!("{}: the matrix is {} x {}, with {} features and {} barcodes", counts_dir.display(), counts.nrow, counts.ncol, features.len(), barcodes.len()));
    }
    let context = match (args.get("from-stages"), args.get("context")) {
        (Some(dir), Some(ctx)) => {
            use doubletfinder::sct::{doublet_pairs, real_cells, sct_input};
            let stages = doubletfinder::stages::StageDir::open(Path::new(dir)).map_err(|e| e.to_string())?;
            let real = real_cells(&stages, ctx, counts.ncol).map_err(|e| e.to_string())?;
            let pairs = doublet_pairs(&stages, ctx).map_err(|e| e.to_string())?;
            counts = sct_input(&counts, &real, &pairs);
            barcodes = real.iter().map(|&c| barcodes[c].clone()).chain((1..=pairs.len()).map(|k| format!("X{k}"))).collect();
            Some((dir.to_string(), ctx.to_string()))
        }
        (None, None) => None,
        _ => return Err("--from-stages and --context go together".into()),
    };
    phases.mark("load");

    // SCTransform seeds R's generator (set.seed(seed.use)) before it draws anything
    let mut rng = RRng::new(params.seed_use.unwrap_or(1448145));
    let sct = sctransform::sctransform(&counts, &params, &mut rng, false, threads).map_err(|e| e.to_string())?;
    phases.mark("sctransform");

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let genes: Vec<sctransform::output::Gene> = features.iter().map(|f| sctransform::output::Gene { id: f.id.clone(), name: f.name.clone() }).collect();
    sctransform::output::write_outputs(&out, &barcodes, &genes, &sct).map_err(|e| e.to_string())?;
    // the counts' per-cell metadata (CreateSeuratObject's nCount_RNA and nFeature_RNA), which
    // `scrust pca` and `scrust neighbors` carry forward for `scrust cluster` (M11)
    if context.is_none() {
        cluster::metadata::RnaMetadata::of(&counts).write(&out.join(cluster::metadata::FILE), &barcodes).map_err(|e| e.to_string())?;
    }
    phases.mark("write");

    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let s = &sct.setup;
    let nr = sct.residuals.top_features.len();
    let summary = Json::obj(vec![
        ("tool", Json::Str("scrust sctransform".into())),
        ("reference", Json::Str(sctransform::REFERENCE.into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        ("parameters", sctransform::output::params_json(&params)),
        ("result", sctransform::output::result_json(&sct, counts.nrow)),
        (
            "run",
            Json::obj(vec![
                ("counts", Json::Str(counts_dir.display().to_string())),
                ("from_stages", context.as_ref().map_or(Json::Null, |c| Json::Str(c.0.clone()))),
                ("context", context.as_ref().map_or(Json::Null, |c| Json::Str(c.1.clone()))),
                ("outdir", Json::Str(out.display().to_string())),
                ("threads", Json::Int(threads as i64)),
                ("start_time", Json::Str(start_time)),
                ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
                ("elapsed_s", Json::Num(phases.elapsed())),
                ("phases", phases.to_json()),
                ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
                ("fingerprint", fingerprint),
            ]),
        ),
    ]);
    write(&out.join("summary.json"), &scio::json::to_file_text(&summary))?;
    eprintln!(
        "scrust sctransform: {} cells, {} of {} genes kept, {} variable features; {:.1} s, {} threads",
        counts.ncol,
        s.genes.len(),
        counts.nrow,
        nr,
        phases.elapsed(),
        threads
    );
    Ok(())
}
