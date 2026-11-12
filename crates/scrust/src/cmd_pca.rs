//! `scrust pca`: Seurat 5.5.1's `RunPCA` (irlba 2.4.1) on the `scale.data` of an SCTransform run
//! (docs/PLAN.md A1 and B3b, milestone M5).

use crate::cli::{self, Args, Phases};
use rcompat::RRng;
use scio::{fmt17, Json, Table};
use std::fmt::Write as _;
use std::path::Path;

pub const USAGE: &str = "\
Usage: scrust pca --sct <dir> --out <dir> [--npcs 50] [--seed 42] [--threads N]

Seurat 5.5.1's RunPCA() with its defaults on an SCT assay, as R computes it (bit for bit under
the same OpenBLAS kernel): PrepDR's zero-variance filter, set.seed(seed.use), irlba 2.4.1's
truncated SVD of t(scale.data), the standard deviations, the loadings and the embeddings
u %*% diag(d). DoubletFinder runs it with --npcs 10.

  --sct <dir>        an SCTransform output directory in `scrust sctransform`'s layout:
                     scale_data.npy (variable features x cells), variable_features.tsv (the
                     rows' id and name, in VariableFeatures order), cell_attr.tsv (barcodes)
  --out <dir>        output directory (must not hold a previous run)
  --npcs N           number of components (default 50, Seurat's; at most features - 1)
  --seed N           seed.use, set before irlba draws its start vector (default 42)
  --threads N        threads (default: SLURM_CPUS_PER_TASK, else 1); results do not depend
                     on it

Outputs:
  embeddings.tsv     barcode, PC_1 .. PC_<npcs>: the cell embeddings
  loadings.tsv       id, name, PC_1 .. PC_<npcs>: the feature loadings (features kept)
  stdev.tsv          pc, stdev
  random_seed.txt    R's .Random.seed after RunPCA, one integer per line
  summary.json       parameters, sizes, irlba's d, iterations and products, total.variance,
                     timings, peak memory
Doubles are written with 17 significant digits.
";

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse(argv)?;
    args.only(&["sct", "out", "npcs", "seed", "threads"])?;
    let sct = args.path("sct")?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let params = pca::PcaParams { npcs: args.number("npcs", 50usize)?, seed_use: Some(args.number("seed", 42i32)?), ..pca::PcaParams::default() };
    if out.join("summary.json").exists() || out.join("embeddings.tsv").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }

    let cells = Table::read(&sct.join("cell_attr.tsv")).map_err(|e| format!("{}: {e}", sct.join("cell_attr.tsv").display()))?;
    let barcodes: Vec<String> = cells.rows.iter().map(|r| r[0].clone()).collect();
    let feats = Table::read(&sct.join("variable_features.tsv")).map_err(|e| format!("{}: {e}", sct.join("variable_features.tsv").display()))?;
    let (rows, nf, nc) = scio::npy::read_f64_row_major(&sct.join("scale_data.npy")).map_err(|e| e.to_string())?;
    if nf != feats.rows.len() || nc != barcodes.len() {
        return Err(format!("{}: scale_data.npy is {nf} x {nc}, with {} variable features and {} cells", sct.display(), feats.rows.len(), barcodes.len()));
    }
    phases.mark("load");

    // RunPCA.default calls set.seed(seed.use) itself, so the generator's earlier state is moot
    let mut rng = RRng::new(params.seed_use.unwrap_or(42));
    let prod_pieces = {
        let p = pca::Products::new(&rows, nc, nf, threads);
        let (r, c) = p.pieces();
        (r.len(), c.len())
    };
    let p = pca::run_pca(&rows, nf, nc, &params, &mut rng, threads).map_err(|e| e.to_string())?;
    drop(rows);
    phases.mark("pca");

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let pcs: Vec<String> = (1..=p.npcs).map(|k| format!("PC_{k}")).collect();
    let mut t = format!("barcode\t{}\n", pcs.join("\t"));
    for (c, b) in barcodes.iter().enumerate() {
        t.push_str(b);
        for k in 0..p.npcs {
            t.push('\t');
            t.push_str(&fmt17(p.embeddings[c + k * nc]));
        }
        t.push('\n');
    }
    write(&out.join("embeddings.tsv"), &t)?;
    let nk = p.features.len();
    let mut t = format!("id\tname\t{}\n", pcs.join("\t"));
    for (i, &f) in p.features.iter().enumerate() {
        let _ = write!(t, "{}\t{}", feats.rows[f][0], feats.rows[f][1]);
        for k in 0..p.npcs {
            t.push('\t');
            t.push_str(&fmt17(p.loadings[i + k * nk]));
        }
        t.push('\n');
    }
    write(&out.join("loadings.tsv"), &t)?;
    let mut t = String::from("pc\tstdev\n");
    for (k, s) in p.sdev.iter().enumerate() {
        let _ = writeln!(t, "{}\t{}", pcs[k], fmt17(*s));
    }
    write(&out.join("stdev.tsv"), &t)?;
    let seed: String = rng.random_seed().iter().map(|v| format!("{v}\n")).collect();
    write(&out.join("random_seed.txt"), &seed)?;
    phases.mark("write");

    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let ir = &p.irlba;
    let summary = Json::obj(vec![
        ("tool", Json::Str("scrust pca".into())),
        ("reference", Json::Str(pca::REFERENCE.into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        (
            "parameters",
            Json::obj(vec![
                ("npcs", Json::Int(params.npcs as i64)),
                ("seed_use", Json::Int(params.seed_use.unwrap_or(0) as i64)),
                ("weight_by_var", Json::Bool(params.weight_by_var)),
                ("irlba_maxit", Json::Int(params.irlba.maxit as i64)),
                ("irlba_tol", Json::Num(params.irlba.tol)),
            ]),
        ),
        (
            "result",
            Json::obj(vec![
                ("n_cells", Json::Int(nc as i64)),
                ("n_features", Json::Int(nf as i64)),
                ("n_pca_features", Json::Int(nk as i64)),
                ("npcs", Json::Int(p.npcs as i64)),
                ("irlba_d", Json::nums(&ir.d)),
                ("irlba_iter", Json::Int(ir.iter as i64)),
                ("irlba_mprod", Json::Int(ir.mprod as i64)),
                ("irlba_converged", Json::Bool(ir.converged)),
                ("irlba_work", Json::Int(ir.work as i64)),
                ("irlba_invariant_draws", Json::Int(ir.invariant_draws as i64)),
                ("total_variance", Json::Num(p.total_variance)),
            ]),
        ),
        (
            "run",
            Json::obj(vec![
                ("sct", Json::Str(sct.display().to_string())),
                ("outdir", Json::Str(out.display().to_string())),
                ("threads", Json::Int(threads as i64)),
                ("product_pieces", Json::obj(vec![("a_v_rows", Json::Int(prod_pieces.0 as i64)), ("t_w_a_cols", Json::Int(prod_pieces.1 as i64))])),
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
        "scrust pca: {nc} cells, {nk} of {nf} features, {} components; irlba {} iterations, {} products; {:.1} s, {threads} threads",
        p.npcs,
        ir.iter,
        ir.mprod,
        phases.elapsed()
    );
    Ok(())
}
