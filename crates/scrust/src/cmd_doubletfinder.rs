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

use crate::cli::{self, Args, Phases};
use doubletfinder::chain::{rng_digest, RustSource};
use doubletfinder::input::{read_counts, read_input};
use doubletfinder::pann::NeighbourDumps;
use doubletfinder::record::{DumpLevel, Recorder};
use doubletfinder::run::{run, Params, PkChoice};
use doubletfinder::singlets::{write_singlets, DIR as SINGLETS};
use doubletfinder::stages::{stage_dumps, StageSource};
use scio::Json;
use std::cell::RefCell;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: scrust doubletfinder --counts <10x dir> --clusters <tsv> --out <dir>
                            [--seed 12345] [--pcs 10] [--pn 0.25] [--nexp-rate 0.15]
                            [--pk auto|<value>] [--threads N] [--cluster-levels <txt>]
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
                          the matrix's column order (the export's seurat_clusters.tsv)
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
                          embeddings from an instrumented R run's stages/ (dumps: S4-S7)

Outputs: cells.tsv (barcode, pANN and classification of every cell), sweep_stats.tsv (pN, pK,
BCreal), bcmvn.tsv (find.pK's table), summary.json (\"result\": pK and its label as R writes it,
nExp_poi, the homotypic proportion, nExp, the doublets made, k, the column names R writes, the
digest of the final .Random.seed; \"run\": parameters, versions, timings and peak memory per
phase), and singlets/: the counts without the doublets, as a new 10x directory (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz): every singlet's column exactly as stored, in input order,
and the input's feature lines. The counts are read again from --counts to write it, so that
they need not stay in memory through the final run.
";

pub fn main(argv: &[String]) -> Result<(), String> {
    let args = Args::parse(argv)?;
    args.only(&["counts", "clusters", "cluster-levels", "from-stages", "out", "threads", "pcs", "pn", "nexp-rate", "pk", "seed", "dump-stages"])?;
    let phases = RefCell::new(Phases::new());
    let start_time = cli::utc_now();
    let counts = args.path("counts")?;
    let clusters_path = args.path("clusters")?;
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
    if out.join("cells.tsv").exists() || out.join("summary.json").exists() || out.join(SINGLETS).exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }

    let (input, levels) = read_input(&counts, &clusters_path, args.get("cluster-levels").map(Path::new)).map_err(|e| e.to_string())?;
    let mut rng_end = None;
    let mut recorder = None;
    let output = match &stages {
        Some(dir) => {
            phases.borrow_mut().mark("load");
            let mut source = StageSource::new(dir).map_err(|e| e.to_string())?;
            run(&input, &params, &mut source, &mut |p| phases.borrow_mut().mark(p)).map_err(|e| e.to_string())?
        }
        None => {
            let m = read_counts(&counts, &input).map_err(|e| e.to_string())?;
            phases.borrow_mut().mark("load");
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

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    output.write_tables(&out).map_err(|e| e.to_string())?;
    match recorder {
        Some(r) => r.finish(&output).map_err(|e| format!("{}: {e}", out.join("stages").display()))?,
        None if dump.is_some() => stage_dumps(&output).write(&out.join("stages")).map_err(|e| e.to_string())?,
        None => {}
    }
    let mut phases = phases.into_inner();
    phases.mark("write");
    let software = format!("scrust {} doubletfinder", env!("CARGO_PKG_VERSION"));
    let singlets = write_singlets(&counts, &output, &out, &software, threads).map_err(|e| format!("{}: {e}", out.join(SINGLETS).display()))?;
    phases.mark("singlets");

    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let mode = match &stages {
        Some(_) => "teacher-forced: R's sweep subsample and PCA embeddings (stages S1, S3)",
        None => "rust: every stage (S1-S7) computed by Rust from the counts, the cluster labels and the seed",
    };
    let pk_text = match params.pk {
        PkChoice::Auto => Json::Str("auto".into()),
        PkChoice::Fixed(v) => Json::Num(v),
    };
    let run_info = Json::obj(vec![
        ("kind", Json::Str("scrust doubletfinder".into())),
        ("mode", Json::Str(mode.into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        (
            "reproduces",
            Json::strs(&[doubletfinder::REFERENCE, sctransform::REFERENCE, pca::REFERENCE]),
        ),
        ("counts", Json::Str(counts.display().to_string())),
        ("clusters", Json::Str(clusters_path.display().to_string())),
        ("cluster_levels", Json::Str(levels.describe())),
        ("from_stages", stages.as_ref().map(|d| Json::Str(d.display().to_string())).unwrap_or(Json::Null)),
        ("outdir", Json::Str(out.display().to_string())),
        (
            "parameters",
            Json::obj(vec![
                ("seed", Json::Int(params.seed as i64)),
                ("pcs", Json::Int(params.pcs as i64)),
                ("pN", Json::Num(params.pn)),
                ("nexp_rate", Json::Num(params.nexp_rate)),
                ("pK", pk_text),
            ]),
        ),
        ("threads", Json::Int(threads as i64)),
        ("dump_stages", Json::Str(dump_arg.into())),
        (
            "singlets",
            Json::obj(vec![
                ("dir", Json::Str(SINGLETS.into())),
                ("features", Json::Int(singlets.nrow as i64)),
                ("cells", Json::Int(singlets.ncol as i64)),
                ("entries", Json::Int(singlets.nnz as i64)),
                ("field", Json::Str(singlets.field.into())),
            ]),
        ),
        ("start_time", Json::Str(start_time)),
        ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
        ("elapsed_s", Json::Num(phases.elapsed())),
        ("phases", phases.to_json()),
        ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
        ("fingerprint", fingerprint),
    ]);
    let digest = rng_end.as_deref().map(rng_digest);
    let summary = Json::obj(vec![("result", output.result_json(digest)), ("run", run_info)]);
    std::fs::write(out.join("summary.json"), scio::json::to_file_text(&summary)).map_err(|e| e.to_string())?;
    eprintln!(
        "scrust doubletfinder: pK {} (BCmetric row {}), nExp {}, {} doublets of {} cells; {}/: {} singlets, {} entries; {:.1} s, {} threads, peak RSS {:.2} GiB",
        output.pk_label,
        output.sweep.as_ref().map_or("-".into(), |s| (s.best_row + 1).to_string()),
        output.nexp,
        output.doublet.iter().filter(|&&d| d).count(),
        output.n_cells,
        SINGLETS,
        singlets.ncol,
        singlets.nnz,
        phases.elapsed(),
        threads,
        cli::peak_rss_kib().unwrap_or(0) as f64 / 1048576.0
    );
    Ok(())
}
