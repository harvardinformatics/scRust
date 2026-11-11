//! `scrust doubletfinder`: DoubletFinder as the workflow calls it (docs/PLAN.md B1 and B3).
//!
//! Milestone M3 runs the DoubletFinder core on R's own embeddings (teacher forcing): with
//! `--from-stages <dir>`, the sweep subsample and the PCA embeddings of every sweep iteration
//! and of the final run come from an instrumented R run's stage dumps, and everything after
//! them is Rust's: neighbours, pANN, the sweep summary and pK, nExp, the calls. The outputs have
//! the reference runs' layout (`cells.tsv`, `sweep_stats.tsv`, `bcmvn.tsv`, `summary.json`), so
//! runs compare file to file.

use crate::cli::{self, Args, Phases};
use doubletfinder::pann::NeighbourDumps;
use doubletfinder::run::{run, Input, Params, PkChoice};
use doubletfinder::stages::{stage_dumps, StageSource};
use scio::{tenx, Json, Table};
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: scrust doubletfinder --counts <10x dir> --clusters <tsv> --out <dir> --from-stages <dir>
                            [--cluster-levels <txt>] [--threads N] [--pcs 10] [--pn 0.25]
                            [--nexp-rate 0.15] [--pk auto|<value>] [--seed 12345]
                            [--dump-stages none|light|full]

DoubletFinder 2.0.6 as the workflow calls it. Milestone M3: teacher forcing only.

  --counts <dir>          10x Matrix Market directory (barcodes.tsv(.gz), features.tsv(.gz))
  --clusters <tsv>        cluster labels: a header line, then barcode and label per cell, in
                          the matrix's column order (the export's seurat_clusters.tsv)
  --cluster-levels <txt>  the factor levels in order, one per line (default:
                          seurat_clusters_levels.txt next to --clusters, else the order of
                          first appearance)
  --from-stages <dir>     an instrumented R run's stages/: the sweep subsample and the PCA
                          embeddings come from there (teacher forcing); required until
                          SCTransform and PCA are ported (M4-M6)
  --out <dir>             output directory (must not hold a previous run)
  --threads N             threads (default: SLURM_CPUS_PER_TASK, else 1); results do not
                          depend on it
  --dump-stages <level>   none (default); light: stage arrays S4 (digests), S5, S6, S7 in the
                          instrumented R runs' layout, in <out>/stages/; full: also the S4
                          order positions
";

/// Reads the cluster labels (second column, or the column named seurat_clusters) and checks
/// the barcodes against the matrix's.
fn read_clusters(path: &Path, barcodes: &[String]) -> Result<Vec<String>, String> {
    let t = Table::read(path).map_err(|e| e.to_string())?;
    let bc = t.col("barcode").unwrap_or(0);
    let cl = t.col("seurat_clusters").unwrap_or(1);
    if t.header.len() < 2 {
        return Err(format!("{}: expected a barcode and a cluster column", path.display()));
    }
    if t.rows.len() != barcodes.len() || t.rows.iter().zip(barcodes).any(|(r, b)| r.get(bc) != Some(b)) {
        return Err(format!("{}: barcodes differ from the matrix's (order or set)", path.display()));
    }
    Ok(t.rows.iter().map(|r| r.get(cl).cloned().unwrap_or_default()).collect())
}

fn read_levels(args: &Args, clusters_path: &Path, labels: &[String]) -> Result<(Vec<String>, String), String> {
    let explicit = args.get("cluster-levels").map(PathBuf::from);
    let beside = clusters_path.parent().map(|d| d.join("seurat_clusters_levels.txt"));
    let path = explicit.or_else(|| beside.filter(|p| p.is_file()));
    match path {
        Some(p) => {
            let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
            Ok((text.lines().map(str::to_owned).collect(), p.display().to_string()))
        }
        None => {
            let mut lv: Vec<String> = Vec::new();
            for l in labels {
                if !lv.contains(l) {
                    lv.push(l.clone());
                }
            }
            Ok((lv, "order of first appearance".into()))
        }
    }
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let args = Args::parse(argv)?;
    args.only(&["counts", "clusters", "cluster-levels", "from-stages", "out", "threads", "pcs", "pn", "nexp-rate", "pk", "seed", "dump-stages"])?;
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let counts = args.path("counts")?;
    let clusters_path = args.path("clusters")?;
    let out = args.path("out")?;
    let Some(stages) = args.get("from-stages").map(PathBuf::from) else {
        return Err("--from-stages is required: Rust's SCTransform and PCA come in milestones M4 to M6".into());
    };
    let threads = cli::threads(&args)?;
    let dump = args.get("dump-stages").unwrap_or("none");
    let neighbour_dumps = match dump {
        "none" => NeighbourDumps::default(),
        "light" => NeighbourDumps { digests: true, positions: false },
        "full" => NeighbourDumps { digests: true, positions: true },
        d => return Err(format!("--dump-stages: unknown level '{d}'")),
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
        keep_stages: dump != "none",
        neighbour_dumps,
    };
    if out.join("cells.tsv").exists() || out.join("summary.json").exists() {
        return Err(format!("{}: holds a previous run; choose another --out", out.display()));
    }

    let barcodes = tenx::read_barcodes(&counts).map_err(|e| e.to_string())?;
    let n_genes = tenx::read_features(&counts).map_err(|e| e.to_string())?.len();
    let clusters = read_clusters(&clusters_path, &barcodes)?;
    let (cluster_levels, levels_source) = read_levels(&args, &clusters_path, &clusters)?;
    let input = Input { barcodes, clusters, cluster_levels, n_genes };
    phases.mark("load");

    let mut source = StageSource::new(&stages).map_err(|e| e.to_string())?;
    let output = run(&input, &params, &mut source, &mut |p| phases.mark(p)).map_err(|e| e.to_string())?;

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    output.write_tables(&out).map_err(|e| e.to_string())?;
    if dump != "none" {
        stage_dumps(&output).write(&out.join("stages")).map_err(|e| e.to_string())?;
    }
    phases.mark("write");

    let fingerprint = Json::parse(&crate::fingerprint::json()).unwrap_or(Json::Null);
    let run_info = Json::obj(vec![
        ("kind", Json::Str("scrust doubletfinder".into())),
        ("mode", Json::Str("teacher-forced: R's sweep subsample and PCA embeddings (stages S1, S3)".into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        ("counts", Json::Str(counts.display().to_string())),
        ("clusters", Json::Str(clusters_path.display().to_string())),
        ("cluster_levels", Json::Str(levels_source)),
        ("from_stages", Json::Str(stages.display().to_string())),
        ("outdir", Json::Str(out.display().to_string())),
        ("threads", Json::Int(threads as i64)),
        ("dump_stages", Json::Str(dump.into())),
        ("start_time", Json::Str(start_time)),
        ("slurm_job_id", cli::env_json("SLURM_JOB_ID")),
        ("elapsed_s", Json::Num(phases.elapsed())),
        ("phases", phases.to_json()),
        ("peak_rss_kib", cli::peak_rss_kib().map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
        ("fingerprint", fingerprint),
    ]);
    let summary = Json::obj(vec![("result", output.result_json(None)), ("run", run_info)]);
    std::fs::write(out.join("summary.json"), scio::json::to_file_text(&summary)).map_err(|e| e.to_string())?;
    eprintln!(
        "scrust doubletfinder: pK {} (BCmetric row {}), nExp {}, {} doublets of {} cells; {:.1} s, {} threads",
        output.pk_label,
        output.sweep.as_ref().map_or("-".into(), |s| (s.best_row + 1).to_string()),
        output.nexp,
        output.doublet.iter().filter(|&&d| d).count(),
        output.n_cells,
        phases.elapsed(),
        threads
    );
    Ok(())
}
