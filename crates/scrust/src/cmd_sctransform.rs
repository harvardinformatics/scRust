//! `scrust sctransform`: SCTransform v2 as Seurat 5.5.1's `SCTransform()` runs it with its
//! defaults (docs/PLAN.md B3a, milestone M4), on a 10x counts directory; with
//! `--corrected-counts`, also the SCT assay's corrected counts and their `log1p` (C3f, M13): the
//! SCT part of `scrust cluster --export`.

use crate::cli::{self, Args, Phases};
use rcompat::RRng;
use scio::manifest::{Entry, Ids, Input, Manifest, Provenance};
use scio::{tenx, Json};
use std::path::Path;

/// The manifest's entry for `cell_metadata.tsv` (the counts' per-cell columns, M11), at `path`.
pub fn rna_entry(path: &str, ncells: usize) -> Entry {
    Entry::new(path, "rna_cell_metadata", "none", "tsv", "the counts' per-cell columns, as CreateSeuratObject computes them: barcode, nCount_RNA, nFeature_RNA")
        .header(true)
        .shape(&[ncells, 3])
        .rows(Ids::named("cells", path, "barcode"))
}

/// What README.txt says first about the SCT part (`scrust sctransform --corrected-counts`, and
/// `sct/` of `scrust cluster --export`).
pub fn readme_intro() -> Vec<String> {
    vec![
        "These files hold SCTransform's model and normalised data as Seurat 5.5.1's SCTransform() makes them (sctransform 0.4.3, glmGamPoi 1.22.0), computed by scrust. NONE OF THE MATRICES HERE ARE RAW COUNTS. The corrected counts are integers, corrected to the median sequencing depth by sctransform's correct(): read as raw counts, they would be taken for counts without complaint. For that reason no file here has a name a 10x reader looks for (matrix.mtx, barcodes.tsv, features.tsv, genes.tsv), and each matrix says what it holds in its comment lines.".into(),
        "summary.json lists every file (its files section): its role, its scale (raw_counts, corrected_counts, log1p_corrected_counts, residuals, ...), its shape, where its row and column names are, its size and SHA-256, and how the output was made. The schema is described in docs/output_manifest.md of the scRust repository.".into(),
        "scRust's loaders (loaders/ in the scRust repository) build objects from these files, each layer in a slot meant for its scale, through the manifest: in R, scrustio::load_seurat(\"<this directory>\") (the SCT assay with its SCTModel) and scrustio::load_sce(); in Python, scrustio.load_anndata(\"<this directory>\"). Each file below also loads with standard readers, as shown.".into(),
    ]
}

pub const USAGE: &str = "\
Usage: scrust sctransform --counts <10x dir> --out <dir> [--threads N] [--seed 1448145]
                          [--variable-features 3000] [--ncells 5000] [--corrected-counts]
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
  --corrected-counts         also the corrected counts of sctransform's correct() (Seurat's SCT
                             counts layer, do.correct.umi = TRUE) and their log1p (its data
                             layer): the SCT part of scrust cluster --export, never in 10x layout
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
  cell_metadata.tsv          barcode, nCount_RNA, nFeature_RNA (CreateSeuratObject's columns)
  summary.json               parameters, sizes, min_variance, bandwidth, timings, peak memory,
                             and the manifest (files): every file's role, scale, shape and
                             checksum, and how the output was made (docs/output_manifest.md)
With --corrected-counts, also (none of them raw counts, none under a 10x name):
  corrected_counts.mtx.gz        the corrected counts (integer Matrix Market): model genes x cells
  log1p_corrected_counts.mtx.gz  their log1p (real Matrix Market)
  sct_features.tsv.gz            the model genes (both matrices' rows), in the counts' order
  sct_cells.tsv.gz               the cells (their columns)
  README.txt                     what each file holds and how to load it
Doubles are written with 17 significant digits.
";

/// Options that change no value of any file the command writes: left out of the command the
/// normalised matrices' comment lines give, so that the files do not depend on them.
const OUTPUT_ONLY: [&str; 2] = ["out", "threads"];

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse_with_flags(argv, &["corrected-counts"])?;
    args.only(&["counts", "out", "threads", "seed", "variable-features", "ncells", "from-stages", "context", "corrected-counts"])?;
    let correct = args.flag("corrected-counts");
    if correct && args.get("from-stages").is_some() {
        return Err("--corrected-counts: from the counts only, not with --from-stages".into());
    }
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

    // the input's checksums: the manifest's provenance and the corrected counts' comment lines
    let input_files = tenx::input_files(&counts_dir).map_err(|e| e.to_string())?;
    let mut inputs = vec![Input::digest("counts", &counts_dir, &input_files.iter().map(|p| p.as_path()).collect::<Vec<_>>()).map_err(|e| e.to_string())?];
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
            // the stage dumps' draws: their directory (a validation input; files not digested)
            inputs.push(Input { option: "from-stages".into(), path: dir.to_string(), files: Vec::new() });
            Some((dir.to_string(), ctx.to_string()))
        }
        (None, None) => None,
        _ => return Err("--from-stages and --context go together".into()),
    };
    phases.mark("load");

    // SCTransform seeds R's generator (set.seed(seed.use)) before it draws anything
    let mut rng = RRng::new(params.seed_use.unwrap_or(1448145));
    let started = sctransform::start(&counts, &params, &mut rng).map_err(|e| e.to_string())?;
    let sct = sctransform::finish_with(started, &params, &mut rng, false, correct, threads).map_err(|e| e.to_string())?;
    phases.mark("sctransform");

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let genes: Vec<sctransform::output::Gene> = features.iter().map(|f| sctransform::output::Gene { id: f.id.clone(), name: f.name.clone() }).collect();
    sctransform::output::write_outputs(&out, &barcodes, &genes, &sct).map_err(|e| e.to_string())?;
    // the counts' per-cell metadata (CreateSeuratObject's nCount_RNA and nFeature_RNA), which
    // `scrust pca` and `scrust neighbors` carry forward for `scrust cluster` (M11)
    if context.is_none() {
        cluster::metadata::RnaMetadata::of(&counts).write(&out.join(cluster::metadata::FILE), &barcodes).map_err(|e| e.to_string())?;
    }
    let command = scio::manifest::canonical_command("sctransform", argv, &OUTPUT_ONLY, &["corrected-counts"]);
    if let Some(c) = &sct.corrected {
        let names: Vec<String> = sct.setup.genes.iter().map(|&g| features[g].name.clone()).collect();
        let (cc, dc) = sctransform::output::corrected_comments(&command, &cli::input_checksums(&inputs), c);
        c.write(&out, &names, &barcodes, &cc, &dc, threads).map_err(|e| e.to_string())?;
    }
    phases.mark("write");

    // the manifest (C3f): every file, its role and scale, and how the output was made
    let s = &sct.setup;
    let mut m = Manifest::new(
        "sctransform",
        Provenance {
            scrust_version: env!("CARGO_PKG_VERSION").into(),
            command_line: cli::command_line("sctransform", argv),
            command,
            inputs,
            reproduces: vec![sctransform::REFERENCE.into()],
            parameters: sctransform::output::params_json(&params),
            seeds: Json::obj(vec![("SCTransform seed.use", params.seed_use.map(|v| Json::Int(v as i64)).unwrap_or(Json::Null))]),
        },
    );
    m.push(Manifest::summary_entry());
    m.entries.extend(sctransform::output::manifest_entries("", counts.ncol, s.genes.len(), s.genes_step1.len(), sct.residuals.top_features.len(), sct.corrected.as_ref()));
    if context.is_none() {
        m.push(rna_entry(cluster::metadata::FILE, counts.ncol));
    }
    if sct.corrected.is_some() {
        m.push(Entry::new("README.txt", "readme", "none", "text", "this file: what each file holds and how to load it"));
        let text = scio::manifest::readme(&m, "scrust sctransform --corrected-counts: SCTransform's model and its normalised data", &readme_intro());
        write(&out.join("README.txt"), &text)?;
    }
    let files = cli::finish_manifest(&out, &mut m, threads)?;

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
        ("files", files),
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
