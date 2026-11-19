//! `scrust sctransform`: SCTransform v2 as Seurat 5.5.1's `SCTransform()` runs it with its
//! defaults (docs/PLAN.md B3a, milestone M4), on a 10x counts directory; with
//! `--corrected-counts`, also the SCT assay's corrected counts and their `log1p` (C3f, M13): the
//! SCT part of `scrust cluster --export`. With `--vars-to-regress` (M15b), SCTransform's second
//! regression of `scale.data`, and with it, `--mt-genes` or `--percent-mt`, the workflow's
//! `percent.mt` in `cell_metadata.tsv`.

use crate::cli::{self, Args, Phases};
use cluster::chain_params::ChainParams;
use rcompat::RRng;
use scio::manifest::{Entry, Ids, Input, Manifest, Provenance};
use scio::{tenx, Json};
use std::path::Path;

/// The manifest's entry for `cell_metadata.tsv` (the counts' per-cell columns, M11), at `path`,
/// with `ncols` columns: 3, or 4 with `percent.mt` (M15b).
pub fn rna_entry(path: &str, ncells: usize, ncols: usize) -> Entry {
    let what = if ncols == 4 {
        "the counts' per-cell columns, as CreateSeuratObject computes them, and percent.mt (PercentageFeatureSet of the mitochondrial features, as the workflow's scripts compute it before SCTransform; summary.json's result.percent_mt says which): barcode, nCount_RNA, nFeature_RNA, percent.mt"
    } else {
        "the counts' per-cell columns, as CreateSeuratObject computes them: barcode, nCount_RNA, nFeature_RNA"
    };
    Entry::new(path, "rna_cell_metadata", "none", "tsv", what).header(true).shape(&[ncells, ncols]).rows(Ids::named("cells", path, "barcode"))
}

/// [`rna_entry`] for a `cell_metadata.tsv` carried forward into `dir` (`scrust pca`, `scrust
/// neighbors`), its columns read from its header.
pub fn rna_entry_of(dir: &Path, path: &str, ncells: usize) -> Entry {
    let ncols = std::fs::read_to_string(dir.join(path)).ok().and_then(|t| t.lines().next().map(|h| h.split('\t').count())).unwrap_or(3);
    rna_entry(path, ncells, ncols)
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
                          [--vars-to-regress <name>[,<name>...]] [--percent-mt]
                          [--mt-genes <file>] [--metadata <tsv>]
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
  --vars-to-regress <names>  vars.to.regress (default: none): the per-cell covariates the second
                             regression removes from scale.data before it is centred (ScaleData's
                             RegressOutMatrix, lm's LINPACK QR through OpenBLAS): percent.mt,
                             nCount_RNA, nFeature_RNA, or numeric columns of --metadata (read as
                             R reads numbers), comma-separated, in the formula's order
  --percent-mt               compute percent.mt (PercentageFeatureSet) without regressing it:
                             cell_metadata.tsv's 4th column (implied by --mt-genes and by
                             regressing it); scrust pca, neighbors and cluster carry it
  --mt-genes <file>          percent.mt's mitochondrial features by id (column 1 of features.tsv):
                             one id per line, blank lines and lines starting with # skipped; ids
                             not found are reported and skipped, none found is an error. Without
                             it: the workflow's \"(?i)^mt-\" on Seurat's row names (column 2 made
                             unique, _ and | replaced by -). Regressing percent.mt when no feature
                             matches is an error that names --mt-genes
  --metadata <tsv>           a table with a barcode column (else the first), whose numeric columns
                             --vars-to-regress may name
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
  cell_metadata.tsv          barcode, nCount_RNA, nFeature_RNA (CreateSeuratObject's columns), and
                             percent.mt when asked for
  summary.json               parameters, sizes, min_variance, bandwidth (with the regression and
                             percent.mt's definition, features and ids not found, when asked for),
                             timings, peak memory,
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

/// The parameters as `summary.json` and the manifest record them: SCTransform's (with
/// `vars_to_regress` when set), and `percent_mt` and `mt_genes` when given (M15b), so that a run
/// without them records what earlier runs did.
fn params_json(cp: &ChainParams) -> Json {
    let mut j = sctransform::output::params_json(&cp.sct);
    if let Json::Object(m) = &mut j {
        if cp.percent_mt {
            m.push(("percent_mt".into(), Json::Bool(true)));
        }
        if let Some(f) = &cp.mt_genes {
            m.push(("mt_genes".into(), Json::Str(f.clone())));
        }
    }
    j
}

fn write(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn main(argv: &[String]) -> Result<(), String> {
    let mut phases = Phases::new();
    let start_time = cli::utc_now();
    let args = Args::parse_with_flags(argv, &["corrected-counts", "percent-mt"])?;
    args.only(&["counts", "out", "threads", "seed", "variable-features", "ncells", "from-stages", "context", "corrected-counts", "vars-to-regress", "percent-mt", "mt-genes", "metadata"])?;
    let correct = args.flag("corrected-counts");
    if correct && args.get("from-stages").is_some() {
        return Err("--corrected-counts: from the counts only, not with --from-stages".into());
    }
    // vars.to.regress and percent.mt (M15b), through the chain's parameters and rules
    let mut flags: Vec<(String, String)> = Vec::new();
    for name in ["vars-to-regress", "mt-genes"] {
        if let Some(v) = args.get(name) {
            flags.push((name.to_string(), v.to_string()));
        }
    }
    if args.flag("percent-mt") {
        flags.push(("percent-mt".to_string(), String::new()));
    }
    let counts_dir = args.path("counts")?;
    let out = args.path("out")?;
    let threads = cli::threads(&args)?;
    let mut params = sctransform::SctParams::default();
    params.seed_use = Some(args.number("seed", 1448145i32)?);
    params.variable_features_n = Some(args.number("variable-features", 3000usize)?);
    params.ncells = args.number("ncells", 5000usize)?;
    params.check()?;
    let cp = {
        let mut c = cluster::chain_params::build(None, None, &flags)?;
        c.sct = sctransform::SctParams { vars_to_regress: c.sct.vars_to_regress.clone(), ..params.clone() };
        c
    };
    params = cp.sct.clone();
    if cp.needs_covariates() && args.get("from-stages").is_some() {
        return Err("--vars-to-regress, --mt-genes and --percent-mt: from the counts only, not with --from-stages".into());
    }
    if args.get("metadata").is_some() && !crate::cmd_cluster::metadata_regressed(&cp) {
        return Err("--metadata: scrust sctransform reads it only for the columns --vars-to-regress names (scrust cluster carries a table's columns into metadata.tsv)".into());
    }
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

    // percent.mt and the covariates (M15b): computed from the counts before SCTransform, as the
    // workflow's scripts do, and handed to it in memory
    let mut rna = cluster::metadata::RnaMetadata::of(&counts);
    let mut percent_mt = None;
    let mut regression = None;
    if cp.needs_covariates() {
        let (cov, more) = crate::cmd_cluster::read_covariates(&cp, args.get("metadata"), &barcodes)?;
        inputs.extend(more);
        let lines = tenx::read_lines(&tenx::find(&counts_dir, "features.tsv").map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let c = scio::Counts::new(counts, barcodes.clone(), lines).map_err(|e| e.to_string())?;
        let (pm, reg, notes) = cluster::chain::covariates_for(&c, &cp, &cov, &mut rna).map_err(|e| e.to_string())?;
        for n in &notes {
            eprintln!("scrust sctransform: {n}");
        }
        counts = c.matrix;
        percent_mt = pm.map(|p| (p, notes));
        regression = reg;
    }
    // SCTransform seeds R's generator (set.seed(seed.use)) before it draws anything
    let mut rng = RRng::new(params.seed_use.unwrap_or(1448145));
    let started = sctransform::start(&counts, &params, &mut rng).map_err(|e| e.to_string())?;
    let sct = sctransform::finish_full(started, &params, &mut rng, false, correct, regression, threads).map_err(|e| e.to_string())?;
    phases.mark("sctransform");

    std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    let genes: Vec<sctransform::output::Gene> = features.iter().map(|f| sctransform::output::Gene { id: f.id.clone(), name: f.name.clone() }).collect();
    sctransform::output::write_outputs(&out, &barcodes, &genes, &sct).map_err(|e| e.to_string())?;
    // the counts' per-cell metadata (CreateSeuratObject's nCount_RNA and nFeature_RNA), which
    // `scrust pca` and `scrust neighbors` carry forward for `scrust cluster` (M11)
    if context.is_none() {
        rna.write(&out.join(cluster::metadata::FILE), &barcodes).map_err(|e| e.to_string())?;
    }
    let command = scio::manifest::canonical_command("sctransform", argv, &OUTPUT_ONLY, &["corrected-counts", "percent-mt"]);
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
            parameters: params_json(&cp),
            seeds: Json::obj(vec![("SCTransform seed.use", params.seed_use.map(|v| Json::Int(v as i64)).unwrap_or(Json::Null))]),
        },
    );
    m.push(Manifest::summary_entry());
    m.entries.extend(sctransform::output::manifest_entries("", counts.ncol, s.genes.len(), s.genes_step1.len(), sct.residuals.top_features.len(), sct.corrected.as_ref()));
    if context.is_none() {
        m.push(rna_entry(cluster::metadata::FILE, counts.ncol, rna.n_columns()));
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
    let mut result = sctransform::output::result_json(&sct, counts.nrow);
    // percent.mt's record (M15b), as scrust cluster's result.sctransform holds it
    if let (Json::Object(m), Some((p, _))) = (&mut result, &percent_mt) {
        m.push(("percent_mt".into(), p.record_json()));
    }
    let mut summary = Json::obj(vec![
        ("tool", Json::Str("scrust sctransform".into())),
        ("reference", Json::Str(sctransform::REFERENCE.into())),
        ("scrust_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        ("parameters", params_json(&cp)),
        ("result", result),
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
    // the warnings about percent.mt and the regression (M15b), only when there are some
    let mut notes: Vec<String> = percent_mt.as_ref().map(|(_, n)| n.clone()).unwrap_or_default();
    if let Some(r) = &sct.regression {
        let (rank, p) = r.rank();
        if rank < p && notes.iter().all(|n| !n.starts_with("vars.to.regress")) {
            notes.push(format!("vars.to.regress: the model matrix has rank {rank} of {p}"));
        }
    }
    if !notes.is_empty() {
        if let Json::Object(m) = &mut summary {
            if let Some((_, Json::Object(run))) = m.iter_mut().find(|(k, _)| k == "run") {
                run.push(("notes".into(), Json::strs(&notes)));
            }
        }
    }
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
