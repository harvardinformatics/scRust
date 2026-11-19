//! The clustering chain's parameters (docs/PLAN.md A7b, "Built for parameter optimisation", and
//! C3, "Parameters"): every Seurat argument A7b lists, with Seurat 5.5.1's default, set by a
//! preset, a TOML file (`scrust cluster --params`) and command-line flags, in that order.
//!
//! - **SCTransform:** `seed.use` (1448145), `ncells` (5000), `variable.features.n` (3000),
//!   `clip.range` (+-sqrt(N/30)).
//! - **RunPCA:** `npcs` (50), `seed.use` (42), `weight.by.var` (TRUE); and the workflow's guard
//!   `require_min_cells_for_pca()` (it stops with 50 cells or fewer), which the preset turns on.
//! - **FindNeighbors:** `dims` (1:10), `k.param` (20), `n.trees` (50), `annoy.metric`
//!   (euclidean), `prune.SNN` (1/15), and AnnoySearch's `search.k` (-1).
//! - **FindClusters:** `resolution` (0.8, several allowed), `algorithm` (1), `n.start` (10),
//!   `n.iter` (10), `random.seed` (0), `modularity.fxn` (1), `group.singletons` (TRUE),
//!   `graph.name` (SCT_snn).
//! - **The scores** (`add_silhouette_to_metadata`): silhouette and purity, both by default, on the
//!   neighbour search's dims; purity's `k` (50) and `weighted` (TRUE). Purity comes with
//!   milestone M12.
//!
//! **`workflow`** is the workflow's calls (C1): seed 12345 in SCTransform, RunPCA and
//! FindClusters, `dims = 1:30`, and the PCA guard.
//!
//! **`vars.to.regress` and `percent.mt`** (milestone M15b; [`crate::covariates`]): SCTransform's
//! `vars.to.regress` (`--vars-to-regress`, `[sctransform] vars_to_regress`; none by default),
//! `percent.mt` asked for without regressing it (`--percent-mt`, `percent_mt`), and its
//! mitochondrial features given by id (`--mt-genes <file>`, `mt_genes`; a relative path in a
//! parameter file is read from the file's directory). `percent.mt` is computed when any of them
//! asks for it.
//!
//! **Other code paths are refused** with a message, as [`sctransform::SctParams::check`] does:
//! another `vst.flavor` or `method`, `do.scale`, `nn.method = "rann"` or `"hnsw"`,
//! `l2.norm = TRUE`, `annoy.metric = "hamming"`, Leiden (`algorithm = 4`), and seeds of `NULL`
//! (which would leave R's stream where an earlier session left it).
//!
//! **Numbers are read as R reads them:** each number's text goes through R's own decimal parser
//! (`as.numeric`), and a fraction written as text, `"1/15"`, is divided in double, as R
//! evaluates `1/15`.

use crate::params::{parse_resolutions, ClusterParams};
use neighbors::{parse_dims, parse_number, Metric, NeighborParams};
use pca::PcaParams;
use scio::toml::{self, Table, Value};
use sctransform::SctParams;

/// The workflow's seed (`SCRNASEQ_PREPROCESS_SEED`).
pub const WORKFLOW_SEED: i32 = 12345;

/// The workflow guard's `npcs` (`require_min_cells_for_pca(seurat_obj, npcs = 50L)`).
pub const WORKFLOW_GUARD_NPCS: usize = 50;

/// Which scores to compute, and purity's arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct Stats {
    pub silhouette: bool,
    pub purity: bool,
    /// `neighborPurity(k = 50)`.
    pub purity_k: usize,
    /// `neighborPurity(weighted = TRUE)`.
    pub purity_weighted: bool,
}

impl Default for Stats {
    fn default() -> Self {
        Stats { silhouette: true, purity: true, purity_k: 50, purity_weighted: true }
    }
}

impl Stats {
    /// The names asked for, in the workflow's order.
    pub fn names(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        if self.silhouette {
            v.push("silhouette");
        }
        if self.purity {
            v.push("purity");
        }
        v
    }

    /// `silhouette,purity`, `silhouette`, `purity` or `none`.
    pub fn parse_list(s: &str) -> Result<(bool, bool), String> {
        if s.trim() == "none" {
            return Ok((false, false));
        }
        let (mut sil, mut pur) = (false, false);
        for t in s.split(',') {
            match t.trim() {
                "silhouette" => sil = true,
                "purity" => pur = true,
                other => return Err(format!("stats: '{other}' is not silhouette, purity or none")),
            }
        }
        Ok((sil, pur))
    }
}

/// Every parameter of the chain.
#[derive(Debug, Clone, PartialEq)]
pub struct ChainParams {
    /// The preset applied (`None` for Seurat's defaults).
    pub preset: Option<String>,
    /// SCTransform's arguments, `vars.to.regress` among them (M15b).
    pub sct: SctParams,
    /// `percent.mt` asked for without regressing it (`--percent-mt`; M15b).
    pub percent_mt: bool,
    /// `--mt-genes`: the file listing the mitochondrial features' ids (M15b); the command reads it
    /// and hands its ids to the chain ([`crate::covariates::CellCovariates`]).
    pub mt_genes: Option<String>,
    pub pca: PcaParams,
    /// The workflow's `require_min_cells_for_pca()`: stop with this many cells or fewer.
    pub pca_guard: Option<usize>,
    pub neighbors: NeighborParams,
    pub cluster: ClusterParams,
    pub stats: Stats,
}

impl Default for ChainParams {
    /// Seurat 5.5.1's defaults throughout.
    fn default() -> Self {
        ChainParams {
            preset: None,
            sct: SctParams::default(),
            percent_mt: false,
            mt_genes: None,
            pca: PcaParams::default(),
            pca_guard: None,
            neighbors: NeighborParams::default(),
            cluster: ClusterParams::default(),
            stats: Stats::default(),
        }
    }
}

impl ChainParams {
    /// The workflow's calls (`--preset workflow`).
    pub fn workflow() -> ChainParams {
        let mut p = ChainParams::default();
        p.apply_preset("workflow").expect("the workflow preset");
        p
    }

    /// Applies a preset by name over the current values.
    pub fn apply_preset(&mut self, name: &str) -> Result<(), String> {
        match name {
            "workflow" => {
                self.sct.seed_use = Some(WORKFLOW_SEED);
                self.pca.seed_use = Some(WORKFLOW_SEED);
                self.pca_guard = Some(WORKFLOW_GUARD_NPCS);
                self.neighbors.dims = (1..=30).collect();
                self.cluster.random_seed = WORKFLOW_SEED;
                self.preset = Some("workflow".into());
                Ok(())
            }
            "seurat" | "none" => Ok(()),
            other => Err(format!("unknown preset '{other}' (workflow)")),
        }
    }

    /// Whether `percent.mt` is computed (M15b): asked for, given a list of features, or regressed.
    pub fn wants_percent_mt(&self) -> bool {
        self.percent_mt || self.mt_genes.is_some() || self.sct.vars_to_regress.iter().any(|v| v == crate::covariates::PERCENT_MT)
    }

    /// Whether the chain needs more than the counts (M15b): `percent.mt`, or covariates.
    pub fn needs_covariates(&self) -> bool {
        self.wants_percent_mt() || !self.sct.vars_to_regress.is_empty()
    }

    /// Checks every stage's parameters and refuses the paths that are not ported.
    pub fn check(&self) -> Result<(), String> {
        self.sct.check().map_err(|e| format!("SCTransform: {e}"))?;
        crate::covariates::check_names(&self.sct.vars_to_regress).map_err(|e| format!("SCTransform: {e}"))?;
        if self.sct.seed_use.is_none() {
            return Err(refused("SCTransform's seed.use = NULL", "it leaves R's stream where an earlier session left it; give a seed"));
        }
        if self.sct.ncells == 0 {
            return Err("SCTransform: ncells must be positive".into());
        }
        if self.sct.variable_features_n == Some(0) {
            return Err("SCTransform: variable.features.n must be positive".into());
        }
        if let Some((lo, hi)) = self.sct.clip_range {
            if !(lo.is_finite() && hi.is_finite() && lo < hi) {
                return Err(format!("SCTransform: clip.range c({lo}, {hi}) must be two finite numbers, the lower first"));
            }
        }
        if self.pca.seed_use.is_none() {
            return Err(refused("RunPCA's seed.use = NULL", "it leaves R's stream where SCTransform left it before irlba draws; give a seed"));
        }
        if self.pca.npcs == 0 {
            return Err("RunPCA: npcs must be positive".into());
        }
        let ncol = self.pca.npcs;
        self.neighbors.check(ncol).map_err(|e| format!("FindNeighbors: {e}"))?;
        self.cluster.check().map_err(|e| format!("FindClusters: {e}"))?;
        if self.stats.purity_k == 0 {
            return Err("neighborPurity: k must be positive".into());
        }
        Ok(())
    }
}

fn refused(what: &str, why: &str) -> String {
    format!("{what} is refused: {why}")
}

/// The flags `scrust cluster` takes for the parameters (each `--<name> <value>`, except
/// [`PARAM_SWITCHES`]).
pub const PARAM_FLAGS: &[&str] = &[
    // SCTransform
    "sct-seed", "ncells", "variable-features", "clip-range", "vars-to-regress", "mt-genes", "vst-flavor", "sct-method", "do-scale",
    // RunPCA
    "npcs", "pca-seed", "weight-by-var", "pca-guard",
    // FindNeighbors
    "dims", "k", "trees", "metric", "prune", "search-k", "nn-method", "l2-norm",
    // FindClusters
    "resolution", "algorithm", "starts", "iterations", "seed", "modularity", "graph-name",
    // the scores
    "stats", "purity-k", "purity-weighted",
];

/// Flags without a value.
pub const PARAM_SWITCHES: &[&str] = &["no-group-singletons", "percent-mt"];

fn int_of(name: &str, s: &str) -> Result<i64, String> {
    s.trim().parse::<i64>().map_err(|_| format!("{name}: '{s}' is not an integer"))
}

fn i32_of(name: &str, s: &str) -> Result<i32, String> {
    crate::params::parse_int(name, s)
}

fn usize_of(name: &str, s: &str) -> Result<usize, String> {
    let v = int_of(name, s)?;
    usize::try_from(v).map_err(|_| format!("{name}: {v} is negative"))
}

fn bool_of(name: &str, s: &str) -> Result<bool, String> {
    match s.trim() {
        "true" | "TRUE" | "T" | "yes" => Ok(true),
        "false" | "FALSE" | "F" | "no" => Ok(false),
        _ => Err(format!("{name}: '{s}' is not true or false")),
    }
}

/// `clip.range` as text: `-3,3` (two numbers read as R reads them).
fn clip_of(name: &str, s: &str) -> Result<(f64, f64), String> {
    let v: Vec<&str> = s.split(',').collect();
    if v.len() != 2 {
        return Err(format!("{name}: expected two numbers, lower,upper"));
    }
    Ok((parse_number(v[0].trim()).map_err(|e| format!("{name}: {e}"))?, parse_number(v[1].trim()).map_err(|e| format!("{name}: {e}"))?))
}

/// Sets one parameter from its text, by its flag name. Every route (flags, TOML) ends here, so a
/// value means the same wherever it is given.
pub fn set(p: &mut ChainParams, name: &str, v: &str) -> Result<(), String> {
    let flag = format!("--{name}");
    let n = flag.as_str();
    match name {
        "sct-seed" => p.sct.seed_use = Some(i32_of(n, v)?),
        "ncells" => p.sct.ncells = usize_of(n, v)?,
        "variable-features" => p.sct.variable_features_n = Some(usize_of(n, v)?),
        "clip-range" => p.sct.clip_range = Some(clip_of(n, v)?),
        "vars-to-regress" => {
            // comma-separated names, in the formula's order; empty for none
            let vars: Vec<String> = v.split(',').map(str::trim).filter(|t| !t.is_empty()).map(str::to_string).collect();
            crate::covariates::check_names(&vars).map_err(|e| format!("{n}: {e}"))?;
            p.sct.vars_to_regress = vars;
        }
        "mt-genes" => {
            if v.trim().is_empty() {
                return Err(format!("{n}: a file is needed"));
            }
            p.mt_genes = Some(v.to_string());
        }
        "percent-mt" => p.percent_mt = bool_of(n, v)?,
        "vst-flavor" => {
            if v.trim() != "v2" {
                return Err(refused(&format!("vst.flavor = \"{v}\""), "only sctransform's v2 flavour (Seurat's default) is ported"));
            }
        }
        "sct-method" => {
            if v.trim() != "glmGamPoi_offset" {
                return Err(refused(&format!("SCTransform's method = \"{v}\""), "only glmGamPoi_offset (the v2 flavour's) is ported"));
            }
        }
        "do-scale" => {
            if bool_of(n, v)? {
                return Err(refused("SCTransform's do.scale = TRUE", "scaling the residuals is not ported"));
            }
        }
        "npcs" => p.pca.npcs = usize_of(n, v)?,
        "pca-seed" => p.pca.seed_use = Some(i32_of(n, v)?),
        "weight-by-var" => p.pca.weight_by_var = bool_of(n, v)?,
        "pca-guard" => {
            p.pca_guard = match v.trim() {
                "none" | "false" | "no" => None,
                "true" | "yes" => Some(WORKFLOW_GUARD_NPCS),
                t => Some(usize_of(n, t)?),
            }
        }
        "dims" => p.neighbors.dims = parse_dims(v)?,
        "k" => p.neighbors.k_param = usize_of(n, v)?,
        "trees" => p.neighbors.n_trees = i32_of(n, v)?,
        "metric" => p.neighbors.metric = Metric::parse(v.trim())?,
        "prune" => p.neighbors.prune_snn = parse_number(v.trim()).map_err(|e| format!("{n}: {e}"))?,
        "search-k" => p.neighbors.search_k = i32_of(n, v)?,
        "nn-method" => {
            if v.trim() != "annoy" {
                return Err(refused(&format!("nn.method = \"{v}\""), "only Annoy (Seurat's default) is ported"));
            }
        }
        "l2-norm" => {
            if bool_of(n, v)? {
                return Err(refused("l2.norm = TRUE", "normalising the embeddings first is not ported"));
            }
        }
        "resolution" => p.cluster.resolutions = parse_resolutions(v)?,
        "algorithm" => p.cluster.algorithm = i32_of(n, v)?,
        "starts" => p.cluster.n_start = i32_of(n, v)?,
        "iterations" => p.cluster.n_iter = i32_of(n, v)?,
        "seed" => p.cluster.random_seed = i32_of(n, v)?,
        "modularity" => p.cluster.modularity_fxn = i32_of(n, v)?,
        "graph-name" => p.cluster.graph_name = v.to_string(),
        "group-singletons" => p.cluster.group_singletons = bool_of(n, v)?,
        "stats" => {
            let (s, q) = Stats::parse_list(v)?;
            p.stats.silhouette = s;
            p.stats.purity = q;
        }
        "purity-k" => p.stats.purity_k = usize_of(n, v)?,
        "purity-weighted" => p.stats.purity_weighted = bool_of(n, v)?,
        _ => return Err(format!("unknown parameter {name}")),
    }
    Ok(())
}

/// The TOML file's sections and keys, and the flag each one sets.
const TOML_KEYS: &[(&str, &str, &str)] = &[
    ("sctransform", "seed", "sct-seed"),
    ("sctransform", "ncells", "ncells"),
    ("sctransform", "variable_features", "variable-features"),
    ("sctransform", "clip_range", "clip-range"),
    ("sctransform", "vars_to_regress", "vars-to-regress"),
    ("sctransform", "mt_genes", "mt-genes"),
    ("sctransform", "percent_mt", "percent-mt"),
    ("sctransform", "vst_flavor", "vst-flavor"),
    ("sctransform", "method", "sct-method"),
    ("sctransform", "do_scale", "do-scale"),
    ("pca", "npcs", "npcs"),
    ("pca", "seed", "pca-seed"),
    ("pca", "weight_by_var", "weight-by-var"),
    ("pca", "guard", "pca-guard"),
    ("neighbors", "dims", "dims"),
    ("neighbors", "k", "k"),
    ("neighbors", "trees", "trees"),
    ("neighbors", "metric", "metric"),
    ("neighbors", "prune", "prune"),
    ("neighbors", "search_k", "search-k"),
    ("neighbors", "nn_method", "nn-method"),
    ("neighbors", "l2_norm", "l2-norm"),
    ("cluster", "resolution", "resolution"),
    ("cluster", "algorithm", "algorithm"),
    ("cluster", "starts", "starts"),
    ("cluster", "iterations", "iterations"),
    ("cluster", "seed", "seed"),
    ("cluster", "modularity", "modularity"),
    ("cluster", "group_singletons", "group-singletons"),
    ("cluster", "graph_name", "graph-name"),
    ("stats", "compute", "stats"),
    ("stats", "purity_k", "purity-k"),
    ("stats", "purity_weighted", "purity-weighted"),
];

/// A TOML value as the text [`set`] reads: numbers as written (so R's parser reads them),
/// strings as they are, booleans as `true`/`false`, arrays comma-separated (dims, resolutions,
/// clip.range, the scores).
fn value_text(key: &str, v: &Value) -> Result<String, String> {
    Ok(match v {
        Value::String(s) => s.clone(),
        Value::Integer(..) | Value::Float(..) => v.number_text().unwrap().to_string(),
        Value::Boolean(b) => b.to_string(),
        Value::Array(a) => {
            let parts: Result<Vec<String>, String> = a.iter().map(|x| value_text(key, x)).collect();
            parts?.join(",")
        }
        Value::Table(_) => return Err(format!("{key}: a table where a value was expected")),
    })
}

/// Reads a parameter file: an optional top-level `preset`, then the sections `sctransform`,
/// `pca`, `neighbors`, `cluster` and `stats` (keys in [`TOML_KEYS`]). Returns the preset named in
/// the file, if any, and the (flag, text) pairs to apply after it, in the file's order.
pub fn read_toml(text: &str) -> Result<(Option<String>, Vec<(String, String)>), String> {
    let t: Table = toml::parse(text).map_err(|e| e.to_string())?;
    let mut preset = None;
    let mut out = Vec::new();
    for (k, v) in &t.entries {
        match (k.as_str(), v) {
            ("preset", Value::String(s)) => preset = Some(s.clone()),
            ("preset", other) => return Err(format!("preset: expected a string, not {}", other.kind())),
            (section, Value::Table(s)) if ["sctransform", "pca", "neighbors", "cluster", "stats"].contains(&section) => {
                for (key, val) in &s.entries {
                    let Some(&(_, _, flag)) = TOML_KEYS.iter().find(|(sec, kk, _)| *sec == section && kk == key) else {
                        return Err(format!("[{section}] {key}: not a parameter of the chain (see docs/validation/clustering.md for the keys)"));
                    };
                    out.push((flag.to_string(), value_text(&format!("{section}.{key}"), val)?));
                }
            }
            (other, _) => return Err(format!("{other}: not a section of a parameter file (sctransform, pca, neighbors, cluster, stats, or preset = \"workflow\")")),
        }
    }
    Ok((preset, out))
}

/// Builds the parameters: Seurat's defaults, then the preset (`cli_preset`, else the file's),
/// then the file's values, then the flags (`(name, value)`, in the order given;
/// `no-group-singletons` and `percent-mt` as switches). Checks the result.
pub fn build(cli_preset: Option<&str>, toml_text: Option<&str>, flags: &[(String, String)]) -> Result<ChainParams, String> {
    build_in(cli_preset, toml_text, None, flags)
}

/// [`build`], with the parameter file's directory (`toml_dir`), from which a relative `mt_genes`
/// path in the file is read (M15b).
pub fn build_in(cli_preset: Option<&str>, toml_text: Option<&str>, toml_dir: Option<&std::path::Path>, flags: &[(String, String)]) -> Result<ChainParams, String> {
    let mut p = ChainParams::default();
    let (file_preset, file_values) = match toml_text {
        Some(t) => read_toml(t).map_err(|e| format!("--params: {e}"))?,
        None => (None, Vec::new()),
    };
    if let Some(name) = cli_preset.map(str::to_string).or(file_preset) {
        p.apply_preset(&name)?;
    }
    for (k, v) in &file_values {
        let v = match (k.as_str(), toml_dir) {
            ("mt-genes", Some(d)) if std::path::Path::new(v).is_relative() => d.join(v).display().to_string(),
            _ => v.clone(),
        };
        set(&mut p, k, &v).map_err(|e| format!("--params: {e}"))?;
    }
    for (k, v) in flags {
        if k == "no-group-singletons" {
            p.cluster.group_singletons = false;
        } else if k == "percent-mt" {
            p.percent_mt = true;
        } else {
            set(&mut p, k, v)?;
        }
    }
    p.check()?;
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect()
    }

    #[test]
    fn seurat_defaults_and_the_workflow_preset() {
        let d = build(None, None, &[]).unwrap();
        assert_eq!(d.sct.seed_use, Some(1448145));
        assert_eq!(d.pca.seed_use, Some(42));
        assert_eq!(d.pca.npcs, 50);
        assert_eq!(d.neighbors.dims, (1..=10).collect::<Vec<_>>());
        assert_eq!(d.neighbors.prune_snn.to_bits(), (1.0f64 / 15.0).to_bits());
        assert_eq!(d.cluster.random_seed, 0);
        assert_eq!(d.cluster.resolutions, vec![0.8]);
        assert_eq!(d.pca_guard, None);
        assert!(d.stats.silhouette && d.stats.purity);
        let w = build(Some("workflow"), None, &[]).unwrap();
        assert_eq!((w.sct.seed_use, w.pca.seed_use, w.cluster.random_seed), (Some(12345), Some(12345), 12345));
        assert_eq!(w.neighbors.dims.len(), 30);
        assert_eq!(w.pca_guard, Some(50));
        assert_eq!(w, ChainParams::workflow());
    }

    #[test]
    fn toml_then_flags() {
        let toml = "preset = \"workflow\"\n[neighbors]\nk = 30\nprune = \"1/15\"\ndims = \"1:20\"\n[cluster]\nresolution = [0.4, 1.2]\nseed = 7\ngroup_singletons = false\n[sctransform]\nclip_range = [-3, 3.5]\n[stats]\ncompute = \"silhouette\"\npurity_k = 30\n";
        let p = build(None, Some(toml), &f(&[("k", "10"), ("resolution", "0.8")])).unwrap();
        assert_eq!(p.preset.as_deref(), Some("workflow"));
        assert_eq!(p.neighbors.k_param, 10, "flags come last");
        assert_eq!(p.neighbors.dims, (1..=20).collect::<Vec<_>>());
        assert_eq!(p.cluster.resolutions, vec![0.8]);
        assert_eq!(p.cluster.random_seed, 7);
        assert!(!p.cluster.group_singletons);
        assert_eq!(p.sct.clip_range, Some((-3.0, 3.5)));
        assert!(p.stats.silhouette && !p.stats.purity);
        assert_eq!(p.stats.purity_k, 30);
        // a number in the file is read as R reads its text
        let p = build(None, Some("[cluster]\nresolution = 0.3\n"), &[]).unwrap();
        assert_eq!(p.cluster.resolutions[0].to_bits(), rcompat::as_numeric("0.3").to_bits());
        // the flag's preset wins over the file's
        let p = build(Some("seurat"), Some("preset = \"workflow\"\n"), &[]).unwrap();
        assert_eq!(p.cluster.random_seed, 0);
    }

    #[test]
    fn refusals() {
        let bad = |toml: &str, flags: &[(&str, &str)], needle: &str| {
            let e = build(None, if toml.is_empty() { None } else { Some(toml) }, &f(flags)).unwrap_err();
            assert!(e.contains(needle), "{toml} {flags:?}: {e}");
        };
        bad("", &[("algorithm", "4")], "Leiden");
        bad("", &[("nn-method", "rann")], "nn.method");
        bad("", &[("nn-method", "hnsw")], "only Annoy");
        bad("", &[("l2-norm", "true")], "l2.norm");
        bad("", &[("vars-to-regress", "percent.mt,percent.mt")], "twice");
        bad("", &[("vars-to-regress", "my-col")], "not a syntactic name");
        bad("", &[("metric", "hamming")], "hamming");
        bad("", &[("vst-flavor", "v1")], "vst.flavor");
        bad("", &[("do-scale", "TRUE")], "do.scale");
        bad("[neighbors]\nl2_norm = true\n", &[], "l2.norm");
        bad("[cluster]\nalgorithm = 4\n", &[], "Leiden");
        bad("[cluster]\nleiden = true\n", &[], "not a parameter");
        bad("[umap]\nseed = 1\n", &[], "not a section");
        bad("", &[("dims", "1:51")], "More dimensions");
        bad("", &[("resolution", "0.8,0.8")], "overwrite");
        bad("", &[("stats", "silhouette,entropy")], "entropy");
        bad("", &[("npcs", "0")], "npcs");
        assert!(build(Some("nonsense"), None, &[]).unwrap_err().contains("unknown preset"));
        // the default code path's own values are accepted
        assert!(build(None, None, &f(&[("nn-method", "annoy"), ("l2-norm", "false"), ("vst-flavor", "v2")])).is_ok());
    }

    #[test]
    fn vars_to_regress_and_percent_mt() {
        let d = build(Some("workflow"), None, &[]).unwrap();
        assert!(d.sct.vars_to_regress.is_empty() && !d.wants_percent_mt() && !d.needs_covariates());
        let p = build(None, None, &f(&[("vars-to-regress", " percent.mt , nCount_RNA")])).unwrap();
        assert_eq!(p.sct.vars_to_regress, vec!["percent.mt", "nCount_RNA"]);
        assert!(p.wants_percent_mt() && p.needs_covariates());
        let p = build(None, None, &f(&[("vars-to-regress", "nCount_RNA")])).unwrap();
        assert!(!p.wants_percent_mt() && p.needs_covariates());
        let p = build(None, None, &f(&[("percent-mt", "")])).unwrap();
        assert!(p.percent_mt && p.wants_percent_mt());
        // the file: a relative mt_genes path is the file's directory's
        let toml = "[sctransform]\nvars_to_regress = [\"percent.mt\", \"score\"]\nmt_genes = \"mt.txt\"\npercent_mt = true\n";
        let p = build_in(None, Some(toml), Some(std::path::Path::new("/some/dir")), &[]).unwrap();
        assert_eq!(p.sct.vars_to_regress, vec!["percent.mt", "score"]);
        assert_eq!(p.mt_genes.as_deref(), Some("/some/dir/mt.txt"));
        assert!(p.percent_mt);
        let p = build_in(None, Some("[sctransform]\nvars_to_regress = []\n"), None, &[]).unwrap();
        assert!(p.sct.vars_to_regress.is_empty());
    }
}
