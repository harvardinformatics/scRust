//! The workflow's DoubletFinder call sequence (docs/PLAN.md B1; `validation/R/lib/doubletfinder.R`)
//! over a source of embeddings:
//!
//! 1. `paramSweep(PCs = 1:10, sct = TRUE)`: the real cells (a 10,000-cell subsample above 10,000
//!    cells), then for each pN the embeddings of real cells plus doublets, and pANN at every pK;
//! 2. `summarizeSweep(GT = FALSE)` and `find.pK`; pK = `as.numeric(as.character(...))` of the
//!    row with the largest BCmetric;
//! 3. nExp with the homotypic adjustment;
//! 4. `doubletFinder(pN = 0.25, pK, nExp)`: the final embeddings, pANN, and the calls.
//!
//! What produces the embeddings is behind [`EmbeddingSource`]: in milestone M3 R's own (stage
//! dumps, [`crate::stages::StageSource`]); from M6 on, doublets, SCTransform and PCA in Rust,
//! drawing from R's random-number stream in the same order as R.

use crate::pann::{classify, final_pann, sweep_iteration, Embeddings, FinalRun, NeighbourDumps, SweepIteration};
use crate::params::{column_name, n_doublets, pk_grid, pn_grid, sweep_name, Labelled};
use crate::summary::{find_pk, model_homotypic, n_exp, summarize_sweep, which_max, BcmvnRow, SweepStat};
use crate::DfError;
use rcompat::{as_character, as_numeric};
use scio::{fmt17, Json, Table};
use std::io;
use std::path::Path;

/// Produces the PCA embeddings DoubletFinder computes distances on. Called in R's order:
/// [`sweep_cells`](Self::sweep_cells) once, [`sweep_embeddings`](Self::sweep_embeddings) for
/// each pN in order, then [`final_embeddings`](Self::final_embeddings), so an implementation
/// that draws random numbers can follow R's single stream.
pub trait EmbeddingSource {
    /// paramSweep's real cells as indices into the input cells, in R's order: all cells, or
    /// `sample(1:n, 10000)` when there are more than 10,000.
    fn sweep_cells(&mut self, n: usize) -> Result<Vec<usize>, DfError>;
    /// The embeddings of sweep iteration `pn_index` (pN = `pn`): the real cells in the order of
    /// `real_cells`, then `n_doublets` artificial doublets.
    fn sweep_embeddings(&mut self, pn_index: usize, pn: f64, real_cells: &[usize], n_doublets: usize) -> Result<Embeddings, DfError>;
    /// The final run's embeddings: every input cell in order, then `n_doublets` doublets.
    fn final_embeddings(&mut self, pn: f64, n_doublets: usize) -> Result<Embeddings, DfError>;
}

/// How pK is chosen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PkChoice {
    /// The workflow's choice: the pK with the largest BCmetric of the sweep.
    Auto,
    /// A given value (R's `pK = <value>`); the sweep is skipped.
    Fixed(f64),
}

/// The workflow's parameters (plan B3: `--pcs`, `--pn`, `--nexp-rate`, `--pk`, `--seed`).
#[derive(Debug, Clone)]
pub struct Params {
    pub pcs: usize,
    pub pn: f64,
    pub nexp_rate: f64,
    pub pk: PkChoice,
    /// Recorded in summary.json; the RNG itself is the source's.
    pub seed: i32,
    pub threads: usize,
    /// Keep per-stage intermediate values (for the stage dumps).
    pub keep_stages: bool,
    /// What to keep of the neighbour ranking.
    pub neighbour_dumps: NeighbourDumps,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            pcs: 10,
            pn: 0.25,
            nexp_rate: 0.15,
            pk: PkChoice::Auto,
            seed: 12345,
            threads: 1,
            keep_stages: false,
            neighbour_dumps: NeighbourDumps::default(),
        }
    }
}

/// What DoubletFinder reads besides the embeddings: the cells, their cluster labels and the
/// factor's levels (for the homotypic adjustment), and the number of genes (reported).
#[derive(Debug, Clone)]
pub struct Input {
    pub barcodes: Vec<String>,
    pub clusters: Vec<String>,
    pub cluster_levels: Vec<String>,
    pub n_genes: usize,
}

/// The sweep and its summaries.
#[derive(Debug, Clone)]
pub struct Sweep {
    pub real_cells: Vec<usize>,
    pub pn: Vec<Labelled>,
    pub pk: Vec<Labelled>,
    pub iterations: Vec<SweepIteration>,
    pub stats: Vec<SweepStat>,
    pub bcmvn: Vec<BcmvnRow>,
    /// 0-based row of bcmvn with the largest BCmetric.
    pub best_row: usize,
}

/// Everything a run produces.
#[derive(Debug, Clone)]
pub struct Output {
    pub n_cells: usize,
    pub n_genes: usize,
    pub params: Params,
    pub sweep: Option<Sweep>,
    /// The pK of the final run, and its label.
    pub pk: f64,
    pub pk_label: String,
    pub nexp_poi: f64,
    pub homotypic_prop: f64,
    pub nexp: f64,
    pub n_doublets: usize,
    pub final_run: FinalRun,
    pub doublet: Vec<bool>,
    pub barcodes: Vec<String>,
    pub pann_column: String,
    pub classification_column: String,
}

/// Runs DoubletFinder as the workflow calls it. `mark(phase)` is called after each phase
/// (for timings and memory).
pub fn run(input: &Input, params: &Params, source: &mut dyn EmbeddingSource, mark: &mut dyn FnMut(&str)) -> Result<Output, DfError> {
    let n = input.barcodes.len();
    if input.clusters.len() != n {
        return Err(DfError::Input(format!("{} cluster labels for {n} cells", input.clusters.len())));
    }
    let threads = params.threads.max(1);

    let (sweep, pk, pk_label) = match params.pk {
        PkChoice::Fixed(pk) => (None, pk, as_character(pk)),
        PkChoice::Auto => {
            let pk_values = pk_grid(n);
            if pk_values.is_empty() {
                return Err(DfError::Input(format!("no pK value passes round(pK * min.cells) >= 1 for {n} cells")));
            }
            let pn: Vec<Labelled> = pn_grid().into_iter().map(Labelled::new).collect();
            let pk: Vec<Labelled> = pk_values.iter().map(|&v| Labelled::new(v)).collect();
            let real_cells = source.sweep_cells(n)?;
            if real_cells.is_empty() || real_cells.iter().any(|&i| i >= n) {
                return Err(DfError::Input("sweep cells out of range".into()));
            }
            let n_real = real_cells.len();
            let mut iterations = Vec::with_capacity(pn.len());
            for (j, p) in pn.iter().enumerate() {
                let nd = n_doublets(n_real, p.value);
                let emb = source.sweep_embeddings(j, p.value, &real_cells, nd)?;
                if emb.n != n_real + nd {
                    return Err(DfError::Input(format!("sweep pN = {}: {} embedded cells, expected {} real + {nd} doublets", p.label, emb.n, n_real)));
                }
                iterations.push(sweep_iteration(&emb, params.pcs, n_real, &pk_values, threads, params.neighbour_dumps)?);
                mark(&format!("sweep_pN{}", p.label));
            }
            mark("paramSweep");
            let pann: Vec<Vec<Vec<f64>>> = iterations.iter().map(|it| it.pann.clone()).collect();
            let stats = summarize_sweep(&pann, &pn, &pk, params.keep_stages, threads)?;
            mark("summarizeSweep");
            let bcmvn = find_pk(&stats);
            let metric: Vec<f64> = bcmvn.iter().map(|r| r.bc_metric).collect();
            let best_row = which_max(&metric).ok_or_else(|| DfError::Input("every BCmetric is NaN: no pK can be chosen".into()))?;
            mark("find.pK");
            // optimal_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
            let pk_final = as_numeric(&bcmvn[best_row].pk_level);
            let label = as_character(pk_final);
            (Some(Sweep { real_cells, pn, pk, iterations, stats, bcmvn, best_row }), pk_final, label)
        }
    };

    let homotypic_prop = model_homotypic(&input.clusters, &input.cluster_levels);
    let (nexp_poi, nexp) = n_exp(n, params.nexp_rate, homotypic_prop);

    let nd = n_doublets(n, params.pn);
    let emb = source.final_embeddings(params.pn, nd)?;
    if emb.n != n + nd {
        return Err(DfError::Input(format!("final run: {} embedded cells, expected {n} + {nd} doublets", emb.n)));
    }
    let final_run = final_pann(&emb, params.pcs, n, pk, threads, params.neighbour_dumps)?;
    drop(emb);
    let doublet = classify(&final_run.pann, nexp)?;
    mark("doubletFinder");

    Ok(Output {
        n_cells: n,
        n_genes: input.n_genes,
        params: params.clone(),
        sweep,
        pk,
        pk_label,
        nexp_poi,
        homotypic_prop,
        nexp,
        n_doublets: nd,
        pann_column: column_name("pANN", params.pn, pk, nexp),
        classification_column: column_name("DF.classifications", params.pn, pk, nexp),
        final_run,
        doublet,
        barcodes: input.barcodes.clone(),
    })
}

/// A length-1 character vector is a JSON scalar in R's `to_json`, longer ones are arrays.
fn r_strs(v: &[String]) -> Json {
    if v.len() == 1 {
        Json::Str(v[0].clone())
    } else {
        Json::strs(v)
    }
}

impl Output {
    /// `cells.tsv`: barcode, pANN and classification of every input cell.
    pub fn cells_table(&self) -> Table {
        Table {
            header: vec!["barcode".into(), "pANN".into(), "classification".into()],
            rows: self
                .barcodes
                .iter()
                .zip(&self.final_run.pann)
                .zip(&self.doublet)
                .map(|((b, &p), &d)| vec![b.clone(), fmt17(p), if d { "Doublet" } else { "Singlet" }.into()])
                .collect(),
        }
    }

    /// `sweep_stats.tsv`: pN and pK as their factor levels, and BCreal.
    pub fn sweep_stats_table(&self) -> Option<Table> {
        let sweep = self.sweep.as_ref()?;
        Some(Table {
            header: vec!["pN".into(), "pK".into(), "BCreal".into()],
            rows: sweep.stats.iter().map(|s| vec![s.pn.level(), s.pk.level(), fmt17(s.bc_real)]).collect(),
        })
    }

    /// `bcmvn.tsv`: find.pK's table.
    pub fn bcmvn_table(&self) -> Option<Table> {
        let sweep = self.sweep.as_ref()?;
        Some(Table {
            header: ["ParamID", "pK", "MeanBC", "VarBC", "BCmetric"].iter().map(|s| s.to_string()).collect(),
            rows: sweep
                .bcmvn
                .iter()
                .map(|r| vec![r.param_id.to_string(), r.pk_level.clone(), fmt17(r.mean_bc), fmt17(r.var_bc), fmt17(r.bc_metric)])
                .collect(),
        })
    }

    /// The `result` section of `summary.json`, with R's fields in R's order. `rng_end_sha256`
    /// (the digest of the final `.Random.seed`) is null when no random numbers were drawn
    /// (teacher forcing).
    pub fn result_json(&self, rng_end_sha256: Option<String>) -> Json {
        let sweep = self.sweep.as_ref();
        let pk_labels: Vec<String> = sweep.map(|s| s.pk.iter().map(|p| p.label.clone()).collect()).unwrap_or_default();
        let pn_labels: Vec<String> = sweep.map(|s| s.pn.iter().map(|p| p.label.clone()).collect()).unwrap_or_default();
        let n = self.n_cells as f64;
        let nd = self.n_doublets as f64;
        Json::obj(vec![
            ("n_cells", Json::Int(self.n_cells as i64)),
            ("n_genes", Json::Int(self.n_genes as i64)),
            ("seed", Json::Int(self.params.seed as i64)),
            ("PCs", Json::Int(self.params.pcs as i64)),
            ("pN", Json::Num(self.params.pn)),
            ("sweep_pK_labels", r_strs(&pk_labels)),
            ("sweep_pN_labels", r_strs(&pn_labels)),
            ("bcmvn_max_row", sweep.map(|s| Json::Int(s.best_row as i64 + 1)).unwrap_or(Json::Null)),
            ("pK", Json::Num(self.pk)),
            ("pK_label", Json::Str(self.pk_label.clone())),
            ("nExp_poi", Json::Num(self.nexp_poi)),
            ("homotypic_prop", Json::Num(self.homotypic_prop)),
            ("nExp", Json::Num(self.nexp)),
            ("n_doublets", Json::Num(nd)),
            ("n_cells_with_doublets", Json::Num(n + nd)),
            ("k", Json::Num(self.final_run.k)),
            ("pANN_column", Json::Str(self.pann_column.clone())),
            ("classification_column", Json::Str(self.classification_column.clone())),
            ("n_doublet_calls", Json::Int(self.doublet.iter().filter(|&&d| d).count() as i64)),
            ("rng_end_sha256", rng_end_sha256.map(Json::Str).unwrap_or(Json::Null)),
        ])
    }

    /// The names R gives the sweep's results, `pN_<pN>_pK_<pK>`.
    pub fn sweep_names(&self) -> Vec<String> {
        self.sweep
            .as_ref()
            .map(|s| s.pn.iter().flat_map(|pn| s.pk.iter().map(move |pk| sweep_name(pn, pk))).collect())
            .unwrap_or_default()
    }

    /// Writes `cells.tsv`, `sweep_stats.tsv` and `bcmvn.tsv` to `dir` (the latter two only
    /// after a sweep).
    pub fn write_tables(&self, dir: &Path) -> io::Result<()> {
        self.cells_table().write(&dir.join("cells.tsv"))?;
        if let Some(t) = self.sweep_stats_table() {
            t.write(&dir.join("sweep_stats.tsv"))?;
        }
        if let Some(t) = self.bcmvn_table() {
            t.write(&dir.join("bcmvn.tsv"))?;
        }
        Ok(())
    }
}
