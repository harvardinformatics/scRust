//! Stage dumps in the layout of the instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`):
//! reading R's (teacher forcing, plan A5b) and writing Rust's, so that the two compare block by
//! block.
//!
//! A stage directory holds `.npy` arrays and `index.tsv`, one row per (array, context): `file`,
//! `ctx`, `row_start`, `nrow`, `ncol`, `dtype`, `stage`. Arrays of one quantity from all
//! contexts are stacked along rows; large matrices have one file per context. Contexts:
//! `sweep_pN<pN>` (pN as `paste()` writes it), `summary:<sweep result name>`, `final`.
//! Cells are 0-based rows of the embeddings (real cells, then doublets); order positions are R's
//! 1-based ones.

use crate::pann::Embeddings;
use crate::run::{EmbeddingSource, Output};
use crate::DfError;
use rcompat::as_character;
use scio::npy::{self, Dtype, Npy};
use scio::Table;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// One block of an array: rows `row_start..row_start + nrow` of `file`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub file: String,
    pub ctx: String,
    pub row_start: usize,
    pub nrow: usize,
    pub ncol: usize,
    pub dtype: String,
    pub stage: String,
}

/// A stage directory written by an instrumented R run (or by Rust).
#[derive(Debug, Clone)]
pub struct StageDir {
    pub dir: PathBuf,
    pub blocks: Vec<Block>,
}

impl StageDir {
    pub fn open(dir: &Path) -> Result<StageDir, DfError> {
        let t = Table::read(&dir.join("index.tsv")).map_err(DfError::Io)?;
        let want = ["file", "ctx", "row_start", "nrow", "ncol", "dtype", "stage"];
        if t.header != want {
            return Err(DfError::Input(format!("{}: unexpected index.tsv header {:?}", dir.display(), t.header)));
        }
        let parse = |s: &str| s.parse::<usize>().map_err(|_| DfError::Input(format!("{}: bad index.tsv number {s}", dir.display())));
        let mut blocks = Vec::with_capacity(t.rows.len());
        for r in &t.rows {
            blocks.push(Block {
                file: r[0].clone(),
                ctx: r[1].clone(),
                row_start: parse(&r[2])?,
                nrow: parse(&r[3])?,
                ncol: parse(&r[4])?,
                dtype: r[5].clone(),
                stage: r[6].clone(),
            });
        }
        Ok(StageDir { dir: dir.to_path_buf(), blocks })
    }

    /// The block of `file` for `ctx`.
    pub fn block(&self, file: &str, ctx: &str) -> Option<&Block> {
        self.blocks.iter().find(|b| b.file == file && b.ctx == ctx)
    }

    /// Every block of `file`, in row order.
    pub fn blocks_of(&self, file: &str) -> Vec<&Block> {
        let mut v: Vec<&Block> = self.blocks.iter().filter(|b| b.file == file).collect();
        v.sort_by_key(|b| b.row_start);
        v
    }

    pub fn read(&self, file: &str) -> Result<Npy, DfError> {
        Npy::read(&self.dir.join(file)).map_err(|e| DfError::Input(e.to_string()))
    }
}

/// Teacher forcing: R's own embeddings (stage S3) and sweep subsample (stage S1) from an
/// instrumented run's stage directory.
pub struct StageSource {
    stages: StageDir,
    embeddings: Option<Npy>,
}

impl StageSource {
    pub fn new(dir: &Path) -> Result<StageSource, DfError> {
        Ok(StageSource { stages: StageDir::open(dir)?, embeddings: None })
    }

    fn embeddings(&mut self, ctx: &str) -> Result<Embeddings, DfError> {
        if self.embeddings.is_none() {
            self.embeddings = Some(self.stages.read("S3_embeddings.npy")?);
        }
        let b = self
            .stages
            .block("S3_embeddings.npy", ctx)
            .ok_or_else(|| DfError::Input(format!("{}: no S3 embeddings for context {ctx}", self.stages.dir.display())))?
            .clone();
        let a = self.embeddings.as_ref().unwrap();
        Embeddings::new(a.block_f64(b.row_start, b.nrow), b.nrow, b.ncol)
    }
}

impl EmbeddingSource for StageSource {
    fn sweep_cells(&mut self, n: usize) -> Result<Vec<usize>, DfError> {
        match self.stages.block("S1_subsample.npy", "sweep") {
            None => Ok((0..n).collect()),
            Some(b) => {
                let b = b.clone();
                let a = self.stages.read("S1_subsample.npy")?;
                let cells: Vec<usize> = a.col_i32(0, b.row_start, b.nrow).into_iter().map(|i| i as usize).collect();
                if cells.iter().any(|&i| i >= n) {
                    return Err(DfError::Input("S1_subsample: cell index out of range".into()));
                }
                Ok(cells)
            }
        }
    }

    fn sweep_embeddings(&mut self, _pn_index: usize, pn: f64, _real_cells: &[usize], _n_doublets: usize) -> Result<Embeddings, DfError> {
        self.embeddings(&format!("sweep_pN{}", as_character(pn)))
    }

    fn final_embeddings(&mut self, _pn: f64, _n_doublets: usize) -> Result<Embeddings, DfError> {
        // the sweep's embeddings are not needed any more
        self.embeddings("final").inspect(|_| self.embeddings = None)
    }
}

/// One context's block of an array: column-major values, its shape and its stage.
struct ArrayBlock<T> {
    ctx: String,
    values: Vec<T>,
    nrow: usize,
    ncol: usize,
    stage: String,
}

/// Collects stage arrays and writes them in R's layout.
#[derive(Default)]
pub struct StageWriter {
    /// name -> blocks of a stacked array
    f64s: BTreeMap<String, Vec<ArrayBlock<f64>>>,
    i32s: BTreeMap<String, Vec<ArrayBlock<i32>>>,
    /// one file per context: (file, block)
    big: Vec<(String, ArrayBlock<i32>)>,
}

/// Stacks the blocks' rows, column by column (Fortran order): (values, total rows, columns).
fn stack<T: Copy>(blocks: &[ArrayBlock<T>]) -> (Vec<T>, usize, usize) {
    let (total, ncol) = (blocks.iter().map(|b| b.nrow).sum::<usize>(), blocks[0].ncol);
    let mut data = Vec::with_capacity(total * ncol);
    for c in 0..ncol {
        for b in blocks {
            data.extend_from_slice(&b.values[c * b.nrow..(c + 1) * b.nrow]);
        }
    }
    (data, total, ncol)
}

impl StageWriter {
    pub fn add_f64(&mut self, name: &str, ctx: &str, values: Vec<f64>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.f64s.entry(name.into()).or_default().push(ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() });
    }

    pub fn add_i32(&mut self, name: &str, ctx: &str, values: Vec<i32>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.i32s.entry(name.into()).or_default().push(ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() });
    }

    /// One matrix file per context (`<name>.<ctx>.npy`), as R writes `S4_positions`.
    pub fn add_big_i32(&mut self, name: &str, ctx: &str, values: Vec<i32>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.big.push((format!("{name}.{ctx}.npy"), ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() }));
    }

    /// Writes every array and `index.tsv` into `dir`.
    pub fn write(&self, dir: &Path) -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let mut index: Vec<Vec<String>> = Vec::new();
        let shape = |nrow: usize, ncol: usize| if ncol == 1 { vec![nrow] } else { vec![nrow, ncol] };
        let mut push_index = |file: &str, blocks: &[(&str, usize, usize, &str)], dtype: Dtype| {
            let mut start = 0;
            for &(ctx, nrow, ncol, stage) in blocks {
                index.push(vec![file.into(), ctx.into(), start.to_string(), nrow.to_string(), ncol.to_string(), dtype.descr().into(), stage.into()]);
                start += nrow;
            }
        };
        for (name, blocks) in &self.f64s {
            let file = format!("{name}.npy");
            let (data, total, ncol) = stack(blocks);
            npy::write_f64(&dir.join(&file), &data, &shape(total, ncol))?;
            let meta: Vec<(&str, usize, usize, &str)> = blocks.iter().map(|b| (b.ctx.as_str(), b.nrow, b.ncol, b.stage.as_str())).collect();
            push_index(&file, &meta, Dtype::F8);
        }
        for (name, blocks) in &self.i32s {
            let file = format!("{name}.npy");
            let (data, total, ncol) = stack(blocks);
            npy::write_i32(&dir.join(&file), &data, &shape(total, ncol))?;
            let meta: Vec<(&str, usize, usize, &str)> = blocks.iter().map(|b| (b.ctx.as_str(), b.nrow, b.ncol, b.stage.as_str())).collect();
            push_index(&file, &meta, Dtype::I4);
        }
        for (file, b) in &self.big {
            npy::write_i32(&dir.join(file), &b.values, &[b.nrow, b.ncol])?;
            push_index(file, &[(b.ctx.as_str(), b.nrow, b.ncol, b.stage.as_str())], Dtype::I4);
        }
        index.sort_by(|a, b| a[0].cmp(&b[0]).then(a[2].parse::<usize>().unwrap().cmp(&b[2].parse::<usize>().unwrap())));
        Table { header: ["file", "ctx", "row_start", "nrow", "ncol", "dtype", "stage"].iter().map(|s| s.to_string()).collect(), rows: index }
            .write(&dir.join("index.tsv"))
    }
}

/// Positions per cell (0-based), as R's L x n matrix of 1-based positions, column-major.
fn positions_matrix(per_cell: &[Vec<u32>]) -> (Vec<i32>, usize) {
    let l = per_cell.first().map_or(0, |p| p.len());
    let mut m = Vec::with_capacity(l * per_cell.len());
    for p in per_cell {
        debug_assert_eq!(p.len(), l);
        m.extend(p.iter().map(|&x| x as i32 + 1));
    }
    (m, l)
}

fn digest_columns(d: &[(u64, u64)]) -> Vec<f64> {
    let mut v: Vec<f64> = d.iter().map(|x| x.0 as f64).collect();
    v.extend(d.iter().map(|x| x.1 as f64));
    v
}

/// Rust's stage arrays S4-S7 for a run, in R's layout and with R's names: `S4_digest` (and
/// `S4_positions.<ctx>` when kept), `S5_pann_sweep`, `S5_pann_final`, `S6_bkde`, `S6_gkde`,
/// `S7_doublet`. S1-S3 (doublets, SCTransform, PCA) are not Rust's in teacher-forced runs and
/// are not written.
pub fn stage_dumps(out: &Output) -> StageWriter {
    let mut w = StageWriter::default();
    if let Some(sweep) = &out.sweep {
        for (it, pn) in sweep.iterations.iter().zip(&sweep.pn) {
            let ctx = format!("sweep_pN{}", pn.label);
            if let Some(d) = &it.digests {
                w.add_f64("S4_digest", &ctx, digest_columns(d), d.len(), 2, "S4");
            }
            if let Some(p) = &it.positions {
                let (m, l) = positions_matrix(p);
                w.add_big_i32("S4_positions", &ctx, m, l, p.len(), "S4");
            }
            let npk = it.pann.len();
            let mut v = Vec::with_capacity(it.n_real * npk);
            for col in &it.pann {
                v.extend_from_slice(col);
            }
            w.add_f64("S5_pann_sweep", &ctx, v, it.n_real, npk, "S5");
        }
        for (s, name) in sweep.stats.iter().zip(out.sweep_names()) {
            let ctx = format!("summary:{name}");
            if let Some(b) = &s.bkde {
                let mut v = b.x.clone();
                v.extend_from_slice(&b.y);
                w.add_f64("S6_bkde", &ctx, v, b.x.len(), 2, "S6");
            }
            if let Some(g) = &s.gkde {
                w.add_f64("S6_gkde", &ctx, g.clone(), g.len(), 1, "S6");
            }
        }
    }
    let f = &out.final_run;
    if let Some(d) = &f.digests {
        w.add_f64("S4_digest", "final", digest_columns(d), d.len(), 2, "S4");
    }
    if let Some(p) = &f.positions {
        let (m, l) = positions_matrix(p);
        w.add_big_i32("S4_positions", "final", m, l, p.len(), "S4");
    }
    w.add_f64("S5_pann_final", "final", f.pann.clone(), f.pann.len(), 1, "S5");
    w.add_i32("S7_doublet", "final", out.doublet.iter().map(|&d| d as i32).collect(), out.doublet.len(), 1, "S7");
    w
}
