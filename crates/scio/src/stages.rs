//! Stage dumps in the layout of the instrumented R runs (`validation/R/lib/doubletfinder_hooks.R`,
//! `validation/R/lib/clustering_hooks.R`): read (teacher forcing, docs/PLAN.md A5b) and written
//! (a Rust run's own dumps), so that a Rust run and an R run compare block by block, checkpoint
//! by checkpoint and scalar by scalar. Shared by DoubletFinder (M6) and the clustering chain
//! (M11).
//!
//! A stage directory holds `.npy` arrays and `index.tsv`, one row per (array, context): `file`,
//! `ctx`, `row_start`, `nrow`, `ncol`, `dtype`, `stage`. Arrays of one quantity from all contexts
//! are stacked along rows; large matrices have one file per context (`<name>.<ctx>.npy`).
//! `.Random.seed` checkpoints are the rows of `rng.npy` (one state per row), labelled by
//! `rng.tsv` (`row`, `ctx`, `stage`, `label`); per-context scalars are in `stages.json`.

use crate::json::Json;
use crate::npy::{self, Dtype, Npy};
use crate::tsv::Table;
use std::collections::BTreeMap;
use std::io;
use std::path::{Path, PathBuf};

/// What a run dumps, as the instrumented R runs' levels: `light` (no residual-sized matrix),
/// `standard` (also `scale.data`, the PCA input), `full` (also the residual matrix of every kept
/// gene, and whatever else a tool keeps only at this level).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DumpLevel {
    Light,
    Standard,
    Full,
}

impl DumpLevel {
    /// `none` (no dumps), `light`, `standard` or `full`.
    pub fn parse(s: &str) -> Result<Option<DumpLevel>, String> {
        match s {
            "none" => Ok(None),
            "light" => Ok(Some(DumpLevel::Light)),
            "standard" => Ok(Some(DumpLevel::Standard)),
            "full" => Ok(Some(DumpLevel::Full)),
            _ => Err(format!("unknown dump level '{s}' (none, light, standard or full)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            DumpLevel::Light => "light",
            DumpLevel::Standard => "standard",
            DumpLevel::Full => "full",
        }
    }
}

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

fn invalid(m: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, m)
}

/// A stage directory written by an instrumented R run (or by Rust).
#[derive(Debug, Clone)]
pub struct StageDir {
    pub dir: PathBuf,
    pub blocks: Vec<Block>,
}

impl StageDir {
    pub fn open(dir: &Path) -> io::Result<StageDir> {
        let t = Table::read(&dir.join("index.tsv")).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", dir.join("index.tsv").display())))?;
        let want = ["file", "ctx", "row_start", "nrow", "ncol", "dtype", "stage"];
        if t.header != want {
            return Err(invalid(format!("{}: unexpected index.tsv header {:?}", dir.display(), t.header)));
        }
        let parse = |s: &str| s.parse::<usize>().map_err(|_| invalid(format!("{}: bad index.tsv number {s}", dir.display())));
        let mut blocks = Vec::with_capacity(t.rows.len());
        for r in &t.rows {
            if r.len() != want.len() {
                return Err(invalid(format!("{}: index.tsv row of {} fields", dir.display(), r.len())));
            }
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

    /// Whether the index lists `file` for `ctx` and the file is there (R's tier-T fixtures list
    /// some large files they do not hold).
    pub fn has(&self, file: &str, ctx: &str) -> bool {
        self.block(file, ctx).is_some() && (self.dir.join(file).is_file() || self.dir.join(format!("{file}.gz")).is_file())
    }

    /// An array of the directory (`<file>` or, as in the committed fixtures, `<file>.gz`).
    pub fn read(&self, file: &str) -> io::Result<Npy> {
        let path = self.dir.join(file);
        let gz = self.dir.join(format!("{file}.gz"));
        let p = if path.exists() || !gz.exists() { path } else { gz };
        let bytes = crate::tenx::read_maybe_gz(&p)?;
        Npy::parse(&bytes).map_err(|e| invalid(format!("{}: {e}", p.display())))
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
    /// files written already, or to be copied: (file, ctx, nrow, ncol, dtype, stage)
    written: Vec<(String, String, usize, usize, Dtype, String)>,
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
    /// A block of doubles, column-major `nrow x ncol`.
    pub fn add_f64(&mut self, name: &str, ctx: &str, values: Vec<f64>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.f64s.entry(name.into()).or_default().push(ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() });
    }

    /// A block of integers, column-major `nrow x ncol`.
    pub fn add_i32(&mut self, name: &str, ctx: &str, values: Vec<i32>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.i32s.entry(name.into()).or_default().push(ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() });
    }

    /// One matrix file per context (`<name>.<ctx>.npy`), as R writes `S4_positions`.
    pub fn add_big_i32(&mut self, name: &str, ctx: &str, values: Vec<i32>, nrow: usize, ncol: usize, stage: &str) {
        assert_eq!(values.len(), nrow * ncol);
        self.big.push((format!("{name}.{ctx}.npy"), ArrayBlock { ctx: ctx.into(), values, nrow, ncol, stage: stage.into() }));
    }

    /// One matrix file per context (`<name>.<ctx>.npy`), written into `dir` right away, as R
    /// writes `S2d_scale_data` and `S2d_y`: `rows` holds the matrix row-major (`nrow x ncol`),
    /// the file is R's column-major layout. Nothing of it is kept.
    #[allow(clippy::too_many_arguments)]
    pub fn write_big_f64_rows(&mut self, dir: &Path, name: &str, ctx: &str, rows: &[f64], nrow: usize, ncol: usize, stage: &str) -> io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let file = format!("{name}.{ctx}.npy");
        npy::write_f64_rows(&dir.join(&file), rows, nrow, ncol)?;
        self.written.push((file, ctx.into(), nrow, ncol, Dtype::F8, stage.into()));
        Ok(())
    }

    /// The files the arrays go to.
    pub fn names(&self) -> Vec<String> {
        let mut v: Vec<String> = self.f64s.keys().chain(self.i32s.keys()).map(|n| format!("{n}.npy")).collect();
        v.extend(self.big.iter().map(|(f, _)| f.clone()));
        v.extend(self.written.iter().map(|w| w.0.clone()));
        v
    }

    /// Writes every array and `index.tsv` into `dir`.
    pub fn write(&self, dir: &Path) -> io::Result<()> {
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
        for (file, ctx, nrow, ncol, dtype, stage) in &self.written {
            push_index(file, &[(ctx.as_str(), *nrow, *ncol, stage.as_str())], *dtype);
        }
        index.sort_by(|a, b| a[0].cmp(&b[0]).then(a[2].parse::<usize>().unwrap().cmp(&b[2].parse::<usize>().unwrap())));
        Table { header: ["file", "ctx", "row_start", "nrow", "ncol", "dtype", "stage"].iter().map(|s| s.to_string()).collect(), rows: index }
            .write(&dir.join("index.tsv"))
    }
}

/// One `.Random.seed` checkpoint: context, stage, label, state.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    pub ctx: String,
    pub stage: String,
    pub label: String,
    pub state: Vec<i32>,
}

/// A run's stage dumps: the arrays ([`StageWriter`]), the `.Random.seed` checkpoints in the order
/// they were taken, the per-context scalars (R's `dfi_scalar`) and the contexts, written into
/// `dir` by [`StageRecorder::finish`] in R's layout.
pub struct StageRecorder {
    pub level: DumpLevel,
    pub dir: PathBuf,
    pub writer: StageWriter,
    pub rng: Vec<Checkpoint>,
    /// per context (in order of first use): scalar name -> value
    scalars: Vec<(String, Vec<(String, Json)>)>,
    contexts: Vec<String>,
}

impl StageRecorder {
    pub fn new(dir: &Path, level: DumpLevel) -> StageRecorder {
        StageRecorder { level, dir: dir.to_path_buf(), writer: StageWriter::default(), rng: Vec::new(), scalars: Vec::new(), contexts: Vec::new() }
    }

    /// A `.Random.seed` checkpoint.
    pub fn rng(&mut self, ctx: &str, stage: &str, label: &str, state: Vec<i32>) {
        self.rng.push(Checkpoint { ctx: ctx.into(), stage: stage.into(), label: label.into(), state });
    }

    /// A context, in run order (listed in `stages.json`).
    pub fn context(&mut self, ctx: &str) {
        if !self.contexts.iter().any(|c| c == ctx) {
            self.contexts.push(ctx.into());
        }
    }

    /// A scalar of context `ctx` (R's `dfi_scalar`); a later value replaces an earlier one.
    pub fn scalar(&mut self, ctx: &str, name: &str, value: Json) {
        let i = match self.scalars.iter().position(|(c, _)| c == ctx) {
            Some(i) => i,
            None => {
                self.scalars.push((ctx.into(), Vec::new()));
                self.scalars.len() - 1
            }
        };
        let m = &mut self.scalars[i].1;
        match m.iter_mut().find(|(k, _)| k == name) {
            Some(slot) => slot.1 = value,
            None => m.push((name.into(), value)),
        }
    }

    /// Writes everything: the arrays, `index.tsv`, `rng.npy` and `rng.tsv` (one row of 626 per
    /// checkpoint, column-major as R writes the matrix), and `stages.json` (the level, the
    /// contexts, the scalars and `producer`).
    pub fn finish(self, producer: &str) -> io::Result<()> {
        let dir = self.dir.clone();
        self.writer.write(&dir)?;
        let len = self.rng.first().map_or(0, |c| c.state.len());
        if self.rng.iter().any(|c| c.state.len() != len) {
            return Err(invalid("RNG checkpoints of different lengths".into()));
        }
        let n = self.rng.len();
        let mut m = Vec::with_capacity(n * len);
        for j in 0..len {
            m.extend(self.rng.iter().map(|c| c.state[j]));
        }
        npy::write_i32(&dir.join("rng.npy"), &m, &[n, len])?;
        Table {
            header: ["row", "ctx", "stage", "label"].iter().map(|s| s.to_string()).collect(),
            rows: self.rng.iter().enumerate().map(|(i, c)| vec![i.to_string(), c.ctx.clone(), c.stage.clone(), c.label.clone()]).collect(),
        }
        .write(&dir.join("rng.tsv"))?;
        let scalars = Json::Object(self.scalars.iter().map(|(c, m)| (c.clone(), Json::Object(m.clone()))).collect());
        let json = Json::obj(vec![
            ("level", Json::Str(self.level.name().into())),
            ("contexts", Json::strs(&self.contexts)),
            ("scalars", scalars),
            ("producer", Json::Str(producer.into())),
        ]);
        std::fs::write(dir.join("stages.json"), crate::json::to_file_text(&json))
    }
}

/// `.Random.seed` checkpoints of a stage directory (`rng.npy`, `rng.tsv`).
pub struct RngCheckpoints {
    states: Npy,
    labels: Table,
}

impl RngCheckpoints {
    pub fn open(stages: &StageDir) -> io::Result<RngCheckpoints> {
        let states = stages.read("rng.npy")?;
        let path = stages.dir.join("rng.tsv");
        let text = String::from_utf8(crate::tenx::read_maybe_gz(&if path.exists() { path.clone() } else { stages.dir.join("rng.tsv.gz") })?)
            .map_err(|_| invalid(format!("{}: not text", path.display())))?;
        let labels = Table::parse(&text).ok_or_else(|| invalid(format!("{}: empty", path.display())))?;
        if labels.header != ["row", "ctx", "stage", "label"] {
            return Err(invalid(format!("{}: unexpected header {:?}", path.display(), labels.header)));
        }
        Ok(RngCheckpoints { states, labels })
    }

    /// Every checkpoint in order: (ctx, stage, label, state).
    pub fn all(&self) -> Vec<(String, String, String, Vec<i32>)> {
        self.labels
            .rows
            .iter()
            .filter_map(|r| {
                let row: usize = r[0].parse().ok()?;
                Some((r[1].clone(), r[2].clone(), r[3].clone(), (0..self.states.ncol()).map(|c| self.states.get_i32(row, c)).collect()))
            })
            .collect()
    }

    /// The state labelled `label` in context `ctx` (the first one, if several).
    pub fn state(&self, ctx: &str, label: &str) -> Option<Vec<i32>> {
        let row: usize = self.labels.rows.iter().find(|r| r[1] == ctx && r[3] == label)?[0].parse().ok()?;
        Some((0..self.states.ncol()).map(|c| self.states.get_i32(row, c)).collect())
    }

    /// Every label of context `ctx`, in order.
    pub fn labels(&self, ctx: &str) -> Vec<String> {
        self.labels.rows.iter().filter(|r| r[1] == ctx).map(|r| r[3].clone()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn written_dumps_read_back() {
        let dir = crate::test_dir("stages");
        let _ = std::fs::remove_dir_all(&dir);
        let mut r = StageRecorder::new(&dir, DumpLevel::Light);
        r.context("chain");
        r.writer.add_f64("A", "chain", vec![1.0, 2.0, 3.0, 4.0], 2, 2, "S1");
        r.writer.add_f64("A", "other", vec![5.0, 6.0], 1, 2, "S1");
        r.writer.add_i32("B", "chain", vec![7, 8, 9], 3, 1, "S2");
        r.rng("chain", "S1", "entry", vec![10403, 624, 1, 2]);
        r.rng("chain", "S2", "exit", vec![10403, 1, 3, 4]);
        r.scalar("chain", "x", Json::Num(0.5));
        r.scalar("chain", "x", Json::Num(1.5));
        r.finish("test").unwrap();
        let s = StageDir::open(&dir).unwrap();
        let b = s.block("A.npy", "other").unwrap().clone();
        assert_eq!((b.row_start, b.nrow, b.ncol, b.stage.as_str()), (2, 1, 2, "S1"));
        let a = s.read("A.npy").unwrap();
        assert_eq!((a.get(0, 1), a.get(2, 1)), (3.0, 6.0));
        assert!(s.has("B.npy", "chain") && !s.has("B.npy", "other"));
        let rng = RngCheckpoints::open(&s).unwrap();
        assert_eq!(rng.state("chain", "exit"), Some(vec![10403, 1, 3, 4]));
        assert_eq!(rng.labels("chain"), vec!["entry", "exit"]);
        assert_eq!(rng.all().len(), 2);
        let j = Json::parse(&std::fs::read_to_string(dir.join("stages.json")).unwrap()).unwrap();
        assert_eq!(j.get("scalars").and_then(|s| s.get("chain")).and_then(|c| c.get("x")).and_then(Json::as_f64), Some(1.5));
    }
}
