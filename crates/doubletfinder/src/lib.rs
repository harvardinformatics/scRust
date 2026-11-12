//! DoubletFinder 2.0.6 (commit 1b244d8), reimplemented without Seurat objects
//! (docs/PLAN.md Part B; milestones M3 to M7).
//!
//! It reproduces the workflow's `doubletfinder.R` call sequence ([`run`]): `paramSweep` over
//! the pN x pK grid, `summarizeSweep` and `find.pK` (the bimodality coefficient), nExp with the
//! homotypic adjustment, and the final `doubletFinder` classification. It writes pANN and the
//! classification for every cell, the sweep tables, a `summary.json` that includes the column
//! names R would write, and the singlets' counts as a new 10x directory.
//!
//! - [`params`]: the pK and pN grids, and the values that pass through strings;
//! - [`pann`]: stages S4 and S5 (neighbours through the [`knn`] crate, pANN) and S7 (calls);
//! - [`summary`]: stage S6 (`summarizeSweep`, `find.pK`) and `modelHomotypic`;
//! - [`run`]: the call sequence over an [`run::EmbeddingSource`];
//! - [`sct`]: SCTransform as DoubletFinder calls it (stage S2): the doublet-augmented input
//!   and the comparison with R's stage dumps;
//! - [`pca`]: RunPCA as DoubletFinder calls it (stage S3) and the comparison with R's dumps,
//!   teacher-forced and chained after Rust's SCTransform;
//! - [`chain`]: the whole pipeline in Rust (M6): doublets drawn from R's random-number stream,
//!   SCTransform and PCA, as an [`run::EmbeddingSource`];
//! - [`input`]: what the command reads: the counts, the cells' cluster labels and the factor's
//!   levels;
//! - [`stages`]: stage dumps in the instrumented R runs' layout, read (teacher forcing) and
//!   written; [`record`]: a Rust run's dumps of every stage, S1 to S7, with the RNG
//!   checkpoints and scalars;
//! - [`singlets`]: the singlets' counts as a new 10x directory, the step's count-level result
//!   (M6b);
//! - [`validate`]: comparison of a Rust run with an R reference run, and with the saved
//!   workflow output.
//!
//! Milestone M3 ran it on R's own embeddings (teacher forcing, plan A5b, still available through
//! [`stages::StageSource`]); from M6 on, [`chain::RustSource`] computes them from the counts.

pub mod chain;
pub mod input;
pub mod pann;
pub mod params;
pub mod pca;
pub mod record;
pub mod run;
pub mod sct;
pub mod singlets;
pub mod stages;
pub mod summary;
pub mod validate;

use std::fmt;

/// The original implementation this crate reproduces (docs/PLAN.md B1).
pub const REFERENCE: &str = "DoubletFinder 2.0.6 (commit 1b244d8)";

/// Why a run stops.
#[derive(Debug)]
pub enum DfError {
    /// Inconsistent or invalid input, or a condition under which R stops.
    Input(String),
    /// `bkde` stops for one sweep result (as R does), e.g. when every pANN is identical.
    Bkde { name: String, source: kde::BkdeError },
    Io(std::io::Error),
}

impl fmt::Display for DfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DfError::Input(m) => f.write_str(m),
            DfError::Bkde { name, source } => write!(f, "summarizeSweep, {name}: {source}"),
            DfError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DfError {}

impl From<std::io::Error> for DfError {
    fn from(e: std::io::Error) -> Self {
        DfError::Io(e)
    }
}

/// `items.map(f)` on `threads` threads; results in item order, whatever the scheduling.
pub(crate) fn par_map<T: Sync, R: Send>(items: &[T], threads: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let threads = threads.max(1).min(items.len().max(1));
    if threads == 1 {
        return items.iter().map(f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<R>>> = items.iter().map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= items.len() {
                    break;
                }
                *slots[i].lock().unwrap() = Some(f(&items[i]));
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
}
