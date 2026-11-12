//! Stage S2d: the Pearson residuals `vst` computes with `residual_type = "pearson"` (Seurat 5.5.1's
//! default `sct.method`), and what `SCTransform.default` makes of them:
//!
//! 1. `min_var = (median(umi@x) / 5)^2` (`min_variance = "umi_median"`), over the stored entries
//!    of the kept genes' counts;
//! 2. for each bin of `bin_size` kept genes, `mu = exp(tcrossprod(model_pars_fit[, -1],
//!    model.matrix(~ log_umi)))` (OpenBLAS `dgemm`, as R calls it), the residuals
//!    `(y - mu) / sqrt(max(mu + mu^2 / theta, min_var))`, clipped to `+-sqrt(ncells)`;
//! 3. per gene, `residual_mean = rowMeans(y)` (base R: long double) and `residual_variance =
//!    matrixStats::rowVars(y)` (double, refined);
//! 4. the variable features: `sort(residual_variance, decreasing = TRUE)` (stable; NA dropped),
//!    the first `variable_features_n`;
//! 5. `scale.data`: their residuals clipped to `+-sqrt(ncells / 30)`, then centred by
//!    `ScaleData` (`FastRowScale`: `mat - rowMeans2(mat, na.rm = TRUE)`).
//!
//! R holds the full residual matrix (genes x cells). Here the bins' products are made one at a
//! time, each by one `tcrossprod` with R's arguments into one of two reused buffers, and each bin's
//! genes are shared among the threads (eight at a time, or one at a time above 16,384 cells:
//! [`gene_group`]); while they turn one bin's product into residuals, one thread makes the next
//! bin's. Only the rows of the current best `variable_features_n` genes are kept, as slots of one
//! buffer that becomes `scale.data` in place (its rows reordered, clipped and centred), so memory
//! stays near `(variable_features_n + 2 * bin_size + 2 * threads) * ncells` doubles: two products,
//! the kept rows, and a residual and a counts row per thread. Results do not depend on the number
//! of threads.

use crate::matrix::GeneRows;
use crate::params::{MinVariance, SctParams};
use crate::setup::bins;
use crate::SctError;
use rcompat::blas::tcrossprod_into;
use rcompat::libm::exp;
use rcompat::matrixstats::{mean2_na_rm, var_refined};
use rcompat::F80;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, VecDeque};
use std::ops::Range;
use std::sync::{Arc, Condvar, Mutex, MutexGuard};

/// Bins' products held at once: the bin being turned into residuals, and the next one.
const PRODUCT_BUFFERS: usize = 2;

/// Genes a thread takes at once, for `ncells` cells. A product is column-major, so eight
/// consecutive genes share their cache lines (64 bytes, two lines per cell where the eight
/// straddle a boundary). Up to 16,384 cells those lines fit in a core's own cache (2 MB of L2 on
/// Sapphire Rapids), so one thread takes eight genes and finds the lines there for seven of them;
/// with more cells they would leave it between genes, and the threads, taking neighbouring genes
/// one at a time, share them through the L3 instead. Only the speed depends on it.
fn gene_group(ncells: usize) -> usize {
    if ncells <= 16_384 {
        8
    } else {
        1
    }
}

/// What S2d reads.
#[derive(Debug, Clone)]
pub struct ResidualInput<'a> {
    /// per kept gene: theta, `(Intercept)`, `log_umi` (`model_pars_fit`)
    pub model_pars_fit: &'a [[f64; 3]],
    /// `log10(colSums(umi))` per cell (`cell_attr$log_umi`)
    pub log_umi: &'a [f64],
}

/// The result of S2d.
#[derive(Debug, Clone)]
pub struct Residuals {
    pub ncells: usize,
    /// `min_var`, and the clips `sqrt(ncells)` and `sqrt(ncells / 30)` (upper bounds)
    pub min_var: f64,
    pub res_clip: (f64, f64),
    pub clip: (f64, f64),
    /// per kept gene
    pub residual_mean: Vec<f64>,
    pub residual_variance: Vec<f64>,
    /// the variable features, positions in the kept genes, highest residual variance first
    pub top_features: Vec<usize>,
    /// `scale.data`: row `k` (cells in order) belongs to `top_features[k]`; row-major
    pub scale_data: Vec<f64>,
    /// the clipped residuals of every kept gene (row-major), when asked for
    pub y: Option<Vec<f64>>,
}

impl Residuals {
    /// Row `k` of `scale.data`.
    pub fn scale_row(&self, k: usize) -> &[f64] {
        &self.scale_data[k * self.ncells..(k + 1) * self.ncells]
    }
}

/// `min_var` for the kept genes' counts.
pub fn min_var(rows: &GeneRows, params: &SctParams) -> f64 {
    match params.min_variance {
        MinVariance::UmiMedian => {
            let m = rows.stored_median() / 5.0;
            m * m // R's x^2 is x * x
        }
        MinVariance::Value(v) => v,
    }
}

/// A candidate variable feature: ordered by residual variance, ties by gene position (the earlier
/// gene first, as R's stable sort keeps it).
#[derive(Debug, Clone, Copy)]
struct Cand {
    var: f64,
    gene: usize,
}

impl Cand {
    /// "better" = earlier in R's decreasing order
    fn better(&self, other: &Cand) -> Ordering {
        match self.var.partial_cmp(&other.var).expect("NaN residual variance in the ranking") {
            Ordering::Equal => other.gene.cmp(&self.gene),
            o => o,
        }
    }
}

impl PartialEq for Cand {
    fn eq(&self, o: &Self) -> bool {
        self.better(o) == Ordering::Equal
    }
}
impl Eq for Cand {}
impl PartialOrd for Cand {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Cand {
    /// reversed, so that a max-heap keeps the worst candidate on top
    fn cmp(&self, o: &Self) -> Ordering {
        o.better(self)
    }
}

/// A bin whose product is ready; its genes are handed out `Run::group` at a time.
struct ActiveBin {
    bin: usize,
    eta: Arc<Vec<f64>>,
    /// next gene to hand out, genes finished (positions in the bin)
    next: usize,
    done: usize,
}

/// A row of the slot buffer: its index, and the row.
type Slot<'s> = (usize, &'s mut [f64]);

/// What the threads share, behind one lock.
struct State<'s> {
    /// the next bin to make the product of; whether a thread is making one
    next_product: usize,
    producing: bool,
    /// product buffers allocated so far (at most `PRODUCT_BUFFERS`), and the unused ones
    buffers: usize,
    free: Vec<Vec<f64>>,
    /// bins with a product, in bin order
    active: VecDeque<ActiveBin>,
    finished: usize,
    /// a thread panicked: the others stop
    failed: bool,
    residual_mean: Vec<f64>,
    residual_variance: Vec<f64>,
    /// the best candidates so far (worst on top) and their residual rows, slots of one buffer;
    /// the slots that hold no kept row and no thread's working row
    heap: BinaryHeap<Cand>,
    kept: HashMap<usize, Slot<'s>>,
    free_slots: Vec<Slot<'s>>,
    all_y: Option<Vec<f64>>,
}

/// One S2d computation: what every thread reads, and the shared state.
struct Run<'a, 's> {
    rows: &'a GeneRows,
    model_pars_fit: &'a [[f64; 3]],
    /// `model.matrix(~ log_umi, cell_attr)`: (Intercept), log_umi; column-major ncells x 2
    xreg: Vec<f64>,
    bins: Vec<Range<usize>>,
    /// genes a thread takes at once ([`gene_group`])
    group: usize,
    min_var: f64,
    res_clip: (f64, f64),
    /// keep the best `k_top` genes at or above `threshold`
    k_top: usize,
    threshold: f64,
    state: Mutex<State<'s>>,
    cv: Condvar,
}

/// Stops the other threads when the one holding it panics, instead of leaving them waiting.
struct StopOnPanic<'r, 'a, 's>(&'r Run<'a, 's>);

impl Drop for StopOnPanic<'_, '_, '_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.lock().failed = true;
            self.0.cv.notify_all();
        }
    }
}

impl<'s> Run<'_, 's> {
    fn lock(&self) -> MutexGuard<'_, State<'s>> {
        match self.state.lock() {
            Ok(st) => st,
            Err(poisoned) => {
                let mut st = poisoned.into_inner();
                st.failed = true;
                st
            }
        }
    }

    fn wait<'g>(&self, st: MutexGuard<'g, State<'s>>) -> MutexGuard<'g, State<'s>> {
        match self.cv.wait(st) {
            Ok(st) => st,
            Err(poisoned) => {
                let mut st = poisoned.into_inner();
                st.failed = true;
                st
            }
        }
    }

    /// One thread's loop: make the next bin's product when a buffer is free, otherwise compute
    /// the residuals of the earliest bin's next genes, otherwise wait.
    fn work(&self) {
        let _stop = StopOnPanic(self);
        let nc = self.rows.ncol;
        let mut counts = vec![0.0f64; nc];
        let mut st = self.lock();
        // the working row: a slot of its own (there is one per thread beyond the kept rows)
        let mut row: Slot<'s> = st.free_slots.pop().expect("a slot per thread");
        loop {
            if st.failed {
                return;
            }
            if !st.producing && st.next_product < self.bins.len() && (st.buffers < PRODUCT_BUFFERS || !st.free.is_empty()) {
                let b = st.next_product;
                st.next_product += 1;
                st.producing = true;
                let mut eta = match st.free.pop() {
                    Some(buf) => buf,
                    None => {
                        st.buffers += 1;
                        Vec::new()
                    }
                };
                drop(st);
                self.product(b, &mut eta);
                st = self.lock();
                st.producing = false;
                st.active.push_back(ActiveBin { bin: b, eta: Arc::new(eta), next: 0, done: 0 });
                self.cv.notify_all();
                continue;
            }
            let job = st.active.iter_mut().find(|a| a.next < self.bins[a.bin].len()).map(|a| {
                let first = a.next;
                a.next = (first + self.group).min(self.bins[a.bin].len());
                (a.bin, first..a.next, Arc::clone(&a.eta))
            });
            if let Some((b, genes, eta)) = job {
                drop(st);
                let nb = self.bins[b].len();
                for r in genes.clone() {
                    let g = self.bins[b].start + r;
                    let (mean, var) = self.gene(&eta, nb, r, g, &mut counts, row.1);
                    self.record(&mut self.lock(), g, mean, var, &mut row);
                }
                drop(eta);
                st = self.lock();
                let i = st.active.iter().position(|a| a.bin == b).expect("the genes' bin is active");
                st.active[i].done += genes.len();
                if st.active[i].done == self.bins[b].len() {
                    // every gene of the bin is done, and each thread dropped its handle on the
                    // product before reporting: the buffer is free for the next product
                    let a = st.active.remove(i).expect("active bin");
                    let buf = Arc::try_unwrap(a.eta).unwrap_or_else(|_| panic!("bin {b}: product still in use"));
                    st.free.push(buf);
                    st.finished += 1;
                    self.cv.notify_all();
                }
                continue;
            }
            if st.finished == self.bins.len() {
                return;
            }
            st = self.wait(st);
        }
    }

    /// The product of bin `b`, `tcrossprod(model_pars_fit[genes_bin, -1], model.matrix(~ log_umi))`:
    /// genes x cells, column-major, into `eta`.
    fn product(&self, b: usize, eta: &mut Vec<f64>) {
        let genes = self.bins[b].clone();
        let nb = genes.len();
        let nc = self.rows.ncol;
        // model_pars_final[genes_bin, -1]: (Intercept) and log_umi, column-major nb x 2
        let mut p = vec![0.0f64; nb * 2];
        for (r, g) in genes.enumerate() {
            p[r] = self.model_pars_fit[g][1];
            p[r + nb] = self.model_pars_fit[g][2];
        }
        eta.resize(nb * nc, 0.0);
        tcrossprod_into(&p, nb, 2, &self.xreg, nc, 2, eta);
    }

    /// The clipped residuals of gene `g`, row `r` of its bin's product `eta` (`nb` rows), into
    /// `row`; returns their `rowMeans` and `rowVars`.
    fn gene(&self, eta: &[f64], nb: usize, r: usize, g: usize, counts: &mut [f64], row: &mut [f64]) -> (f64, f64) {
        let nc = self.rows.ncol;
        let (min_var, clip) = (self.min_var, self.res_clip);
        self.rows.dense_row(g, counts);
        let theta = self.model_pars_fit[g][0];
        for c in 0..nc {
            let mu = exp(eta[r + c * nb]);
            // pearson_residual(): model_var <- mu + mu^2 / theta; floored at min_var
            let mut model_var = mu + mu * mu / theta;
            if model_var < min_var {
                model_var = min_var;
            }
            let mut v = (counts[c] - mu) / model_var.sqrt();
            // clip_matrix_values(y, res_clip_range)
            if v < clip.0 {
                v = clip.0;
            }
            if v > clip.1 {
                v = clip.1;
            }
            row[c] = v;
        }
        // rowMeans(y): long-double sum in column order, divided by ncol in long double
        let mut s = F80::ZERO;
        for &v in row.iter() {
            s = s.add_f64(v);
        }
        ((s / F80::from_u64(nc as u64)).to_f64(), var_refined(row))
    }

    /// Stores gene `g`'s statistics, and its row if it is among the best so far: the row's slot
    /// moves into `kept`, and `row` becomes the slot of the row it displaces (or a free slot).
    fn record(&self, st: &mut State<'s>, g: usize, mean: f64, var: f64, row: &mut Slot<'s>) {
        let nc = self.rows.ncol;
        st.residual_mean[g] = mean;
        st.residual_variance[g] = var;
        if let Some(all) = st.all_y.as_mut() {
            all[g * nc..(g + 1) * nc].copy_from_slice(row.1);
        }
        if var.is_nan() || self.k_top == 0 || !(var >= self.threshold) {
            return; // sort() drops NA
        }
        let cand = Cand { var, gene: g };
        if st.heap.len() < self.k_top {
            // fewer kept rows than k_top: a slot is free (the buffer has k_top plus one per thread)
            let fresh = st.free_slots.pop().expect("a free slot while fewer rows are kept");
            st.heap.push(cand);
            st.kept.insert(g, std::mem::replace(row, fresh));
        } else if cand.better(st.heap.peek().expect("full heap")) == Ordering::Greater {
            let worst = st.heap.pop().expect("full heap");
            let old = st.kept.remove(&worst.gene).expect("a kept row per candidate");
            st.heap.push(cand);
            st.kept.insert(g, std::mem::replace(row, old));
        }
    }
}

/// Stage S2d. `rows` are the kept genes' counts, `threads` the threads sharing each bin's genes.
pub fn residuals(rows: &GeneRows, input: &ResidualInput, params: &SctParams, keep_y: bool, threads: usize) -> Result<Residuals, SctError> {
    residuals_in_groups(rows, input, params, keep_y, threads, gene_group(rows.ncol))
}

/// [`residuals`], the threads taking `group` genes at a time.
fn residuals_in_groups(rows: &GeneRows, input: &ResidualInput, params: &SctParams, keep_y: bool, threads: usize, group: usize) -> Result<Residuals, SctError> {
    let ng = rows.nrow();
    let nc = rows.ncol;
    if input.model_pars_fit.len() != ng || input.log_umi.len() != nc {
        return Err(SctError::Input("residuals: model_pars_fit or log_umi does not match the counts".into()));
    }
    let min_var = min_var(rows, params);
    let res_clip = params.res_clip(nc);
    let clip = params.sct_clip(nc);
    // model.matrix(~ log_umi, cell_attr): (Intercept), log_umi; column-major ncells x 2
    let mut xreg = vec![1.0f64; nc * 2];
    xreg[nc..].copy_from_slice(input.log_umi);
    // variable.features.n, or every gene at or above variable.features.rv.th
    let (k_top, threshold) = match params.variable_features_n {
        Some(n) => (n.min(ng), f64::NEG_INFINITY),
        None => (usize::MAX, params.variable_features_rv_th),
    };
    // the threads share each bin's genes: more threads than genes would idle
    let threads = threads.max(1).min(ng.max(1));
    // One buffer of rows: at most k_top kept rows (every gene with the threshold) and one working
    // row per thread. Its pages are touched as rows are written; it becomes scale.data in place.
    let n_slots = k_top.min(ng) + threads;
    let mut slab = vec![0.0f64; n_slots * nc];
    let (residual_mean, residual_variance, top_features, slot_of_top, all_y) = {
        let mut slots: Vec<Slot> = Vec::with_capacity(n_slots);
        let mut rest: &mut [f64] = &mut slab;
        for i in 0..n_slots {
            let (row, tail) = std::mem::take(&mut rest).split_at_mut(nc);
            slots.push((i, row));
            rest = tail;
        }
        // popped from the end: the lowest slots first
        slots.reverse();
        let run = Run {
            rows,
            model_pars_fit: input.model_pars_fit,
            xreg,
            bins: bins(ng, params.bin_size),
            group: group.max(1),
            min_var,
            res_clip,
            k_top,
            threshold,
            state: Mutex::new(State {
                next_product: 0,
                producing: false,
                buffers: 0,
                free: Vec::new(),
                active: VecDeque::new(),
                finished: 0,
                failed: false,
                residual_mean: vec![0.0f64; ng],
                residual_variance: vec![0.0f64; ng],
                heap: BinaryHeap::with_capacity(k_top.min(ng) + 1),
                kept: HashMap::with_capacity(k_top.min(ng) + 1),
                free_slots: slots,
                all_y: if keep_y { Some(vec![0.0f64; ng * nc]) } else { None },
            }),
            cv: Condvar::new(),
        };
        std::thread::scope(|scope| {
            for _ in 1..threads {
                scope.spawn(|| run.work());
            }
            run.work();
        });
        // (a thread's panic has been propagated by the scope)
        let State { residual_mean, residual_variance, heap, kept, all_y, .. } = run.state.into_inner().unwrap_or_else(|e| e.into_inner());
        // the variable features in R's order, and the slot of each one's row
        let mut top: Vec<Cand> = heap.into_vec();
        top.sort_by(|a, b| b.better(a));
        let top_features: Vec<usize> = top.iter().map(|c| c.gene).collect();
        let slot_of_top: Vec<usize> = top_features.iter().map(|g| kept.get(g).expect("a kept row per variable feature").0).collect();
        (residual_mean, residual_variance, top_features, slot_of_top, all_y)
    };
    // scale.data in place: the kept rows in R's order, then clipped to +-sqrt(ncells / 30) and
    // (do.center) centred by rowMeans2(na.rm = TRUE), as ScaleData does
    gather_rows(&mut slab, nc, &slot_of_top);
    slab.truncate(top_features.len() * nc);
    slab.shrink_to_fit();
    let mut scale_data = slab;
    for k in 0..top_features.len() {
        let row = &mut scale_data[k * nc..(k + 1) * nc];
        for d in row.iter_mut() {
            let mut v = *d;
            if v < clip.0 {
                v = clip.0;
            }
            if v > clip.1 {
                v = clip.1;
            }
            *d = v;
        }
        if params.do_center {
            let rm = mean2_na_rm(row);
            for d in row.iter_mut() {
                *d -= rm;
            }
        }
    }
    Ok(Residuals { ncells: nc, min_var, res_clip, clip, residual_mean, residual_variance, top_features, scale_data, y: all_y })
}

/// Reorders the rows of `m` (each `nc` long) in place so that row `k` holds what row `src[k]`
/// held, for every `k < src.len()`; the rows after them are left in no particular state. `src`
/// must not repeat a row. Rows whose old content no target needs start chains, which end at a
/// row beyond the targets; what remains are cycles, rotated through one spare row.
fn gather_rows(m: &mut [f64], nc: usize, src: &[usize]) {
    let k = src.len();
    if nc == 0 || k == 0 {
        return;
    }
    let n = m.len() / nc;
    let mut needed = vec![false; n];
    for &s in src {
        assert!(s < n && !needed[s], "gather_rows: a source row out of range or repeated");
        needed[s] = true;
    }
    let mut done = vec![false; k];
    let copy = |m: &mut [f64], from: usize, to: usize| {
        if from != to {
            m.copy_within(from * nc..(from + 1) * nc, to * nc);
        }
    };
    for t in 0..k {
        if needed[t] || done[t] {
            continue;
        }
        // row t's old content is not needed: fill it, then the row it read from, and so on
        let mut cur = t;
        loop {
            let s = src[cur];
            copy(m, s, cur);
            done[cur] = true;
            if s >= k || done[s] {
                break;
            }
            cur = s;
        }
    }
    let mut spare: Vec<f64> = Vec::new();
    for t in 0..k {
        if done[t] {
            continue;
        }
        // a cycle through row t
        if spare.len() != nc {
            spare = vec![0.0f64; nc];
        }
        spare.copy_from_slice(&m[t * nc..(t + 1) * nc]);
        let mut cur = t;
        loop {
            let s = src[cur];
            done[cur] = true;
            if s == t {
                m[cur * nc..(cur + 1) * nc].copy_from_slice(&spare);
                break;
            }
            copy(m, s, cur);
            cur = s;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcompat::blas::tcrossprod;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn unif(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 31) as f64
        }
    }

    /// `ng` genes x `nc` cells, about a quarter of the entries stored; gene `dup.1` repeats gene
    /// `dup.0` (with its model parameters), so that two residual variances tie across bins.
    fn data(ng: usize, nc: usize, dup: (usize, usize), seed: u64) -> (GeneRows, Vec<[f64; 3]>, Vec<f64>) {
        let mut rng = Lcg(seed);
        let mut rows = GeneRows { ncol: nc, p: vec![0], j: Vec::new(), x: Vec::new() };
        let mut pars = Vec::with_capacity(ng);
        for g in 0..ng {
            if g == dup.1 {
                let (a, b) = (rows.p[dup.0], rows.p[dup.0 + 1]);
                let (j, x) = (rows.j[a..b].to_vec(), rows.x[a..b].to_vec());
                rows.j.extend(j);
                rows.x.extend(x);
                pars.push(pars[dup.0]);
            } else {
                for c in 0..nc {
                    if rng.next() % 4 == 0 {
                        rows.j.push(c as u32);
                        rows.x.push((rng.next() % 25) as f64);
                    }
                }
                // theta = Inf (a Poisson gene) now and then
                let theta = if g % 9 == 4 { f64::INFINITY } else { 0.2 + 30.0 * rng.unif() };
                pars.push([theta, -4.0 + 3.0 * rng.unif(), std::f64::consts::LN_10]);
            }
            rows.p.push(rows.j.len());
        }
        let log_umi = (0..nc).map(|_| 1.0 + 2.0 * rng.unif()).collect();
        (rows, pars, log_umi)
    }

    /// The same computation without threads or selection on the fly: every bin's product, every
    /// gene's residual row, then R's stable decreasing sort.
    fn reference(rows: &GeneRows, input: &ResidualInput, params: &SctParams) -> Residuals {
        let (ng, nc) = (rows.nrow(), rows.ncol);
        let min_var = min_var(rows, params);
        let res_clip = params.res_clip(nc);
        let clip = params.sct_clip(nc);
        let mut xreg = vec![1.0f64; nc * 2];
        xreg[nc..].copy_from_slice(input.log_umi);
        let mut y = vec![0.0f64; ng * nc];
        let (mut mean, mut var) = (vec![0.0f64; ng], vec![0.0f64; ng]);
        let mut counts = vec![0.0f64; nc];
        for bin in bins(ng, params.bin_size) {
            let nb = bin.len();
            let mut p = vec![0.0f64; nb * 2];
            for (r, g) in bin.clone().enumerate() {
                p[r] = input.model_pars_fit[g][1];
                p[r + nb] = input.model_pars_fit[g][2];
            }
            let eta = tcrossprod(&p, nb, 2, &xreg, nc, 2);
            for (r, g) in bin.enumerate() {
                rows.dense_row(g, &mut counts);
                let theta = input.model_pars_fit[g][0];
                let row = &mut y[g * nc..(g + 1) * nc];
                for c in 0..nc {
                    let mu = exp(eta[r + c * nb]);
                    let model_var = (mu + mu * mu / theta).max(min_var);
                    row[c] = ((counts[c] - mu) / model_var.sqrt()).clamp(res_clip.0, res_clip.1);
                }
                let mut s = F80::ZERO;
                for &v in row.iter() {
                    s = s.add_f64(v);
                }
                mean[g] = (s / F80::from_u64(nc as u64)).to_f64();
                var[g] = var_refined(row);
            }
        }
        let mut order: Vec<usize> = (0..ng).filter(|&g| !var[g].is_nan()).collect();
        order.sort_by(|&a, &b| var[b].partial_cmp(&var[a]).unwrap());
        let top: Vec<usize> = match params.variable_features_n {
            Some(n) => order.into_iter().take(n).collect(),
            None => order.into_iter().filter(|&g| var[g] >= params.variable_features_rv_th).collect(),
        };
        let mut scale_data = Vec::with_capacity(top.len() * nc);
        for &g in &top {
            let mut r: Vec<f64> = y[g * nc..(g + 1) * nc].iter().map(|v| v.clamp(clip.0, clip.1)).collect();
            if params.do_center {
                let rm = mean2_na_rm(&r);
                r.iter_mut().for_each(|d| *d -= rm);
            }
            scale_data.extend(r);
        }
        Residuals { ncells: nc, min_var, res_clip, clip, residual_mean: mean, residual_variance: var, top_features: top, scale_data, y: Some(y) }
    }

    fn same(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
    }

    #[test]
    fn residuals_match_the_plain_computation_for_any_number_of_threads() {
        // 43 genes in bins of 4: ten full bins and a one-gene bin (R's dgemv route). Genes 5 and
        // 29 are the same gene at the same place of their bins, so their variances tie.
        let (rows, pars, log_umi) = data(43, 61, (5, 29), 7);
        let input = ResidualInput { model_pars_fit: &pars, log_umi: &log_umi };
        let base = SctParams { bin_size: 4, ..SctParams::default() };
        let all = reference(&rows, &input, &SctParams { variable_features_n: Some(43), ..base.clone() });
        assert_eq!(all.residual_variance[5].to_bits(), all.residual_variance[29].to_bits());
        // the tie at the cut: gene 5 is the last feature, gene 29 the first one left out
        let at_tie = all.top_features.iter().position(|&g| g == 5).unwrap() + 1;
        assert_eq!(all.top_features[at_tie], 29);
        let rv_th = {
            let mut v = all.residual_variance.clone();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v[20]
        };
        let cases = [
            SctParams { variable_features_n: Some(9), ..base.clone() },
            SctParams { variable_features_n: Some(at_tie), ..base.clone() },
            SctParams { variable_features_n: Some(43), ..base.clone() },
            SctParams { variable_features_n: Some(0), ..base.clone() },
            SctParams { variable_features_n: None, variable_features_rv_th: rv_th, ..base.clone() },
            SctParams { variable_features_n: Some(9), bin_size: 500, ..base.clone() },
            SctParams { variable_features_n: Some(9), do_center: false, ..base.clone() },
        ];
        for params in &cases {
            let want = reference(&rows, &input, params);
            // group 0: the public entry point, which picks the group for the cell count
            for (threads, keep_y, group) in [1, 2, 3, 4, 7, 64].into_iter().flat_map(|t| [false, true].into_iter().flat_map(move |k| [0, 1, 3, 8].map(|g| (t, k, g)))) {
                let got = if group == 0 {
                    residuals(&rows, &input, params, keep_y, threads).unwrap()
                } else {
                    residuals_in_groups(&rows, &input, params, keep_y, threads, group).unwrap()
                };
                let what = format!("{params:?}, {threads} threads, keep_y {keep_y}, group {group}");
                assert!(same(&got.residual_mean, &want.residual_mean), "residual_mean: {what}");
                assert!(same(&got.residual_variance, &want.residual_variance), "residual_variance: {what}");
                assert_eq!(got.top_features, want.top_features, "features: {what}");
                assert!(same(&got.scale_data, &want.scale_data), "scale.data: {what}");
                assert_eq!(got.min_var.to_bits(), want.min_var.to_bits());
                match (&got.y, keep_y) {
                    (Some(y), true) => assert!(same(y, want.y.as_ref().unwrap()), "y: {what}"),
                    (None, false) => {}
                    _ => panic!("y kept when not asked, or not kept when asked: {what}"),
                }
            }
        }
    }

    #[test]
    fn rows_gathered_in_place() {
        // every way the kept rows can sit in the slots: chains from rows nobody reads, cycles,
        // rows already in place, sources beyond the targets
        let mut rng = Lcg(3);
        for case in 0..400 {
            let n = 1 + (rng.next() % 12) as usize;
            let k = (rng.next() % (n as u64 + 1)) as usize;
            let nc = (rng.next() % 4) as usize;
            // a random injective map from k targets to n rows (a shuffled prefix)
            let mut perm: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                perm.swap(i, (rng.next() % (i as u64 + 1)) as usize);
            }
            let src = &perm[..k];
            let m: Vec<f64> = (0..n * nc).map(|v| v as f64).collect();
            let mut got = m.clone();
            gather_rows(&mut got, nc, src);
            for (t, &s) in src.iter().enumerate() {
                assert_eq!(got[t * nc..(t + 1) * nc], m[s * nc..(s + 1) * nc], "case {case}: n {n}, k {k}, src {src:?}");
            }
        }
    }

    #[test]
    fn residuals_of_degenerate_inputs() {
        // one cell: every variance is NA, so there are no variable features
        let (rows, pars, log_umi) = data(5, 1, (0, 3), 11);
        let input = ResidualInput { model_pars_fit: &pars, log_umi: &log_umi };
        let params = SctParams { bin_size: 2, ..SctParams::default() };
        for threads in [1, 3] {
            let got = residuals(&rows, &input, &params, false, threads).unwrap();
            assert!(got.residual_variance.iter().all(|v| v.is_nan()));
            assert!(got.top_features.is_empty() && got.scale_data.is_empty());
        }
        // no genes
        let rows = GeneRows { ncol: 4, p: vec![0], j: Vec::new(), x: Vec::new() };
        let got = residuals(&rows, &ResidualInput { model_pars_fit: &[], log_umi: &[2.0; 4] }, &params, true, 4).unwrap();
        assert!(got.residual_mean.is_empty() && got.top_features.is_empty() && got.y.unwrap().is_empty());
    }
}
