//! Stage S6: `summarizeSweep(GT = FALSE)` and `find.pK`, and the homotypic adjustment of nExp
//! (`modelHomotypic`).

use crate::params::Labelled;
use crate::DfError;
use kde::{bimodality_coefficient, bkde, density_function, Bkde};
use rcompat::arith::{round, seq_length_out};
use rcompat::{mean, sd, sum};

/// One (pN, pK) of the sweep: its bimodality coefficient, and the intermediate values the
/// instrumented R run dumps (when kept).
#[derive(Debug, Clone)]
pub struct SweepStat {
    pub pn: Labelled,
    pub pk: Labelled,
    /// `BCreal`: the bimodality coefficient of the pANN density.
    pub bc_real: f64,
    /// `bkde(pANN)`, if kept.
    pub bkde: Option<Bkde>,
    /// `gkde(x)` with `x = seq(min(pANN), max(pANN), length.out = n)`, if kept.
    pub gkde: Option<Vec<f64>>,
}

/// `summarizeSweep`'s statistic for one pANN vector: `approxfun(bkde(pANN, kernel = "normal"))`
/// evaluated at `seq(min(pANN), max(pANN), length.out = length(pANN))`, then
/// `bimodality_coefficient` of those values. Returns BCreal, the estimate and the values.
pub fn bc_of_pann(pann: &[f64]) -> Result<(f64, Bkde, Vec<f64>), kde::BkdeError> {
    let est = bkde(pann)?;
    let fun = density_function(&est).map_err(|_| kde::BkdeError::NotFinite)?;
    let (lo, hi) = pann.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
    let xs = seq_length_out(lo, hi, pann.len() as f64).map_err(|_| kde::BkdeError::NotFinite)?;
    let g = fun.eval_all(&xs);
    Ok((bimodality_coefficient(&g), est, g))
}

/// `summarizeSweep(sweep.list, GT = FALSE)`: BCreal for every (pN, pK), pN-major as the sweep
/// list. `pann[j][k]` is the pANN vector of pN `pn[j]` and pK `pk[k]`. pN and pK are carried as
/// the strings in the list names and the values R parses back from them.
pub fn summarize_sweep(pann: &[Vec<Vec<f64>>], pn: &[Labelled], pk: &[Labelled], keep: bool, threads: usize) -> Result<Vec<SweepStat>, DfError> {
    let items: Vec<(usize, usize)> = (0..pn.len()).flat_map(|j| (0..pk.len()).map(move |k| (j, k))).collect();
    let compute = |&(j, k): &(usize, usize)| -> Result<SweepStat, DfError> {
        let (bc, est, g) = bc_of_pann(&pann[j][k]).map_err(|e| DfError::Bkde {
            name: crate::params::sweep_name(&pn[j], &pk[k]),
            source: e,
        })?;
        Ok(SweepStat {
            pn: pn[j].clone(),
            pk: pk[k].clone(),
            bc_real: bc,
            bkde: keep.then_some(est),
            gkde: keep.then_some(g),
        })
    };
    // independent items, results in item order: the same output for any thread count
    crate::par_map(&items, threads, compute).into_iter().collect()
}

/// One row of `find.pK`'s result (`bcmvn`).
#[derive(Debug, Clone, PartialEq)]
pub struct BcmvnRow {
    /// `ParamID`: 1, 2, ...
    pub param_id: usize,
    /// The pK factor level: `as.character()` of the parsed pK.
    pub pk_level: String,
    /// `mean(BCreal)` over pN.
    pub mean_bc: f64,
    /// `sd(BCreal)^2`: the square of the standard deviation, not `var()`.
    pub var_bc: f64,
    /// `mean(BCreal) / (sd(BCreal)^2)`.
    pub bc_metric: f64,
}

/// `find.pK(sweep.stats)` without ground truth: for each pK level, in order of first
/// appearance, the mean and squared standard deviation of BCreal over pN, and their ratio.
pub fn find_pk(stats: &[SweepStat]) -> Vec<BcmvnRow> {
    let mut levels: Vec<String> = Vec::new();
    for s in stats {
        let l = s.pk.level();
        if !levels.contains(&l) {
            levels.push(l);
        }
    }
    levels
        .into_iter()
        .enumerate()
        .map(|(i, level)| {
            let bc: Vec<f64> = stats.iter().filter(|s| s.pk.level() == level).map(|s| s.bc_real).collect();
            let m = mean(&bc, false);
            let s = sd(&bc, false);
            let s2 = s * s; // sd(...)^2
            BcmvnRow { param_id: i + 1, pk_level: level, mean_bc: m, var_bc: s2, bc_metric: mean(&bc, false) / s2 }
        })
        .collect()
}

/// `which.max(x)`: the first index of the maximum, ignoring NaN and NA; None if all are.
pub fn which_max(x: &[f64]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, &v) in x.iter().enumerate() {
        if !v.is_nan() && best.is_none_or(|b| v > x[b]) {
            best = Some(i);
        }
    }
    best
}

/// `modelHomotypic(annotations)`: `sum((table(annotations)/length(annotations))^2)`, the table
/// over the factor levels in level order (empty levels add 0), labels outside the levels being
/// NA (not tabulated, but counted in the length). R's long-double `sum`.
pub fn model_homotypic(labels: &[String], levels: &[String]) -> f64 {
    let n = labels.len() as f64;
    let counts: Vec<usize> = levels.iter().map(|l| labels.iter().filter(|x| *x == l).count()).collect();
    let sq: Vec<f64> = counts
        .iter()
        .map(|&c| {
            let f = c as f64 / n;
            f * f
        })
        .collect();
    sum(&sq, false)
}

/// The workflow's expected doublet numbers: `nExp_poi = round(rate * n)` and
/// `nExp_poi.adj = round(nExp_poi * (1 - homotypic.prop))`.
pub fn n_exp(n: usize, rate: f64, homotypic_prop: f64) -> (f64, f64) {
    let poi = round(rate * n as f64);
    (poi, round(poi * (1.0 - homotypic_prop)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn which_max_like_r() {
        assert_eq!(which_max(&[1.0, f64::NAN, 3.0, 3.0, 2.0]), Some(2));
        assert_eq!(which_max(&[f64::NAN, f64::NAN]), None);
        assert_eq!(which_max(&[f64::NAN, f64::INFINITY, 5.0]), Some(1));
    }

    #[test]
    fn homotypic_and_nexp() {
        let labels: Vec<String> = ["0", "1", "0", "2", "0", "1"].iter().map(|s| s.to_string()).collect();
        let levels: Vec<String> = ["0", "1", "2", "3"].iter().map(|s| s.to_string()).collect();
        let p = model_homotypic(&labels, &levels);
        assert!((p - (0.25 + 1.0 / 9.0 + 1.0 / 36.0)).abs() < 1e-15);
        // the tier-T object soupx_tenx_test: 300 cells in clusters of 165, 121 and 14
        let l: Vec<String> = std::iter::repeat_n("0", 165).chain(std::iter::repeat_n("1", 121)).chain(std::iter::repeat_n("2", 14)).map(String::from).collect();
        let lv: Vec<String> = ["0", "1", "2"].iter().map(|s| s.to_string()).collect();
        let prop = model_homotypic(&l, &lv);
        assert_eq!(scio::fmt17(prop), "0.46735555555555558"); // R's homotypic_prop
        assert_eq!(n_exp(300, 0.15, prop), (45.0, 24.0));
    }
}
