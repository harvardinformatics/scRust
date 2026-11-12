//! R 4.5.3's Sheather-Jones bandwidth, `bw.SJ(x)` with its defaults (`nb = 1000`, method "ste";
//! `src/library/stats/R/bandwidths.R`, C routines `bw_den`, `bw_den_binned`, `bw_phi4`,
//! `bw_phi6` in `src/library/stats/src/bandwidths.c`), solved by [`crate::zeroin::uniroot`].

use crate::arith::r_pow;
use crate::libm::{exp, pow};
use crate::quantile::iqr;
use crate::summary::sd;
use crate::zeroin::uniroot;

const DELMAX: f64 = 1000.0;
const M_1_SQRT_2PI: f64 = 0.398942280401432677939946059934;

/// `bw_den(nb, x)`: the binned pairwise distances, from all pairs (`n <= nb / 2`).
fn bw_den(nb: usize, x: &[f64]) -> Result<(f64, Vec<f64>), String> {
    let n = x.len();
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    for (i, &v) in x.iter().enumerate() {
        if !v.is_finite() {
            return Err(format!("non-finite x[{}] in bandwidth calculation", i + 1));
        }
        if v < xmin {
            xmin = v;
        }
        if v > xmax {
            xmax = v;
        }
    }
    let rang = (xmax - xmin) * 1.01;
    if rang == 0.0 {
        return Err("data are constant in bandwidth calculation".into());
    }
    let dd = rang / nb as f64;
    let mut cnt = vec![0.0f64; nb];
    // (int) truncates toward zero
    let as_int = |v: f64| v as i32 as i64;
    if xmin / dd < i32::MIN as f64 || xmax / dd > i32::MAX as f64 {
        for i in 1..n {
            let ii = as_int((x[i] - xmin) / dd);
            for j in 0..i {
                let jj = as_int((x[j] - xmin) / dd);
                cnt[(ii - jj).unsigned_abs() as usize] += 1.0;
            }
        }
    } else {
        for i in 1..n {
            let ii = as_int(x[i] / dd);
            for j in 0..i {
                let jj = as_int(x[j] / dd);
                cnt[(ii - jj).unsigned_abs() as usize] += 1.0;
            }
        }
    }
    Ok((dd, cnt))
}

/// `bw_den_binned(xxx)`: the same counts from the bin counts.
fn bw_den_binned(counts: &[i32]) -> Vec<f64> {
    let nb = counts.len();
    let mut cnt = vec![0.0f64; nb];
    for ii in 0..nb {
        let w = counts[ii] as f64;
        cnt[0] += w * (w - 1.0); // no distances to self
        for jj in 0..ii {
            cnt[ii - jj] += w * counts[jj] as f64;
        }
    }
    cnt[0] *= 0.5; // pairs in the same bin were counted twice
    cnt
}

/// `bw_pair_cnts(x, nb, binned)` (R code for the binned case).
fn bw_pair_cnts(x: &[f64], nb: usize, binned: bool) -> Result<(f64, Vec<f64>), String> {
    if !binned {
        return bw_den(nb, x);
    }
    // r <- range(x); d <- diff(r) * 1.01/nb
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in x {
        if v.is_nan() {
            return Err("NA in x".into());
        }
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let d = (hi - lo) * 1.01 / nb as f64;
    // xx <- trunc(abs(x)/d) * sign(x); xx <- xx - min(xx) + 1
    let sign = |v: f64| if v > 0.0 { 1.0 } else if v < 0.0 { -1.0 } else { 0.0 };
    let xx: Vec<f64> = x.iter().map(|&v| (v.abs() / d).trunc() * sign(v)).collect();
    let m = xx.iter().copied().fold(f64::INFINITY, f64::min);
    // tabulate(xx, nb): as.integer() (truncation), bins 1..nb, others ignored
    let mut tab = vec![0i32; nb];
    for &v in &xx {
        let b = v - m + 1.0;
        if !(b.abs() < 2147483648.0) {
            continue; // NA_INTEGER
        }
        let b = b as i32;
        if b > 0 && b as usize <= nb {
            tab[b as usize - 1] += 1;
        }
    }
    Ok((d, bw_den_binned(&tab)))
}

/// `bw_phi4(n, d, cnt, h)`.
fn bw_phi4(n: usize, d: f64, cnt: &[f64], h: f64) -> f64 {
    let mut sum = 0.0f64;
    for (i, &c) in cnt.iter().enumerate() {
        let mut delta = i as f64 * d / h;
        delta *= delta;
        if delta >= DELMAX {
            break;
        }
        let term = exp(-delta / 2.0) * (delta * delta - 6.0 * delta + 3.0);
        sum += term * c;
    }
    let nf = n as f64;
    sum = 2.0 * sum + nf * 3.0; // the diagonal
    sum / (nf * (n - 1) as f64 * pow(h, 5.0)) * M_1_SQRT_2PI
}

/// `bw_phi6(n, d, cnt, h)`.
fn bw_phi6(n: usize, d: f64, cnt: &[f64], h: f64) -> f64 {
    let mut sum = 0.0f64;
    for (i, &c) in cnt.iter().enumerate() {
        let mut delta = i as f64 * d / h;
        delta *= delta;
        if delta >= DELMAX {
            break;
        }
        let term = exp(-delta / 2.0) * (delta * delta * delta - 15.0 * delta * delta + 45.0 * delta - 15.0);
        sum += term * c;
    }
    let nf = n as f64;
    sum = 2.0 * sum - 15.0 * nf; // the diagonal
    sum / (nf * (n - 1) as f64 * pow(h, 7.0)) * M_1_SQRT_2PI
}

/// `bw.SJ(x)` (method "ste", `nb = 1000`, the default search interval and tolerance).
pub fn bw_sj(x: &[f64]) -> Result<f64, String> {
    let n = x.len();
    if n < 2 {
        return Err("need at least 2 data points".into());
    }
    let nb = 1000usize;
    let (d, cnt) = bw_pair_cnts(x, nb, n > nb / 2)?;
    let sdh = |h: f64| bw_phi4(n, d, &cnt, h);
    let tdh = |h: f64| bw_phi6(n, d, &cnt, h);
    // scale <- min(sd(x), IQR(x)/1.349)
    let (s1, s2) = (sd(x, false), iqr(x) / 1.349);
    let scale = if s1.is_nan() || s2.is_nan() { f64::NAN } else if s2 < s1 { s2 } else { s1 };
    let nf = n as f64;
    let a = 1.24 * scale * r_pow(nf, -1.0 / 7.0);
    let b = 1.23 * scale * r_pow(nf, -1.0 / 9.0);
    let c1 = 1.0 / (2.0 * std::f64::consts::PI.sqrt() * nf);
    let td = -tdh(b);
    if !td.is_finite() || td <= 0.0 {
        return Err("sample is too sparse to find TD".into());
    }
    let hmax = 1.144 * scale * r_pow(nf, -1.0 / 5.0);
    let mut lower = 0.1 * hmax;
    let mut upper = hmax;
    let alph2 = 1.357 * r_pow(sdh(a) / td, 1.0 / 7.0);
    if !alph2.is_finite() {
        return Err("sample is too sparse to find alph2".into());
    }
    let fsd = |h: f64| r_pow(c1 / sdh(alph2 * r_pow(h, 5.0 / 7.0)), 1.0 / 5.0) - h;
    let mut itry = 1;
    while fsd(lower) * fsd(upper) > 0.0 {
        if itry > 99 {
            return Err("no solution in the specified range of bandwidths".into());
        }
        if itry % 2 == 1 {
            upper *= 1.2;
        } else {
            lower /= 1.2;
        }
        itry += 1;
    }
    // tol = 0.1 * lower, a promise forced only now, after the search interval was widened
    let tol = 0.1 * lower;
    Ok(uniroot(&mut |h| fsd(h), lower, upper, tol, 1000)?.root)
}
