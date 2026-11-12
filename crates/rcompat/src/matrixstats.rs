//! matrixStats 1.5.0 as glmGamPoi and sctransform reach it (directly or through MatrixGenerics
//! and DelayedMatrixStats, which call matrixStats for ordinary matrices), and R's `R_qsort_I`
//! (`src/main/qsort.c`), which `weightedMedian` sorts with.
//!
//! - [`row_means2`], [`col_means2`]: `rowMeans2`/`colMeans2` (`refine = TRUE`): the sum in
//!   index order, the mean, then the mean of the residuals added on. matrixStats declares these
//!   accumulators `LDOUBLE`, but that is `long double` only under `HAVE_LONG_DOUBLE`, which R
//!   defines in its private `config.h`, not in the headers packages are built with: the installed
//!   package computes in double (base R's own `rowMeans` uses long double);
//! - [`row_vars`]: `rowVars` (`refine = TRUE`): the same three passes in double;
//! - [`weighted_median`]: `weightedMedian(x, w, interpolate = FALSE)` (ties "weighted");
//! - [`weighted_median_interpolated`]: `weightedMedian(x, w)`, whose `interpolate` defaults to
//!   TRUE (`is.null(ties)`).
//!
//! Matrices are column-major (`nrow x ncol`), as R stores them.

fn mean_refined<I: Iterator<Item = f64> + Clone>(values: I) -> f64 {
    // matrixStats' LDOUBLE is `long double` only when HAVE_LONG_DOUBLE is defined, which is in
    // R's private config.h, not in the headers packages see: the package is built with double.
    let mut sum = 0.0f64;
    let mut count: u64 = 0;
    for v in values.clone() {
        sum += v;
        count += 1;
    }
    if sum > f64::MAX {
        return f64::INFINITY;
    }
    if sum < -f64::MAX {
        return f64::NEG_INFINITY;
    }
    let n = count as f64;
    let avg = sum / n;
    let mut s = 0.0f64;
    for v in values {
        s += v - avg;
    }
    avg + s / n
}

/// `rowMeans2(x, cols = cols)`: one mean per row over the columns `cols` (all when `None`), in
/// that order.
pub fn row_means2(x: &[f64], nrow: usize, ncol: usize, cols: Option<&[usize]>) -> Vec<f64> {
    assert_eq!(x.len(), nrow * ncol);
    let all: Vec<usize>;
    let cols = match cols {
        Some(c) => c,
        None => {
            all = (0..ncol).collect();
            &all
        }
    };
    (0..nrow).map(|i| mean_refined(cols.iter().map(move |&j| x[i + j * nrow]))).collect()
}

/// `rowMeans2(x, na.rm = TRUE)` of one row's values (in column order), as Seurat's
/// `FastRowScale` centres: NaN values are skipped in the sum, the count and the refinement.
pub fn mean2_na_rm(values: &[f64]) -> f64 {
    let mut sum = 0.0f64;
    let mut count: u64 = 0;
    for &v in values {
        if !v.is_nan() {
            sum += v;
            count += 1;
        }
    }
    if sum > f64::MAX {
        return f64::INFINITY;
    }
    if sum < -f64::MAX {
        return f64::NEG_INFINITY;
    }
    let n = count as f64;
    let avg = sum / n;
    let mut s = 0.0f64;
    for &v in values {
        if !v.is_nan() {
            s += v - avg;
        }
    }
    avg + s / n
}

/// `colMeans2(x)`.
pub fn col_means2(x: &[f64], nrow: usize, ncol: usize) -> Vec<f64> {
    assert_eq!(x.len(), nrow * ncol);
    (0..ncol).map(|j| mean_refined(x[j * nrow..(j + 1) * nrow].iter().copied())).collect()
}

/// `rowVars` of one row's values (in column order), `refine = TRUE`; NA for fewer than 2 values.
pub fn var_refined(values: &[f64]) -> f64 {
    let kk = values.len();
    if kk <= 1 {
        return crate::summary::NA_REAL;
    }
    let mut sum = 0.0;
    for &v in values {
        sum += v;
    }
    let mut mu = sum / kk as f64;
    sum = 0.0;
    for &v in values {
        sum += v - mu;
    }
    mu += sum / kk as f64;
    let mut sigma2 = 0.0;
    for &v in values {
        let d = v - mu;
        sigma2 += d * d;
    }
    sigma2 / (kk - 1) as f64
}

/// `rowVars(x)`.
pub fn row_vars(x: &[f64], nrow: usize, ncol: usize) -> Vec<f64> {
    assert_eq!(x.len(), nrow * ncol);
    let mut row = vec![0.0; ncol];
    (0..nrow)
        .map(|i| {
            for j in 0..ncol {
                row[j] = x[i + j * nrow];
            }
            var_refined(&row)
        })
        .collect()
}

/// R's `R_qsort_I(v, I, 1, n)`: sorts `v` increasingly (Singleton's quicksort, not stable) and
/// permutes `idx` alongside.
pub fn r_qsort_i(v: &mut [f64], idx: &mut [i32]) {
    let n = v.len();
    assert_eq!(n, idx.len());
    if n < 2 {
        return;
    }
    // 1-based views
    let mut vv = Vec::with_capacity(n + 1);
    vv.push(0.0);
    vv.extend_from_slice(v);
    let mut ii_ = Vec::with_capacity(n + 1);
    ii_.push(0);
    ii_.extend_from_slice(idx);
    let (v1, i1) = (&mut vv, &mut ii_);
    let mut il = [0i64; 40];
    let mut iu = [0i64; 40];
    let mut r = 0.375f64;
    let (mut i, mut j) = (1i64, n as i64);
    let ii = i;
    let mut m = 1usize;
    #[derive(PartialEq)]
    enum L {
        L10,
        L20,
        L80,
        After,
    }
    let mut lab = L::L10;
    loop {
        match lab {
            L::L10 => {
                if i < j {
                    if r < 0.5898437 {
                        r += 0.0390625;
                    } else {
                        r -= 0.21875;
                    }
                    lab = L::L20;
                } else {
                    lab = L::L80;
                }
                continue;
            }
            L::L20 => {
                let mut k = i;
                let ij = i + ((j - i) as f64 * r) as i64;
                let mut it = i1[ij as usize];
                let mut vt = v1[ij as usize];
                if v1[i as usize] > vt {
                    i1[ij as usize] = i1[i as usize];
                    i1[i as usize] = it;
                    it = i1[ij as usize];
                    v1[ij as usize] = v1[i as usize];
                    v1[i as usize] = vt;
                    vt = v1[ij as usize];
                }
                let mut l = j;
                if v1[j as usize] < vt {
                    i1[ij as usize] = i1[j as usize];
                    i1[j as usize] = it;
                    it = i1[ij as usize];
                    v1[ij as usize] = v1[j as usize];
                    v1[j as usize] = vt;
                    vt = v1[ij as usize];
                    if v1[i as usize] > vt {
                        i1[ij as usize] = i1[i as usize];
                        i1[i as usize] = it;
                        v1[ij as usize] = v1[i as usize];
                        v1[i as usize] = vt;
                        vt = v1[ij as usize];
                    }
                }
                loop {
                    loop {
                        l -= 1;
                        if !(v1[l as usize] > vt) {
                            break;
                        }
                    }
                    let tt = i1[l as usize];
                    let vtt = v1[l as usize];
                    loop {
                        k += 1;
                        if !(v1[k as usize] < vt) {
                            break;
                        }
                    }
                    if k > l {
                        break;
                    }
                    i1[l as usize] = i1[k as usize];
                    i1[k as usize] = tt;
                    v1[l as usize] = v1[k as usize];
                    v1[k as usize] = vtt;
                }
                m += 1;
                if l - i <= j - k {
                    il[m] = k;
                    iu[m] = j;
                    j = l;
                } else {
                    il[m] = i;
                    iu[m] = l;
                    i = k;
                }
                lab = L::After;
                continue;
            }
            L::L80 => {
                if m == 1 {
                    break;
                }
                i = il[m];
                j = iu[m];
                m -= 1;
                lab = L::After;
                continue;
            }
            L::After => {
                if j - i > 10 {
                    lab = L::L20;
                    continue;
                }
                if i == ii {
                    lab = L::L10;
                    continue;
                }
                // insertion sort of v[i..j]
                i -= 1;
                loop {
                    // L100
                    let mut it;
                    let mut vt;
                    loop {
                        i += 1;
                        if i == j {
                            break;
                        }
                        it = i1[(i + 1) as usize];
                        vt = v1[(i + 1) as usize];
                        if !(v1[i as usize] <= vt) {
                            let mut k = i;
                            loop {
                                i1[(k + 1) as usize] = i1[k as usize];
                                v1[(k + 1) as usize] = v1[k as usize];
                                k -= 1;
                                if !(vt < v1[k as usize]) {
                                    break;
                                }
                            }
                            i1[(k + 1) as usize] = it;
                            v1[(k + 1) as usize] = vt;
                        }
                    }
                    break;
                }
                lab = L::L80;
                continue;
            }
        }
    }
    v.copy_from_slice(&v1[1..]);
    idx.copy_from_slice(&i1[1..]);
}

/// `weightedMedian(x, w, interpolate = FALSE)` (ties weighted) for finite `x` and finite positive
/// weights, as glmGamPoi's `loc_median_fit` calls it.
pub fn weighted_median(x: &[f64], w: &[f64]) -> f64 {
    assert_eq!(x.len(), w.len());
    // drop non-positive weights (and NaN data, which the callers never pass)
    let mut xt = Vec::with_capacity(x.len());
    let mut wt = Vec::with_capacity(x.len());
    let mut wtotal = 0.0;
    for (&xi, &wi) in x.iter().zip(w) {
        if wi > 0.0 {
            assert!(wi.is_finite() && !xi.is_nan(), "weighted_median: infinite weight or NaN value");
            xt.push(xi);
            wt.push(wi);
            wtotal += wi;
        }
    }
    let nx = xt.len();
    if nx == 0 {
        return crate::summary::NA_REAL;
    }
    if nx == 1 {
        return xt[0];
    }
    let mut idx: Vec<i32> = (0..nx as i32).collect();
    r_qsort_i(&mut xt, &mut idx);
    let mut wcum = vec![0.0; nx];
    let mut tmp = 0.0;
    let mut half = nx + 1;
    for ii in 0..nx {
        tmp += wt[idx[ii] as usize] / wtotal;
        wcum[ii] = tmp;
        if tmp > 0.5 {
            half = ii;
            break;
        }
    }
    if half == 0 || half == nx {
        return xt[half];
    }
    if half > nx {
        // the cumulative weight never passed 1/2 (rounding): C reads x[nx + 1], out of bounds
        panic!("weighted_median: cumulative weight never exceeded 1/2");
    }
    let wlow = wcum[half - 1];
    let whigh = 1.0 - wlow;
    if whigh > 0.5 {
        xt[half]
    } else {
        wlow * xt[half - 1] + whigh * xt[half]
    }
}

/// `weightedMedian(x, w)` with its default `interpolate = TRUE`, for finite `x` and finite
/// positive weights: glmGamPoi's `loc_median_fit` calls it when its window is wider than the data.
pub fn weighted_median_interpolated(x: &[f64], w: &[f64]) -> f64 {
    assert_eq!(x.len(), w.len());
    let mut xt = Vec::with_capacity(x.len());
    let mut wt = Vec::with_capacity(x.len());
    let mut wtotal = 0.0;
    for (&xi, &wi) in x.iter().zip(w) {
        if wi > 0.0 {
            assert!(wi.is_finite() && !xi.is_nan(), "weighted_median_interpolated: infinite weight or NaN value");
            xt.push(xi);
            wt.push(wi);
            wtotal += wi;
        }
    }
    let nx = xt.len();
    if nx == 0 {
        return crate::summary::NA_REAL;
    }
    if nx == 1 {
        return xt[0];
    }
    let mut idx: Vec<i32> = (0..nx as i32).collect();
    r_qsort_i(&mut xt, &mut idx);
    let mut wcum = vec![0.0; nx];
    let mut tmp2 = 0.0;
    let mut half = nx + 1;
    for ii in 0..nx {
        let tmp = wt[idx[ii] as usize] / wtotal;
        tmp2 += tmp;
        wcum[ii] = tmp2 - (tmp / 2.0);
        if wcum[ii] >= 0.5 {
            half = ii;
            break;
        }
    }
    if half == 0 {
        return xt[0];
    }
    if half > nx {
        // the adjusted cumulative weight never reached 1/2 (rounding): C reads x[nx], out of bounds
        panic!("weighted_median_interpolated: cumulative weight never reached 1/2");
    }
    // the width and height of the rectangle, then of the triangle up to half the total weight
    let dx = xt[half] - xt[half - 1];
    let dy_big = wcum[half] - wcum[half - 1];
    let dy = 0.5 - wcum[half];
    let dx = (dy / dy_big) * dx;
    dx + xt[half]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qsort_sorts() {
        let mut v = vec![5.0, 1.0, 4.0, 1.0, 3.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 8.0, 9.0, 7.0, 9.0, 3.0];
        let mut idx: Vec<i32> = (0..v.len() as i32).collect();
        let orig = v.clone();
        r_qsort_i(&mut v, &mut idx);
        for k in 1..v.len() {
            assert!(v[k - 1] <= v[k]);
        }
        for k in 0..v.len() {
            assert_eq!(orig[idx[k] as usize], v[k]);
        }
    }

    #[test]
    fn medians() {
        assert_eq!(weighted_median(&[1.0, 2.0, 3.0], &[1.0, 1.0, 1.0]), 2.0);
        assert_eq!(row_means2(&[1.0, 2.0, 3.0, 4.0], 2, 2, None), vec![2.0, 3.0]);
        assert_eq!(var_refined(&[1.0, 2.0, 3.0]), 1.0);
    }
}
