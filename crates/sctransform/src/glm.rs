//! glmGamPoi 1.22.0's `glm_gp(Y, design = ~ 1, offset = log(umi), size_factors = FALSE)` as
//! sctransform's `fit_glmGamPoi_offset` calls it (stage S2b): the intercept-only model takes the
//! group-wise path (one group), so the fit is
//!
//! 1. `estimate_dispersions_roughly`: moment estimates `(rowVars(Y) - xim * rowMeans2(Y)) /
//!    rowMeans2(Y)^2` with `xim = 1 / mean(colMeans2(exp(offset)))`, negatives and NA set to 0;
//! 2. `estimate_betas_roughly_group_wise`: `log(rowMeans2(Y / exp(offset)))`;
//! 3. `fitBeta_one_group` (C++): Newton-Raphson for each gene's intercept at those dispersions,
//!    and where a step makes it NaN (or 100 iterations pass), `optimize()` over the negative
//!    binomial log-likelihood (R's `dnbinom` and Brent's method);
//! 4. `overdispersion_mle` per gene: the Cox-Reid adjusted profile likelihood in `log(theta)`,
//!    maximised by `nlminb` with analytic score and Hessian (`conventional_*_fast`, C++), over a
//!    count table built in a libstdc++ `unordered_map` (whose iteration order sets the summation
//!    order, and whose integer keys truncate half-integer counts);
//! 5. `overdispersion_shrinkage`: of its results only the dispersion trend is used, a weighted
//!    running median (`loc_median_fit`) of the MLEs against the genes' mean `Mu`;
//! 6. `fitBeta_one_group` again at the trend.
//!
//! The function returns the second fit's intercepts and the MLE dispersions (not the trend).
//! Every reduction is in the order of the C++ or R code; matrixStats computes in double;
//! `mean()` and `var()` are R's long-double ones; the Cox-Reid sums are OpenBLAS `dgemv` calls,
//! as Armadillo makes them.

use rcompat::blas;
use rcompat::gamma::{digamma, lgamma_glibc, lgammafn, trigamma};
use rcompat::libm::{exp, log};
use rcompat::matrixstats::{var_refined, weighted_median, weighted_median_interpolated};
use rcompat::nbinom::dnbinom_mu_log;
use rcompat::optimize::{optimize_maximum, OPTIMIZE_TOL};
use rcompat::port::{nlminb, Problem};
use rcompat::{dnorm, mean, order, seq_length_out, var};

/// `cr_correction_factor` of glmGamPoi's overdispersion.cpp.
const CR_CORRECTION_FACTOR: f64 = 0.99;

// --- the count table ----------------------------------------------------------------------------

/// libstdc++'s `std::__detail::__prime_list` (64-bit), read from the R env's `libstdc++.so.6`:
/// the bucket counts `_Prime_rehash_policy::_M_next_bkt` chooses from.
const PRIME_LIST: [u64; 304] = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 103, 109, 113, 127, 137, 139, 149,
    157, 167, 179, 193, 199, 211, 227, 241,
    257, 277, 293, 313, 337, 359, 383, 409,
    439, 467, 503, 541, 577, 619, 661, 709,
    761, 823, 887, 953, 1031, 1109, 1193, 1289,
    1381, 1493, 1613, 1741, 1879, 2029, 2179, 2357,
    2549, 2753, 2971, 3209, 3469, 3739, 4027, 4349,
    4703, 5087, 5503, 5953, 6427, 6949, 7517, 8123,
    8783, 9497, 10273, 11113, 12011, 12983, 14033, 15173,
    16411, 17749, 19183, 20753, 22447, 24281, 26267, 28411,
    30727, 33223, 35933, 38873, 42043, 45481, 49201, 53201,
    57557, 62233, 67307, 72817, 78779, 85229, 92203, 99733,
    107897, 116731, 126271, 136607, 147793, 159871, 172933, 187091,
    202409, 218971, 236897, 256279, 277261, 299951, 324503, 351061,
    379787, 410857, 444487, 480881, 520241, 562841, 608903, 658753,
    712697, 771049, 834181, 902483, 976369, 1056323, 1142821, 1236397,
    1337629, 1447153, 1565659, 1693859, 1832561, 1982627, 2144977, 2320627,
    2510653, 2716249, 2938679, 3179303, 3439651, 3721303, 4026031, 4355707,
    4712381, 5098259, 5515729, 5967347, 6456007, 6984629, 7556579, 8175383,
    8844859, 9569143, 10352717, 11200489, 12117689, 13109983, 14183539, 15345007,
    16601593, 17961079, 19431899, 21023161, 22744717, 24607243, 26622317, 28802401,
    31160981, 33712729, 36473443, 39460231, 42691603, 46187573, 49969847, 54061849,
    58488943, 63278561, 68460391, 74066549, 80131819, 86693767, 93793069, 101473717,
    109783337, 118773397, 128499677, 139022417, 150406843, 162723577, 176048909, 190465427,
    206062531, 222936881, 241193053, 260944219, 282312799, 305431229, 330442829, 357502601,
    386778277, 418451333, 452718089, 489790921, 529899637, 573292817, 620239453, 671030513,
    725980837, 785430967, 849749479, 919334987, 994618837, 1076067617, 1164186217, 1259520799,
    1362662261, 1474249943, 1594975441, 1725587117, 1866894511, 2019773507, 2185171673, 2364114217,
    2557710269, 2767159799, 2993761039, 3238918481, 3504151727, 3791104843, 4101556399, 4294967291,
    6442450933, 8589934583, 12884901857, 17179869143, 25769803693, 34359738337, 51539607367, 68719476731,
    103079215087, 137438953447, 206158430123, 274877906899, 412316860387, 549755813881, 824633720731, 1099511627689,
    1649267441579, 2199023255531, 3298534883309, 4398046511093, 6597069766607, 8796093022151, 13194139533241, 17592186044399,
    26388279066581, 35184372088777, 52776558133177, 70368744177643, 105553116266399, 140737488355213, 211106232532861, 281474976710597,
    562949953421231, 1125899906842597, 2251799813685119, 4503599627370449, 9007199254740881, 18014398509481951, 36028797018963913, 72057594037927931,
    144115188075855859, 288230376151711717, 576460752303423433, 1152921504606846883, 2305843009213693951, 4611686018427387847, 9223372036854775783, 18446744073709551557,
];

/// `_Prime_rehash_policy::_M_next_bkt(n)` (max load factor 1).
fn next_bkt(n: u64) -> u64 {
    const FAST_BKT: [u64; 14] = [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13];
    if n < 14 {
        return if n == 0 { 1 } else { FAST_BKT[n as usize] };
    }
    // lower_bound(__prime_list + 6, last prime, n), the last prime as the past-the-end value
    let last = PRIME_LIST.len() - 2;
    let slice = &PRIME_LIST[6..last];
    let k = slice.partition_point(|&p| p < n);
    PRIME_LIST[6 + k]
}

/// glmGamPoi's `make_table_if_small(x, stop_if_larger)`: the frequency of each count `(long) x`,
/// in the iteration order of a libstdc++ `std::unordered_map<long, size_t>` reserved for
/// `stop_if_larger` elements; empty if more than `stop_if_larger` distinct keys occur.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CountTable {
    pub keys: Vec<f64>,
    pub freqs: Vec<f64>,
}

impl CountTable {
    pub fn new(x: &[f64], stop_if_larger: i32) -> CountTable {
        let stop = stop_if_larger.max(0) as usize;
        let nb = next_bkt(stop as u64) as usize;
        // singly linked list of nodes; buckets hold the node *before* their first node
        // (`None` = empty bucket, `Some(usize::MAX)` = the list head `_M_before_begin`)
        const HEAD: usize = usize::MAX;
        let mut buckets: Vec<Option<usize>> = vec![None; nb];
        let mut head_next: Option<usize> = None;
        let mut next: Vec<Option<usize>> = Vec::new();
        let mut keys: Vec<i64> = Vec::new();
        let mut counts: Vec<u64> = Vec::new();
        let bucket_of = |k: i64| (k as u64 % nb as u64) as usize;
        for &v in x {
            let key = v as i64; // (long) v: truncation toward zero
            let b = bucket_of(key);
            // find the key in its bucket
            let mut found = None;
            if let Some(before) = buckets[b] {
                let mut cur = if before == HEAD { head_next } else { next[before] };
                while let Some(c) = cur {
                    if bucket_of(keys[c]) != b {
                        break;
                    }
                    if keys[c] == key {
                        found = Some(c);
                        break;
                    }
                    cur = next[c];
                }
            }
            match found {
                Some(c) => counts[c] += 1,
                None => {
                    let node = keys.len();
                    keys.push(key);
                    counts.push(1);
                    next.push(None);
                    // _M_insert_bucket_begin
                    match buckets[b] {
                        Some(before) => {
                            let after = if before == HEAD { head_next } else { next[before] };
                            next[node] = after;
                            if before == HEAD {
                                head_next = Some(node);
                            } else {
                                next[before] = Some(node);
                            }
                        }
                        None => {
                            next[node] = head_next;
                            head_next = Some(node);
                            if let Some(nx) = next[node] {
                                buckets[bucket_of(keys[nx])] = Some(node);
                            }
                            buckets[b] = Some(HEAD);
                        }
                    }
                    if keys.len() > stop {
                        return CountTable::default();
                    }
                }
            }
        }
        let mut out = CountTable { keys: Vec::with_capacity(keys.len()), freqs: Vec::with_capacity(keys.len()) };
        let mut cur = head_next;
        while let Some(c) = cur {
            // the transform's lambda takes std::pair<int, size_t>: the key through int
            out.keys.push(keys[c] as i32 as f64);
            out.freqs.push(counts[c] as f64);
            cur = next[c];
        }
        out
    }

    fn present(&self) -> bool {
        !self.keys.is_empty()
    }
}

// --- the profile likelihood and its derivatives ---------------------------------------------------

/// Armadillo's `model_matrix.t() * (model_matrix.each_col() % v)` for the intercept-only design
/// (a column of ones): `gemv` with the transpose, i.e. OpenBLAS `dgemv("T", n, 1, 1, v, n, 1, 1,
/// 0, y, 1)`, for n > 1 (a 1 x 1 matrix is multiplied without BLAS).
fn cr_sum(v: &[f64], ones: &[f64]) -> f64 {
    let n = v.len();
    if n == 1 {
        // gemv_emul_tinysq: y = A * x
        return v[0] * ones[0];
    }
    let mut y = [0.0f64];
    blas::dgemv(b'T', n, 1, 1.0, v, n, ones, 1, 0.0, &mut y, 1);
    y[0]
}

/// `conventional_loglikelihood_fast(y, mu, log_theta, model_matrix, do_cr_adj, keys, freqs)`.
pub fn loglik(y: &[f64], mu: &[f64], log_theta: f64, cr: bool, tab: &CountTable, ones: &[f64]) -> f64 {
    let theta = exp(log_theta);
    let mut cr_term = 0.0;
    if cr {
        let w: Vec<f64> = mu.iter().map(|&m| 1.0 / (1.0 / m + theta)).collect();
        let b = cr_sum(&w, ones);
        // lu(L, U, P, b) of the 1 x 1 matrix: L = 1, U = b
        let mut ld = 0.0 + log(1.0); // sum(log(diagvec(L)))
        ld += if b < 1e-50 { log(1e-50) } else { log(b) };
        cr_term = -0.5 * ld * CR_CORRECTION_FACTOR;
    }
    let theta_neg1 = 1.0 / theta; // R_pow_di(theta, -1)
    let n = y.len() as f64;
    let mut lgamma_term = 0.0;
    if tab.present() {
        for (&u, &f) in tab.keys.iter().zip(&tab.freqs) {
            lgamma_term += f * lgamma_glibc(u + theta_neg1);
        }
    } else {
        // Rcpp sugar: sum(lgamma(y + theta_neg1)), R's lgammafn
        for &v in y {
            lgamma_term += lgammafn(v + theta_neg1);
        }
    }
    lgamma_term -= n * lgamma_glibc(theta_neg1);
    let mut ll_part = 0.0;
    for (&v, &m) in y.iter().zip(mu) {
        ll_part += (-v - theta_neg1) * log(m + theta_neg1);
    }
    ll_part -= n * theta_neg1 * log(theta);
    lgamma_term + ll_part + cr_term
}

/// `conventional_score_function_fast(...)`.
pub fn score(y: &[f64], mu: &[f64], log_theta: f64, cr: bool, tab: &CountTable, ones: &[f64]) -> f64 {
    let theta = exp(log_theta);
    let theta_neg1 = 1.0 / theta;
    let mut cr_term = 0.0;
    if cr {
        let w: Vec<f64> = mu.iter().map(|&m| 1.0 / (1.0 / m + theta)).collect();
        let dw: Vec<f64> = w.iter().map(|&wi| -1.0 * wi * wi).collect();
        let b = cr_sum(&w, ones);
        let db = cr_sum(&dw, ones);
        // inv_sympd(b + eye * 1e-6) of a 1 x 1 matrix: 1 / a
        let b_inv = 1.0 / (b + 1.0 * 1e-6);
        // trace(b_inv * db): acc1 = 0 + b_inv * db; acc = acc1 + 0
        let tr = (0.0 + b_inv * db) + 0.0;
        cr_term = -0.5 * tr * CR_CORRECTION_FACTOR;
    }
    let n = y.len() as f64;
    let mut digamma_term = 0.0;
    let (mut max_y, mut sum_y, mut sum_prod_y) = (0.0f64, 0.0f64, 0.0f64);
    if tab.present() {
        for (&u, &f) in tab.keys.iter().zip(&tab.freqs) {
            digamma_term += f * digamma(u + theta_neg1);
            sum_y += f * u;
            sum_prod_y += f * (u - 1.0) * u;
            max_y = if max_y < u { u } else { max_y }; // std::max(max_y, u)
        }
    } else {
        for &v in y {
            digamma_term += digamma(v + theta_neg1);
            sum_y += v;
            sum_prod_y += (v - 1.0) * v;
            max_y = if max_y < v { v } else { max_y };
        }
    }
    let corr = if theta_neg1 > 1e5 { sum_prod_y / (2.0 * theta_neg1) } else { 0.0 };
    if max_y * 1e6 < theta_neg1 {
        digamma_term = sum_y - corr;
    } else {
        digamma_term -= n * digamma(theta_neg1);
        digamma_term *= theta_neg1;
        let cap = sum_y - corr;
        // std::min(digamma_term, cap): (cap < digamma_term) ? cap : digamma_term
        digamma_term = if cap < digamma_term { cap } else { digamma_term };
    }
    let mut ll_part = 0.0;
    for (&v, &m) in y.iter().zip(mu) {
        let mu_theta = m * theta;
        if mu_theta < 1e-10 {
            ll_part += mu_theta * mu_theta * (1.0 / (1.0 + mu_theta) - 0.5);
        } else if mu_theta < 1e-4 {
            let inv = 1.0 / (1.0 + mu_theta);
            let upper = mu_theta * mu_theta * inv;
            let lower = mu_theta * mu_theta * (inv - 0.5);
            let suggest = log(1.0 + mu_theta) - m / (m + theta_neg1);
            // std::max(std::min(suggest, upper), lower)
            let mn = if upper < suggest { upper } else { suggest };
            ll_part += if mn < lower { lower } else { mn };
        } else {
            ll_part += log(1.0 + mu_theta) - m / (m + theta_neg1);
        }
        ll_part += v / (m + theta_neg1);
    }
    ll_part *= theta_neg1;
    ll_part - digamma_term + cr_term * theta
}

/// `conventional_deriv_score_function_fast(...)`.
pub fn deriv_score(y: &[f64], mu: &[f64], log_theta: f64, cr: bool, tab: &CountTable, ones: &[f64]) -> f64 {
    let theta = exp(log_theta);
    let (mut cr_term, mut cr_term2) = (0.0, 0.0);
    if cr {
        let w: Vec<f64> = mu.iter().map(|&m| 1.0 / (1.0 / m + theta)).collect();
        let dw: Vec<f64> = w.iter().map(|&wi| -1.0 * wi * wi).collect();
        let d2w: Vec<f64> = dw.iter().zip(&w).map(|(&dwi, &wi)| -2.0 * dwi * wi).collect();
        let b = cr_sum(&w, ones);
        let db = cr_sum(&dw, ones);
        let d2b = cr_sum(&d2w, ones);
        let b_inv = 1.0 / (b + 1.0 * 1e-6);
        // d_i_db = b_inv * db (1 x 1: gemv_emul_tinysq)
        let d_i_db = b_inv * db;
        let ddetb = (0.0 + d_i_db) + 0.0; // trace(d_i_db)
        let tr_sq = (0.0 + d_i_db * d_i_db) + 0.0; // trace(d_i_db * d_i_db)
        let tr_b_d2b = (0.0 + b_inv * d2b) + 0.0; // trace(b_inv * d2b)
        let d2detb = ddetb * ddetb - tr_sq + tr_b_d2b;
        cr_term = (0.5 * (ddetb * ddetb) - 0.5 * d2detb) * CR_CORRECTION_FACTOR;
        cr_term2 = -0.5 * ddetb * CR_CORRECTION_FACTOR;
    }
    let theta_neg1 = 1.0 / theta;
    let theta_neg2 = 1.0 / (theta * theta); // R_pow_di(theta, -2)
    let n = y.len() as f64;
    let (mut digamma_term, mut trigamma_term) = (0.0f64, 0.0f64);
    if tab.present() {
        for (&u, &f) in tab.keys.iter().zip(&tab.freqs) {
            digamma_term += f * digamma(u + theta_neg1);
            trigamma_term += f * trigamma(u + theta_neg1);
        }
        trigamma_term *= theta_neg2;
        digamma_term -= n * digamma(theta_neg1);
        trigamma_term -= theta_neg2 * n * trigamma(theta_neg1);
    } else {
        let mut s = 0.0;
        for &v in y {
            s += digamma(v + theta_neg1);
        }
        digamma_term = s;
        digamma_term -= n * digamma(theta_neg1);
        let mut s2 = 0.0;
        for &v in y {
            s2 += trigamma(v + theta_neg1);
        }
        trigamma_term = theta_neg2 * s2;
        trigamma_term -= theta_neg2 * n * trigamma(theta_neg1);
    }
    let (mut ll_part_1, mut ll_part_2) = (0.0, 0.0);
    for (&v, &m) in y.iter().zip(mu) {
        ll_part_1 += log(1.0 + m * theta) + (v - m) / (m + theta_neg1);
        ll_part_2 += (m * m * theta + v) / (1.0 + m * theta) / (1.0 + m * theta);
    }
    let ll_part = -2.0 * theta_neg1 * (ll_part_1 - digamma_term) + (ll_part_2 + trigamma_term);
    ll_part + cr_term * (theta * theta) + (ll_part_1 - digamma_term) * theta_neg1 + cr_term2 * theta
}

// --- the overdispersion MLE ---------------------------------------------------------------------

/// The result of [`overdispersion_mle`].
#[derive(Debug, Clone, PartialEq)]
pub struct Overdispersion {
    pub estimate: f64,
    pub iterations: i32,
    /// nlminb's convergence code (`None` when no optimisation ran)
    pub convergence: Option<i32>,
}

struct Profile<'a> {
    y: &'a [f64],
    mu: &'a [f64],
    cr: bool,
    tab: &'a CountTable,
    ones: &'a [f64],
}

impl Problem for Profile<'_> {
    fn objective(&mut self, x: &[f64]) -> f64 {
        -loglik(self.y, self.mu, x[0], self.cr, self.tab, self.ones)
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) {
        g[0] = -score(self.y, self.mu, x[0], self.cr, self.tab, self.ones);
    }
    fn hessian(&mut self, x: &[f64], h: &mut [f64]) {
        h[0] = -deriv_score(self.y, self.mu, x[0], self.cr, self.tab, self.ones);
    }
}

/// `overdispersion_mle_impl(y, mu, model_matrix = ~1, do_cox_reid_adjustment = TRUE,
/// n_subsamples = length(y), max_iter = 200)` for one gene (glmGamPoi's
/// `conventional_overdispersion_mle`).
pub fn overdispersion_mle(y: &[f64], mu_in: &[f64], ones: &[f64]) -> Result<Overdispersion, String> {
    let n = y.len();
    let tab = CountTable::new(y, (n as f64 / 2.0) as i32);
    if y.iter().all(|&v| v == 0.0) {
        return Ok(Overdispersion { estimate: 0.0, iterations: 0, convergence: None });
    }
    let mu: Vec<f64> = mu_in.iter().map(|&m| if m == 0.0 { 1e-6 } else { m }).collect();
    let far_left = score(y, &mu, log(1e-8), true, &tab, ones);
    if far_left < 0.0 {
        return Ok(Overdispersion { estimate: 0.0, iterations: 0, convergence: None });
    }
    let m = mean(y, false);
    let mut start = (var(y, false) - m) / (m * m);
    if start.is_nan() || start <= 0.0 {
        start = 0.5;
    }
    let (lower, upper) = (log(1e-16), log(1e16));
    let x0 = [log(start)];
    let mut prob = Profile { y, mu: &mu, cr: true, tab: &tab, ones };
    let mut res = nlminb(&x0, &mut prob, &[lower], &[upper], Some(200), true).map_err(|e| e.to_string())?;
    if res.convergence != 0 {
        // "Do the same thing again with numerical hessian": nlminb without the Hessian
        res = nlminb(&x0, &mut prob, &[lower], &[upper], Some(200), false).map_err(|e| e.to_string())?;
        if res.convergence != 0 {
            // "Still problematic result: do the same thing without Cox-Reid adjustment"
            let mut prob = Profile { y, mu: &mu, cr: false, tab: &tab, ones };
            res = nlminb(&x0, &mut prob, &[lower], &[upper], Some(200), false).map_err(|e| e.to_string())?;
        }
    }
    Ok(Overdispersion { estimate: exp(res.par[0]), iterations: res.iterations, convergence: Some(res.convergence) })
}

// --- the Newton fit -----------------------------------------------------------------------------

/// glmGamPoi's `compute_gp_deviance(y, mu, theta)`.
pub fn gp_deviance(y: f64, mu: f64, theta: f64) -> f64 {
    let max0 = |a: f64| if a < 0.0 { 0.0 } else { a }; // std::max(a, 0.0)
    if theta < 1e-6 {
        if y == 0.0 {
            2.0 * mu
        } else {
            max0(2.0 * (y * log(y / mu) - (y - mu)))
        }
    } else if y == 0.0 {
        2.0 / theta * log(1.0 + mu * theta)
    } else {
        let s1 = y * log((mu + y * mu * theta) / (y + y * mu * theta));
        let s2 = 1.0 / theta * log((1.0 + mu * theta) / (1.0 + y * theta));
        max0(-2.0 * (s1 - s2))
    }
}

/// The result of `fitBeta_one_group` for one gene.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaFit {
    pub beta: f64,
    pub iterations: i32,
    pub deviance: f64,
}

/// `fitBeta_one_group` for one gene: Newton-Raphson on the intercept, tolerance 1e-8, at most 100
/// iterations; when a step makes beta NaN or the iterations run out, [`optimize_helper`] (the
/// reported iterations stay Newton's). Errors only where that fallback is not ported (theta = 0).
pub fn fit_beta_one_group(counts: &[f64], off: &[f64], theta: f64, beta_start: f64) -> Result<BetaFit, String> {
    const TOL: f64 = 1e-8;
    const MAX_ITER: i32 = 100;
    if beta_start.is_nan() || theta.is_nan() {
        return Ok(BetaFit { beta: rcompat::NA_REAL, iterations: 0, deviance: rcompat::NA_REAL });
    }
    let mut beta = beta_start;
    let mut iter = 0;
    while iter < MAX_ITER {
        let (mut dl, mut ddl) = (0.0f64, 0.0f64);
        let mut all_zero = true;
        for (&count, &o) in counts.iter().zip(off) {
            all_zero = all_zero && count == 0.0;
            let mu = exp(beta + o);
            let denom = 1.0 + mu * theta;
            dl += (count - mu) / denom;
            ddl += mu * (1.0 + count * theta) / denom / denom;
        }
        if all_zero {
            beta = f64::NEG_INFINITY;
            break;
        }
        let step = dl / ddl;
        beta += step;
        if step.abs() < TOL {
            break;
        } else if beta.is_nan() {
            break;
        }
        iter += 1;
    }
    if iter == MAX_ITER || beta.is_nan() {
        // not converged: glmGamPoi tries again with optimize() and still reports `iter`
        beta = optimize_helper(counts, off, theta).map_err(|e| format!("fitBeta_one_group (theta {theta:e}, start {beta_start:e}, {iter} Newton iterations): {e}"))?;
    }
    let mut dev = 0.0;
    for (&count, &o) in counts.iter().zip(off) {
        dev += gp_deviance(count, exp(beta + o), theta);
    }
    Ok(BetaFit { beta, iterations: iter, deviance: dev })
}

/// glmGamPoi's `estimate_betas_group_wise_optimize_helper(y, offset, theta)`, the fallback of
/// [`fit_beta_one_group`]:
/// `optimize(function(beta) sum(dnbinom(y, mu = exp(beta + offset), size = 1/theta, log = TRUE)),
/// lower = -30, upper = 30, maximum = TRUE)$maximum`. With a half-integer count (a doublet)
/// `dnbinom` is -Inf for every beta, so `optimize()` sees a constant (`DBL_MAX` after negation)
/// and returns where its golden-section steps end, as in R. `theta = 0` (`size = Inf`, the
/// Poisson limit, which needs nmath's `dpois_raw`) is not ported and gives an error.
pub fn optimize_helper(counts: &[f64], off: &[f64], theta: f64) -> Result<f64, String> {
    let size = 1.0 / theta;
    if !size.is_finite() {
        return Err(format!("optimize() fallback with theta = {theta:e}: dnbinom's Poisson limit (size = Inf) is not ported"));
    }
    let mut lp = vec![0.0f64; counts.len()];
    let mut loglik = |beta: f64| -> f64 {
        for ((l, &y), &o) in lp.iter_mut().zip(counts).zip(off) {
            *l = dnbinom_mu_log(y, size, exp(beta + o));
        }
        rcompat::sum(&lp, false)
    };
    Ok(optimize_maximum(&mut loglik, -30.0, 30.0, OPTIMIZE_TOL))
}

// --- the dispersion trend -------------------------------------------------------------------------

/// glmGamPoi's `loc_median_fit(x, y, npoints, weighted = TRUE)`.
pub fn loc_median_fit(x: &[f64], y: &[f64], npoints: f64) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let npoints = if 1.0 > npoints { 1.0 } else { npoints }; // max(1, npoints)
    let npoints = if (n as f64) < npoints { n as f64 } else { npoints }; // min(length(x), npoints)
    let ord = order(x, false);
    let ordered_y: Vec<f64> = ord.iter().map(|&i| y[i]).collect();
    let half_points = (npoints / 2.0).floor() as usize;
    let start = half_points + 1;
    let end = n as i64 - half_points as i64;
    if end < start as i64 {
        let w: Vec<f64> = seq_length_out(-3.0, 3.0, n as f64).unwrap().iter().map(|&v| dnorm(v, 0.0, 1.0, false)).collect();
        let wm = weighted_median_interpolated(&ordered_y, &w); // weightedMedian(ordered_y, w = dnorm(...)): interpolate = TRUE
        return vec![wm; n];
    }
    let end = end as usize;
    let weights: Vec<f64> = seq_length_out(-3.0, 3.0, (half_points * 2 + 1) as f64).unwrap().iter().map(|&v| dnorm(v, 0.0, 1.0, false)).collect();
    let mut res = vec![f64::NAN; n];
    for idx in start..=end {
        // 1-based window idx - half_points ..= idx + half_points
        let sel = &ordered_y[idx - half_points - 1..idx + half_points];
        res[idx - 1] = weighted_median(sel, &weights);
    }
    let first = res[start - 1];
    for r in res.iter_mut().take(start.saturating_sub(1).max(1)) {
        *r = first;
    }
    let last = res[end - 1];
    for r in res.iter_mut().skip(n.min(end + 1) - 1) {
        *r = last;
    }
    // res[order(order(x))]
    let oo = order(&ord.iter().map(|&i| i as f64).collect::<Vec<_>>(), false);
    oo.iter().map(|&i| res[i]).collect()
}

// --- glm_gp for one bin -------------------------------------------------------------------------

/// Per-gene results of [`glm_gp_offset`], with the intermediate values the R hooks dump
/// (`S2b_glm`).
#[derive(Debug, Clone, PartialEq)]
pub struct GeneGlm {
    pub disp_init: f64,
    pub beta_init: f64,
    pub beta1: f64,
    pub iter1: i32,
    pub disp_mle: f64,
    pub disp_trend: f64,
    pub beta2: f64,
    pub iter2: i32,
    pub deviance2: f64,
}

/// `glm_gp(Y, design = ~ 1, offset = offset, size_factors = FALSE)` for one bin: `y` holds the
/// genes' counts gene by gene (`ngenes x ncells`, row-major: gene `g` at `g * ncells`), `offset`
/// one value per cell. `threads` parallelises the per-gene work (results do not depend on it).
pub fn glm_gp_offset(y: &[f64], ngenes: usize, ncells: usize, offset: &[f64], threads: usize) -> Result<Vec<GeneGlm>, String> {
    assert_eq!(y.len(), ngenes * ncells);
    assert_eq!(offset.len(), ncells);
    let row = |g: usize| &y[g * ncells..(g + 1) * ncells];
    let exp_off: Vec<f64> = offset.iter().map(|&o| exp(o)).collect();
    let ones = vec![1.0f64; ncells];
    // 1. dispersions by moments: xim = 1 / mean(colMeans2(exp(offset_matrix)))
    let col_means: Vec<f64> = exp_off.iter().map(|&v| col_mean_of_constant(v, ngenes)).collect();
    let xim = 1.0 / mean(&col_means, false);
    // per gene: rowVars, rowMeans2, rough beta
    let per_gene: Vec<(f64, f64)> = par_map(ngenes, threads, |g| {
        let r = row(g);
        let bv = var_refined(r);
        let bm = mean_refined(r.iter().copied());
        let d = (bv - xim * bm) / (bm * bm);
        let disp_init = if d.is_nan() || d < 0.0 { 0.0 } else { d };
        // 2. log(rowMeans2(Y / exp(offset)))
        let beta_init = log(mean_refined(r.iter().zip(&exp_off).map(|(&c, &e)| c / e)));
        (disp_init, beta_init)
    });
    // 3. first Newton fit
    let fit1: Vec<Result<BetaFit, String>> = par_map(ngenes, threads, |g| fit_beta_one_group(row(g), offset, per_gene[g].0, per_gene[g].1));
    let mut beta1 = Vec::with_capacity(ngenes);
    let mut iter1 = Vec::with_capacity(ngenes);
    for (g, f) in fit1.into_iter().enumerate() {
        let f = f.map_err(|e| format!("first fit, gene {} of the bin: {e}", g + 1))?;
        // pmax(Beta, -1e8)
        let b = if !f.beta.is_nan() && -1e8 > f.beta { -1e8 } else { f.beta };
        beta1.push(b);
        iter1.push(f.iterations);
    }
    // Mu = exp(Beta %*% t(model_matrix) + offset); the overdispersion MLE; rowMeans2(Mu)
    let mle: Vec<Result<(f64, f64), String>> = par_map(ngenes, threads, |g| {
        let mu: Vec<f64> = offset.iter().map(|&o| exp(beta1[g] * 1.0 + o)).collect();
        let od = overdispersion_mle(row(g), &mu, &ones)?;
        Ok((od.estimate, mean_refined(mu.iter().copied())))
    });
    let mut disp_mle = Vec::with_capacity(ngenes);
    let mut gene_means = Vec::with_capacity(ngenes);
    for m in mle {
        let (d, gm) = m?;
        disp_mle.push(d);
        gene_means.push(gm);
    }
    // 5. the dispersion trend
    let npoints = { let v = 0.1 * ngenes as f64; if 100.0 > v { 100.0 } else { v } };
    let est: Vec<usize> = (0..ngenes).filter(|&g| !(disp_mle[g].is_nan() || gene_means[g].is_nan())).collect();
    let trend_est = loc_median_fit(&est.iter().map(|&g| gene_means[g]).collect::<Vec<_>>(), &est.iter().map(|&g| disp_mle[g]).collect::<Vec<_>>(), npoints);
    let mut trend = vec![rcompat::NA_REAL; ngenes];
    for (k, &g) in est.iter().enumerate() {
        trend[g] = trend_est[k];
    }
    // 6. second Newton fit at the trend, from the first fit's Beta (times 1, exactly)
    let fit2: Vec<Result<BetaFit, String>> = par_map(ngenes, threads, |g| fit_beta_one_group(row(g), offset, trend[g], beta1[g] * 1.0));
    let mut out = Vec::with_capacity(ngenes);
    for (g, f) in fit2.into_iter().enumerate() {
        let f = f.map_err(|e| format!("second fit, gene {} of the bin: {e}", g + 1))?;
        let b = if !f.beta.is_nan() && -1e8 > f.beta { -1e8 } else { f.beta };
        out.push(GeneGlm {
            disp_init: per_gene[g].0,
            beta_init: per_gene[g].1,
            beta1: beta1[g],
            iter1: iter1[g],
            disp_mle: disp_mle[g],
            disp_trend: trend[g],
            beta2: b,
            iter2: f.iterations,
            deviance2: f.deviance,
        });
    }
    Ok(out)
}

/// matrixStats' refined mean (double) of an iterator.
fn mean_refined<I: Iterator<Item = f64> + Clone>(it: I) -> f64 {
    let v: Vec<f64> = it.collect();
    rcompat::matrixstats::row_means2(&v, 1, v.len(), None)[0]
}

/// `colMeans2` of a column holding `k` copies of `v` (each row of the offset matrix is the
/// offset vector), computed as matrixStats does.
fn col_mean_of_constant(v: f64, k: usize) -> f64 {
    let col = vec![v; k];
    rcompat::matrixstats::col_means2(&col, k, 1)[0]
}

/// `(0..n).map(f)` on `threads` threads, results in order.
pub(crate) fn par_map<R: Send>(n: usize, threads: usize, f: impl Fn(usize) -> R + Sync) -> Vec<R> {
    let threads = threads.max(1).min(n.max(1));
    if threads == 1 {
        return (0..n).map(f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<R>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= n {
                    break;
                }
                *slots[i].lock().unwrap() = Some(f(i));
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bucket_counts() {
        assert_eq!(next_bkt(158), 167);
        assert_eq!(next_bkt(2500), 2549);
        assert_eq!(next_bkt(13), 13);
    }

    #[test]
    fn table_truncates_and_orders() {
        let t = CountTable::new(&[1.5, 1.0, 0.0, 2.0, 0.0], 2);
        // more than 2 distinct keys: no table
        assert!(t.keys.is_empty());
        let t = CountTable::new(&[1.5, 1.0, 0.0, 0.5, 0.0, 1.0], 3);
        // keys 1 (x3: 1.5, 1, 1) and 0 (x3); inserted 1 then 0, each into an empty bucket:
        // each goes to the front of the list
        assert_eq!(t.keys, vec![0.0, 1.0]);
        assert_eq!(t.freqs, vec![3.0, 3.0]);
    }
}
