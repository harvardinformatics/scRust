//! R's `sample()` and `sample.int()`, ported from R 4.5.3 (`src/library/base/R/sample.R`,
//! `do_sample` in `src/main/random.c`, `do_sample2` in `src/main/unique.c`, and `revsort` in
//! `src/main/sort.c`), with `sample.kind = "Rejection"`.
//!
//! Results are 0-based: R's indices minus one. Each function draws exactly the uniforms R draws.

use crate::rng::RRng;
use std::collections::HashSet;
use std::fmt;

/// The conditions under which R's `sample()` stops with an error.
#[derive(Debug, Clone, PartialEq)]
pub enum SampleError {
    /// "invalid first argument": a population of 0 (or more than 4.5e15) with size > 0.
    InvalidPopulation,
    /// "cannot take a sample larger than the population when 'replace = FALSE'"
    TooLarge,
    /// "NA in probability vector"
    NaProbability,
    /// "negative probability"
    NegativeProbability,
    /// "too few positive probabilities"
    TooFewPositive,
}

impl fmt::Display for SampleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            SampleError::InvalidPopulation => "invalid first argument",
            SampleError::TooLarge => {
                "cannot take a sample larger than the population when 'replace = FALSE'"
            }
            SampleError::NaProbability => "NA in probability vector",
            SampleError::NegativeProbability => "negative probability",
            SampleError::TooFewPositive => "too few positive probabilities",
        })
    }
}

impl std::error::Error for SampleError {}

/// `sample.int(n, size, replace)` without weights: 0-based indices into a population of `n`.
///
/// Like R, this uses the hash-based algorithm of `sample2` when `n > 1e7`, `replace` is false and
/// `size <= n/2`, and otherwise `do_sample`: one `R_unif_index(n)` per draw with replacement (and
/// for `size < 2`), or a partial Fisher-Yates shuffle that swaps the last remaining element into
/// each drawn slot ("swap-remove").
pub fn sample_int(rng: &mut RRng, n: usize, size: usize, replace: bool) -> Result<Vec<usize>, SampleError> {
    let dn = n as f64;
    if !replace && dn > 1e7 && size as f64 <= dn / 2.0 {
        return sample2(rng, n, size);
    }
    if dn > 4.5e15 || (size > 0 && n == 0) {
        return Err(SampleError::InvalidPopulation);
    }
    if !replace && size > n {
        return Err(SampleError::TooLarge);
    }
    if replace || size < 2 {
        return Ok((0..size).map(|_| rng.unif_index(dn) as usize).collect());
    }
    let mut x: Vec<usize> = (0..n).collect();
    let mut remaining = n;
    let mut out = Vec::with_capacity(size);
    for _ in 0..size {
        let j = rng.unif_index(remaining as f64) as usize;
        out.push(x[j]);
        remaining -= 1;
        x[j] = x[remaining];
    }
    Ok(out)
}

/// `.Internal(sample2(n, size))`: draws `R_unif_index(n)` and redraws values already taken, up
/// to 100 tries per element (after which R keeps the duplicate).
fn sample2(rng: &mut RRng, n: usize, size: usize) -> Result<Vec<usize>, SampleError> {
    let dn = n as f64;
    if dn > 4.5e15 || (size > 0 && n == 0) {
        return Err(SampleError::InvalidPopulation);
    }
    let mut seen = HashSet::with_capacity(size);
    let mut out = Vec::with_capacity(size);
    for _ in 0..size {
        let mut v = 0;
        for _ in 0..100 {
            v = rng.unif_index(dn) as usize;
            if seen.insert(v) {
                break;
            }
        }
        out.push(v);
    }
    Ok(out)
}

/// `sample.int(length(prob), size, replace, prob)`: weighted sampling, 0-based indices.
///
/// The probabilities are checked and normalised as `FixupProb` does. With replacement, R uses
/// Walker's alias method when more than 200 entries have `n * p > 0.1`, and otherwise inversion
/// against the cumulative probabilities sorted by `revsort`. Without replacement, it uses
/// `ProbSampleNoReplace`.
pub fn sample_int_weighted(rng: &mut RRng, prob: &[f64], size: usize, replace: bool) -> Result<Vec<usize>, SampleError> {
    let n = prob.len();
    if size > 0 && n == 0 {
        return Err(SampleError::InvalidPopulation);
    }
    if !replace && size > n {
        return Err(SampleError::TooLarge);
    }
    let mut p = prob.to_vec();
    fixup_prob(&mut p, size, replace)?;
    let ans = if replace {
        let nc = p.iter().filter(|&&pi| n as f64 * pi > 0.1).count();
        if nc > 200 {
            walker_prob_sample_replace(rng, &p, size)
        } else {
            prob_sample_replace(rng, &mut p, size)
        }
    } else {
        prob_sample_no_replace(rng, &mut p, size)
    };
    Ok(ans.into_iter().map(|i| i - 1).collect())
}

/// `FixupProb`: reject non-finite and negative weights, require enough positive ones, and divide
/// by their sum (accumulated in double, in index order).
fn fixup_prob(p: &mut [f64], require_k: usize, replace: bool) -> Result<(), SampleError> {
    let mut sum = 0.0;
    let mut npos = 0;
    for &pi in p.iter() {
        if !pi.is_finite() {
            return Err(SampleError::NaProbability);
        }
        if pi < 0.0 {
            return Err(SampleError::NegativeProbability);
        }
        if pi > 0.0 {
            npos += 1;
            sum += pi;
        }
    }
    if npos == 0 || (!replace && require_k > npos) {
        return Err(SampleError::TooFewPositive);
    }
    for pi in p.iter_mut() {
        *pi /= sum;
    }
    Ok(())
}

/// `revsort(a, ib, n)`: sorts `a` into descending order by heapsort, permuting `ib` alongside.
/// Heapsort is not stable, so the order of tied values is part of R's behaviour.
pub fn revsort(a: &mut [f64], ib: &mut [i32]) {
    let n = a.len();
    assert_eq!(n, ib.len());
    if n <= 1 {
        return;
    }
    // The C code indexes from 1 (a--, ib--).
    let at = |i: usize| i - 1;
    let mut l = (n >> 1) + 1;
    let mut ir = n;
    loop {
        let (ra, ii);
        if l > 1 {
            l -= 1;
            ra = a[at(l)];
            ii = ib[at(l)];
        } else {
            ra = a[at(ir)];
            ii = ib[at(ir)];
            a[at(ir)] = a[at(1)];
            ib[at(ir)] = ib[at(1)];
            ir -= 1;
            if ir == 1 {
                a[at(1)] = ra;
                ib[at(1)] = ii;
                return;
            }
        }
        let mut i = l;
        let mut j = l << 1;
        while j <= ir {
            if j < ir && a[at(j)] > a[at(j + 1)] {
                j += 1;
            }
            if ra > a[at(j)] {
                a[at(i)] = a[at(j)];
                ib[at(i)] = ib[at(j)];
                i = j;
                j += i;
            } else {
                j = ir + 1;
            }
        }
        a[at(i)] = ra;
        ib[at(i)] = ii;
    }
}

/// `ProbSampleReplace`: inversion against the cumulative sums of the probabilities sorted by
/// `revsort`; one uniform per draw. Returns 1-based indices.
fn prob_sample_replace(rng: &mut RRng, p: &mut [f64], nans: usize) -> Vec<usize> {
    let n = p.len();
    let mut perm: Vec<i32> = (1..=n as i32).collect();
    revsort(p, &mut perm);
    for i in 1..n {
        p[i] += p[i - 1];
    }
    let nm1 = n - 1;
    (0..nans)
        .map(|_| {
            let ru = rng.unif_rand();
            let mut j = 0;
            while j < nm1 {
                if ru <= p[j] {
                    break;
                }
                j += 1;
            }
            perm[j] as usize
        })
        .collect()
}

/// `walker_ProbSampleReplace` (Walker's alias method, Ripley 1987, Alg 3.13B) with
/// `sample.kind = "Rejection"`: each draw takes `k = R_unif_index(n)` and then one more uniform
/// for `rU = k + u`. Returns 1-based indices.
fn walker_prob_sample_replace(rng: &mut RRng, p: &[f64], nans: usize) -> Vec<usize> {
    let n = p.len();
    // HL[0..H) holds the entries with q < 1, HL[L..n) those with q >= 1 (filled from the end).
    let mut hl = vec![0usize; n];
    let mut q = vec![0f64; n];
    // R leaves the alias of an entry it never processes uninitialised; such an entry can only be
    // reached through rounding error. Zero (R's index 1) is what fresh memory holds.
    let mut a = vec![0usize; n];
    let (mut h, mut l) = (0usize, n);
    for i in 0..n {
        q[i] = p[i] * n as f64;
        if q[i] < 1.0 {
            hl[h] = i;
            h += 1;
        } else {
            l -= 1;
            hl[l] = i;
        }
    }
    if h > 0 && l < n {
        for k in 0..n - 1 {
            let i = hl[k];
            let j = hl[l];
            a[i] = j;
            q[j] += q[i] - 1.0;
            if q[j] < 1.0 {
                l += 1;
            }
            if l >= n {
                break;
            }
        }
    }
    for (i, qi) in q.iter_mut().enumerate() {
        *qi += i as f64;
    }
    (0..nans)
        .map(|_| {
            let k = rng.unif_index(n as f64) as usize;
            let ru = k as f64 + rng.unif_rand();
            if ru < q[k] {
                k + 1
            } else {
                a[k] + 1
            }
        })
        .collect()
}

/// `ProbSampleNoReplace`: sort by `revsort`, then for each draw find the first cumulative mass
/// reaching `totalmass * u`, and remove that entry. Returns 1-based indices.
fn prob_sample_no_replace(rng: &mut RRng, p: &mut [f64], nans: usize) -> Vec<usize> {
    let n = p.len();
    let mut perm: Vec<i32> = (1..=n as i32).collect();
    revsort(p, &mut perm);
    let mut totalmass = 1.0;
    let mut out = Vec::with_capacity(nans);
    let mut n1 = n as isize - 1;
    for _ in 0..nans {
        let rt = totalmass * rng.unif_rand();
        let mut mass = 0.0;
        let mut j = 0usize;
        while (j as isize) < n1 {
            mass += p[j];
            if rt <= mass {
                break;
            }
            j += 1;
        }
        out.push(perm[j] as usize);
        totalmass -= p[j];
        let n1u = n1 as usize;
        for k in j..n1u {
            p[k] = p[k + 1];
            perm[k] = perm[k + 1];
        }
        n1 -= 1;
    }
    out
}

/// `sample(x, size, replace)` for a vector `x`: `x[sample.int(length(x), size, replace)]`.
/// (R treats a single number `x >= 1` as `1:x`; use [`sample_int`] for that.)
pub fn sample<T: Clone>(rng: &mut RRng, x: &[T], size: usize, replace: bool) -> Result<Vec<T>, SampleError> {
    Ok(sample_int(rng, x.len(), size, replace)?.into_iter().map(|i| x[i].clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revsort_sorts_descending_with_permutation() {
        let mut a = [0.3, 0.1, 0.4, 0.1, 0.5, 0.9, 0.2, 0.6];
        let orig = a;
        let mut ib: Vec<i32> = (1..=8).collect();
        revsort(&mut a, &mut ib);
        assert!(a.windows(2).all(|w| w[0] >= w[1]));
        for (v, &i) in a.iter().zip(&ib) {
            assert_eq!(*v, orig[(i - 1) as usize]);
        }
    }

    #[test]
    fn errors_like_r() {
        let mut r = RRng::new(1);
        assert_eq!(sample_int(&mut r, 3, 4, false), Err(SampleError::TooLarge));
        assert_eq!(sample_int(&mut r, 0, 1, true), Err(SampleError::InvalidPopulation));
        assert_eq!(sample_int_weighted(&mut r, &[1.0, 0.0], 2, false), Err(SampleError::TooFewPositive));
        assert_eq!(sample_int_weighted(&mut r, &[1.0, -1.0], 1, true), Err(SampleError::NegativeProbability));
        assert_eq!(sample_int_weighted(&mut r, &[1.0, f64::NAN], 1, true), Err(SampleError::NaProbability));
        assert_eq!(sample_int(&mut r, 0, 0, false), Ok(vec![]));
    }

    #[test]
    fn without_replacement_is_a_permutation() {
        let mut r = RRng::new(5);
        let mut s = sample_int(&mut r, 1000, 1000, false).unwrap();
        s.sort();
        assert_eq!(s, (0..1000).collect::<Vec<_>>());
    }
}
