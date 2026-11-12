//! `nlminb()`: R 4.5.3's driver (`src/library/stats/R/nlminb.R` and `port.c`) and the PORT library
//! routines it reaches for a bounded problem with an analytic gradient, with and without an
//! analytic Hessian (`DRMNHB`, `DRMNGB` and their call trees in
//! `src/library/stats/src/portsrc.f`), ported statement by statement.
//!
//! glmGamPoi 1.22.0 maximises its overdispersion likelihood with
//! `nlminb(start, objective, gradient, hessian, lower = log(1e-16), upper = log(1e16),
//! control = list(iter.max = 200))` (`DRMNHB`) and, when that does not converge, again without
//! the Hessian (`DRMNGB`, BFGS updates and a double-dogleg step). Both are here, with their call
//! trees.
//!
//! **Layout.** PORT works through reverse communication on two arrays, `IV` (integers) and `V`
//! (doubles), and passes sub-arrays of `V` to its subroutines. The port keeps one 1-based double
//! workspace `w` holding `V(1..LV)` followed by the caller's `B`, `D`, `G`, `H` and `X`, and one
//! 1-based integer workspace `iw` holding `IV(1..LIV)`; every routine takes these and the 1-based
//! positions of its array arguments, so Fortran's aliasing of sub-arrays is reproduced exactly.
//! Labels become loops over a state enum where the Fortran jumps backwards.
//!
//! **Arithmetic.** gfortran without `-ffast-math` and without FMA (`-march=nocona`) evaluates every
//! expression left to right as written, so plain `f64` code in the same order gives the same bits;
//! `X**2` is `X*X`; `DD7TPR` and `DV2NRM` are R's C replacements, OpenBLAS `ddot` and `dnrm2`,
//! called through [`crate::blas`].

use crate::blas;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum PortError {
    /// A code path of PORT that this port does not implement.
    NotImplemented(&'static str),
    /// R's driver stops: "NA/NaN gradient evaluation", "NA/NaN Hessian evaluation".
    NanDerivative(&'static str),
}

impl fmt::Display for PortError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PortError::NotImplemented(w) => write!(f, "nlminb: {w} is not implemented"),
            PortError::NanDerivative(w) => write!(f, "NA/NaN {w} evaluation"),
        }
    }
}

impl std::error::Error for PortError {}

// --- machine constants ------------------------------------------------------------------------

/// `DR7MDC(K)` through R's `d1mach`.
fn dr7mdc(k: i32) -> f64 {
    let big = f64::MAX;
    let eta = f64::MIN_POSITIVE;
    let machep = f64::EPSILON;
    match k {
        1 => eta,
        2 => (256.0 * eta).sqrt() / 16.0,
        3 => machep,
        4 => machep.sqrt(),
        5 => (big / 256.0).sqrt() * 16.0,
        6 => big,
        _ => 0.0,
    }
}

// --- IV and V subscripts ----------------------------------------------------------------------

mod ivx {
    pub const TOOBIG: usize = 2;
    pub const IVNEED: usize = 3;
    pub const VNEED: usize = 4;
    pub const MODEL: usize = 5;
    pub const NFCALL: usize = 6;
    pub const NFGCAL: usize = 7;
    pub const RADINC: usize = 8;
    pub const RESTOR: usize = 9;
    pub const STAGE: usize = 10;
    pub const STGLIM: usize = 11;
    pub const SWITCH: usize = 12;
    pub const XIRC: usize = 13;
    pub const DTYPE: usize = 16;
    pub const MXFCAL: usize = 17;
    pub const MXITER: usize = 18;
    pub const OUTLEV: usize = 19;
    pub const PARPRT: usize = 20;
    pub const PRUNIT: usize = 21;
    pub const SOLPRT: usize = 22;
    pub const STATPR: usize = 23;
    pub const X0PRT: usize = 24;
    pub const INITS: usize = 25;
    pub const IRC: usize = 29;
    pub const NGCALL: usize = 30;
    pub const NITER: usize = 31;
    pub const MLSTGD: usize = 32;
    pub const KAGQT: usize = 33;
    pub const MODE: usize = 35;
    pub const DG: usize = 37;
    pub const OLDN: usize = 38;
    pub const STEP: usize = 40;
    pub const N0: usize = 41;
    pub const LMAT: usize = 42;
    pub const X0: usize = 43;
    pub const LASTIV: usize = 44;
    pub const LASTV: usize = 45;
    pub const NEXTIV: usize = 46;
    pub const NEXTV: usize = 47;
    pub const NC: usize = 48;
    pub const PARSAV: usize = 49;
    pub const NVDFLT: usize = 50;
    pub const ALGSAV: usize = 51;
    pub const NFCOV: usize = 52;
    pub const NGCOV: usize = 53;
    pub const DTYPE0: usize = 54;
    pub const CNVCOD: usize = 55;
    pub const PERM: usize = 58;
    pub const DTOL: usize = 59;
    pub const W: usize = 34;
}

mod vx {
    pub const DGNORM: usize = 1;
    pub const DSTNRM: usize = 2;
    pub const DST0: usize = 3;
    pub const GTSTEP: usize = 4;
    pub const STPPAR: usize = 5;
    pub const NREDUC: usize = 6;
    pub const PREDUC: usize = 7;
    pub const RADIUS: usize = 8;
    pub const RAD0: usize = 9;
    pub const F: usize = 10;
    pub const FDIF: usize = 11;
    pub const FLSTGD: usize = 12;
    pub const F0: usize = 13;
    pub const GTSLST: usize = 14;
    pub const PLSTGD: usize = 15;
    pub const RADFAC: usize = 16;
    pub const RELDX: usize = 17;
    pub const DSTSAV: usize = 18;
    pub const EPSLON: usize = 19;
    pub const PHMNFC: usize = 20;
    pub const PHMXFC: usize = 21;
    pub const DECFAC: usize = 22;
    pub const INCFAC: usize = 23;
    pub const RDFCMN: usize = 24;
    pub const RDFCMX: usize = 25;
    pub const TUNER1: usize = 26;
    pub const TUNER2: usize = 27;
    pub const TUNER3: usize = 28;
    pub const TUNER4: usize = 29;
    pub const TUNER5: usize = 30;
    pub const AFCTOL: usize = 31;
    pub const RFCTOL: usize = 32;
    pub const XCTOL: usize = 33;
    pub const XFTOL: usize = 34;
    pub const LMAX0: usize = 35;
    pub const LMAXS: usize = 36;
    pub const SCTOL: usize = 37;
    pub const DINIT: usize = 38;
    pub const DTINIT: usize = 39;
    pub const D0INIT: usize = 40;
    pub const DFAC: usize = 41;
    pub const ETA0: usize = 42;
    pub const BIAS: usize = 43;
}

use ivx::*;

// --- vector helpers (R's C replacements and short Fortran routines) ----------------------------

/// `DD7TPR(P, X, Y)`: OpenBLAS `ddot`.
fn dd7tpr(p: i32, w: &[f64], x: usize, y: usize) -> f64 {
    if p <= 0 {
        // ddot with n <= 0 returns 0
        return 0.0;
    }
    let p = p as usize;
    blas::ddot(&w[x..x + p], &w[y..y + p])
}

/// `DV2NRM(P, X)`: OpenBLAS `dnrm2`.
fn dv2nrm(p: i32, w: &[f64], x: usize) -> f64 {
    if p <= 0 {
        return 0.0;
    }
    blas::dnrm2(&w[x..x + p as usize])
}

/// `DV2AXY(P, W, A, X, Y)`: `W = A*X + Y`, element by element in order (the C replacement).
fn dv2axy(p: i32, w: &mut [f64], out: usize, a: f64, x: usize, y: usize) {
    for i in 0..p.max(0) as usize {
        w[out + i] = a * w[x + i] + w[y + i];
    }
}

/// `DV7CPY(P, Y, X)`: copy (memmove).
fn dv7cpy(p: i32, w: &mut [f64], dest: usize, src: usize) {
    if p > 0 {
        w.copy_within(src..src + p as usize, dest);
    }
}

/// `DV7SCP(P, Y, S)`: set to a constant.
fn dv7scp(p: i32, w: &mut [f64], dest: usize, c: f64) {
    for i in 0..p.max(0) as usize {
        w[dest + i] = c;
    }
}

/// `DV7VMP(N, X, Y, Z, K)`: `X = Y * Z` (K >= 0) or `Y / Z` (K < 0).
fn dv7vmp(n: i32, w: &mut [f64], x: usize, y: usize, z: usize, k: i32) {
    for i in 0..n.max(0) as usize {
        w[x + i] = if k >= 0 { w[y + i] * w[z + i] } else { w[y + i] / w[z + i] };
    }
}

/// `DV7IPR(N, IP, X)`: `X(i) := X(IP(i))` (the C replacement; IP 1-based).
fn dv7ipr(n: i32, iw: &[i32], ip: usize, w: &mut [f64], x: usize) {
    let n = n.max(0) as usize;
    let xcp: Vec<f64> = (0..n).map(|i| w[x + iw[ip + i] as usize - 1]).collect();
    w[x..x + n].copy_from_slice(&xcp);
}

/// `I7PNVR(N, X, Y)`: inverse permutation, `X(Y(i)) = i`.
fn i7pnvr(n: i32, iw: &mut [i32], x: usize, y: usize) {
    for i in 0..n.max(0) as usize {
        let yi = iw[y + i] as usize;
        iw[x + yi - 1] = i as i32 + 1;
    }
}

/// `I7SHFT(N, K, X)`: circular shift.
fn i7shft(n: i32, k: i32, iw: &mut [i32], x: usize) {
    let at = |i: i32| x + i as usize - 1;
    if k < 0 {
        let k1 = -k;
        if k1 >= n {
            return;
        }
        let t = iw[at(n)];
        let nm1 = n - k1;
        for ii in 1..=nm1 {
            let i = n - ii;
            iw[at(i + 1)] = iw[at(i)];
        }
        iw[at(k1)] = t;
    } else {
        if k >= n {
            return;
        }
        let nm1 = n - 1;
        let t = iw[at(k)];
        for i in k..=nm1 {
            iw[at(i)] = iw[at(i + 1)];
        }
        iw[at(n)] = t;
    }
}

/// `DV7SHF(N, K, X)`: shift `X(K..N)` left by one, `X(K)` to the end.
fn dv7shf(n: i32, k: i32, w: &mut [f64], x: usize) {
    if k >= n {
        return;
    }
    let at = |i: i32| x + i as usize - 1;
    let t = w[at(k)];
    for i in k..=n - 1 {
        w[at(i)] = w[at(i + 1)];
    }
    w[at(n)] = t;
}

/// `DH2RFG(A, B, X, Y, Z)`: 2x2 Householder reflection; returns the new A.
fn dh2rfg(a: f64, b: f64) -> (f64, f64, f64, f64) {
    if b == 0.0 {
        return (a, 0.0, 0.0, 0.0);
    }
    let t = a.abs() + b.abs();
    let mut a1 = a / t;
    let b1 = b / t;
    let mut c = (a1 * a1 + b1 * b1).sqrt();
    if a1 > 0.0 {
        c = -c;
    }
    a1 -= c;
    let z = b1 / a1;
    let x = a1 / c;
    let y = b1 / c;
    (t * c, x, y, z)
}

/// `DH2RFA(N, A, B, X, Y, Z)` for N = 1 (the only use here).
fn dh2rfa1(w: &mut [f64], a: usize, b: usize, x: f64, y: f64, z: f64) {
    let t = w[a] * x + w[b] * y;
    w[a] += t;
    w[b] += t * z;
}

/// `DQ7RSH(K, P, HAVQTR, QTR, R, W)` with `HAVQTR = .FALSE.`.
fn dq7rsh(k: i32, p: i32, w: &mut [f64], r: usize, wv: usize) {
    if k >= p {
        return;
    }
    let rr = |i: i32| r + i as usize - 1;
    let wa = |i: i32| wv + i as usize - 1;
    let km1 = k - 1;
    let mut k1 = k * km1 / 2;
    dv7cpy(k, w, wa(1), rr(k1 + 1));
    let mut wj = w[wa(k)];
    let pm1 = p - 1;
    let mut j1 = k1 + km1;
    for j in k..=pm1 {
        let jm1 = j - 1;
        let jp1 = j + 1;
        if jm1 > 0 {
            dv7cpy(jm1, w, rr(k1 + 1), rr(j1 + 2));
        }
        j1 += jp1;
        k1 += j;
        let a = w[rr(j1)];
        let b = w[rr(j1 + 1)];
        let (x, z);
        if b == 0.0 {
            w[rr(k1)] = a;
            x = 0.0;
            z = 0.0;
        } else {
            let (newa, xx, yy, zz) = dh2rfg(a, b);
            w[rr(k1)] = newa;
            x = xx;
            z = zz;
            if j != pm1 {
                let mut i1 = j1;
                for i in jp1..=pm1 {
                    i1 += i;
                    dh2rfa1(w, rr(i1), rr(i1 + 1), xx, yy, zz);
                }
            }
        }
        let t = x * wj;
        w[wa(j)] = wj + t;
        wj = t * z;
    }
    w[wa(p)] = wj;
    dv7cpy(p, w, rr(k1 + 1), wa(1));
}

/// `DS7IPR(P, IP, H)`: apply the permutation IP to the rows and columns of the packed lower
/// triangle H.
fn ds7ipr(p: i32, iw: &mut [i32], ip: usize, w: &mut [f64], h: usize) {
    let ipa = |i: i32| ip + i as usize - 1;
    let ha = |i: i32| h + i as usize - 1;
    for i in 1..=p {
        let mut j = iw[ipa(i)];
        if j == i {
            continue;
        }
        iw[ipa(i)] = j.abs();
        if j < 0 {
            continue;
        }
        let mut k = i;
        loop {
            // 10
            let (mut j1, mut k1) = (j, k);
            if j > k {
                j1 = k;
                k1 = j;
            }
            let kmj = k1 - j1;
            let mut l = j1 - 1;
            let mut jm = j1 * l / 2;
            let mut km = k1 * (k1 - 1) / 2;
            if l > 0 {
                for _ in 1..=l {
                    jm += 1;
                    let t = w[ha(jm)];
                    km += 1;
                    w[ha(jm)] = w[ha(km)];
                    w[ha(km)] = t;
                }
            }
            km += 1;
            let kk = km + kmj;
            jm += 1;
            let t = w[ha(jm)];
            w[ha(jm)] = w[ha(kk)];
            w[ha(kk)] = t;
            j1 = l;
            l = kmj - 1;
            if l > 0 {
                for m in 1..=l {
                    jm = jm + j1 + m;
                    let t = w[ha(jm)];
                    km += 1;
                    w[ha(jm)] = w[ha(km)];
                    w[ha(km)] = t;
                }
            }
            if k1 < p {
                let l = p - k1;
                let k1m = k1 - 1;
                let mut km = kk;
                for m in 1..=l {
                    km = km + k1m + m;
                    let jm = km - kmj;
                    let t = w[ha(jm)];
                    w[ha(jm)] = w[ha(km)];
                    w[ha(km)] = t;
                }
            }
            // 80
            k = j;
            j = iw[ipa(k)];
            iw[ipa(k)] = -j;
            if j <= i {
                break;
            }
        }
    }
}

/// `DL7IVM(N, X, L, Y)`: solve `L*X = Y` (L packed lower triangular, by rows).
fn dl7ivm(n: i32, w: &mut [f64], x: usize, l: usize, y: usize) {
    let xa = |i: i32| x + i as usize - 1;
    let la = |i: i32| l + i as usize - 1;
    let ya = |i: i32| y + i as usize - 1;
    let mut k = 1;
    while k <= n {
        if w[ya(k)] != 0.0 {
            break;
        }
        w[xa(k)] = 0.0;
        k += 1;
    }
    if k > n {
        return;
    }
    let mut j = k * (k + 1) / 2;
    w[xa(k)] = w[ya(k)] / w[la(j)];
    if k >= n {
        return;
    }
    let k = k + 1;
    for i in k..=n {
        let t = dd7tpr(i - 1, w, la(j + 1), xa(1));
        j += i;
        w[xa(i)] = (w[ya(i)] - t) / w[la(j)];
    }
}

/// `DL7ITV(N, X, L, Y)`: solve `L**T * X = Y`.
fn dl7itv(n: i32, w: &mut [f64], x: usize, l: usize, y: usize) {
    let xa = |i: i32| x + i as usize - 1;
    let la = |i: i32| l + i as usize - 1;
    for i in 1..=n {
        w[xa(i)] = w[y + i as usize - 1];
    }
    let np1 = n + 1;
    let mut i0 = n * (n + 1) / 2;
    for ii in 1..=n {
        let i = np1 - ii;
        let xi = w[xa(i)] / w[la(i0)];
        w[xa(i)] = xi;
        if i <= 1 {
            return;
        }
        i0 -= i;
        if xi == 0.0 {
            continue;
        }
        let im1 = i - 1;
        for j in 1..=im1 {
            let ij = i0 + j;
            w[xa(j)] -= xi * w[la(ij)];
        }
    }
}

/// `DL7SRT(N1, N, L, A, IRC)`: Cholesky factor rows N1..N of A into L; returns IRC.
fn dl7srt(n1: i32, n: i32, w: &mut [f64], l: usize, a: usize) -> i32 {
    let la = |i: i32| l + i as usize - 1;
    let aa = |i: i32| a + i as usize - 1;
    let mut i0 = n1 * (n1 - 1) / 2;
    for i in n1..=n {
        let mut td = 0.0;
        if i != 1 {
            let mut j0 = 0;
            let im1 = i - 1;
            for j in 1..=im1 {
                let mut t = 0.0;
                if j != 1 {
                    let jm1 = j - 1;
                    for k in 1..=jm1 {
                        let ik = i0 + k;
                        let jk = j0 + k;
                        t += w[la(ik)] * w[la(jk)];
                    }
                }
                let ij = i0 + j;
                j0 += j;
                t = (w[aa(ij)] - t) / w[la(j0)];
                w[la(ij)] = t;
                td += t * t;
            }
        }
        i0 += i;
        let t = w[aa(i0)] - td;
        if t <= 0.0 {
            w[la(i0)] = t;
            return i;
        }
        w[la(i0)] = t.sqrt();
    }
    0
}

/// `DL7SVN(P, L, X, Y)`: estimate of the smallest singular value of the packed lower triangle L.
fn dl7svn(p: i32, w: &mut [f64], l: usize, x: usize, y: usize) -> f64 {
    let la = |i: i32| l + i as usize - 1;
    let xa = |i: i32| x + i as usize - 1;
    let ya = |i: i32| y + i as usize - 1;
    let mut ix = 2;
    let pm1 = p - 1;
    let mut ii = 0;
    let mut j0 = p * pm1 / 2;
    let mut jj = j0 + p;
    if w[la(jj)] == 0.0 {
        return 0.0;
    }
    ix = (3432 * ix) % 9973;
    let mut b = 0.5 * (1.0 + ix as f64 / 9973.0);
    let mut xplus = b / w[la(jj)];
    w[xa(p)] = xplus;
    if p > 1 {
        for i in 1..=pm1 {
            ii += i;
            if w[la(ii)] == 0.0 {
                return 0.0;
            }
            let ji = j0 + i;
            w[xa(i)] = xplus * w[la(ji)];
        }
        for jjj in 1..=pm1 {
            let j = p - jjj;
            ix = (3432 * ix) % 9973;
            b = 0.5 * (1.0 + ix as f64 / 9973.0);
            xplus = b - w[xa(j)];
            let mut xminus = -b - w[xa(j)];
            let mut splus = xplus.abs();
            let mut sminus = xminus.abs();
            let jm1 = j - 1;
            j0 = j * jm1 / 2;
            jj = j0 + j;
            xplus /= w[la(jj)];
            xminus /= w[la(jj)];
            if jm1 != 0 {
                for i in 1..=jm1 {
                    let ji = j0 + i;
                    splus += (w[xa(i)] + w[la(ji)] * xplus).abs();
                    sminus += (w[xa(i)] + w[la(ji)] * xminus).abs();
                }
            }
            if sminus > splus {
                xplus = xminus;
            }
            w[xa(j)] = xplus;
            if jm1 > 0 {
                dv2axy(jm1, w, xa(1), xplus, la(j0 + 1), xa(1));
            }
        }
    }
    let t = 1.0 / dv2nrm(p, w, xa(1));
    for i in 1..=p {
        w[xa(i)] *= t;
    }
    for j in 1..=p {
        let jm1 = j - 1;
        let j0 = j * jm1 / 2;
        let jj = j0 + j;
        let mut t = 0.0;
        if jm1 > 0 {
            t = dd7tpr(jm1, w, la(j0 + 1), ya(1));
        }
        w[ya(j)] = (w[xa(j)] - t) / w[la(jj)];
    }
    1.0 / dv2nrm(p, w, ya(1))
}

/// `DS7LVM(P, Y, S, X)`: `Y = S * X` for the packed lower triangle of the symmetric S.
fn ds7lvm(p: i32, w: &mut [f64], y: usize, s: usize, x: usize) {
    let ya = |i: i32| y + i as usize - 1;
    let sa = |i: i32| s + i as usize - 1;
    let xa = |i: i32| x + i as usize - 1;
    let mut j = 1;
    for i in 1..=p {
        w[ya(i)] = dd7tpr(i, w, sa(j), xa(1));
        j += i;
    }
    if p <= 1 {
        return;
    }
    j = 1;
    for i in 2..=p {
        let xi = w[xa(i)];
        let im1 = i - 1;
        j += 1;
        for k in 1..=im1 {
            w[ya(k)] += w[sa(j)] * xi;
            j += 1;
        }
    }
}

/// `DRLDST(P, D, X, X0)`: relative difference between X and X0.
fn drldst(p: i32, w: &[f64], d: usize, x: usize, x0: usize) -> f64 {
    let mut emax = 0.0;
    let mut xmax = 0.0;
    for i in 0..p.max(0) as usize {
        let t = (w[d + i] * (w[x + i] - w[x0 + i])).abs();
        if emax < t {
            emax = t;
        }
        let t = w[d + i] * (w[x + i].abs() + w[x0 + i].abs());
        if xmax < t {
            xmax = t;
        }
    }
    if xmax > 0.0 { emax / xmax } else { 0.0 }
}

/// Fortran's `DMAX1(A, B)` / `DMIN1(A, B)` (gfortran's `MAX`/`MIN` for non-NaN arguments).
fn dmax1(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}
fn dmin1(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

// --- defaults ---------------------------------------------------------------------------------

/// `DV7DFL(ALG, LV, V)`: default values of V; `v` is the 1-based position of `V(1)`.
fn dv7dfl(alg: i32, w: &mut [f64], v: usize) {
    use vx::*;
    let at = |i: usize| v + i - 1;
    let machep = dr7mdc(3);
    w[at(AFCTOL)] = if machep > 1.0e-10 { machep * machep } else { 1.0e-20 };
    w[at(DECFAC)] = 0.5;
    let sqteps = dr7mdc(4);
    w[at(DFAC)] = 0.6;
    w[at(DTINIT)] = 1.0e-6;
    let mepcrt = crate::libm::pow(machep, 1.0 / 3.0);
    w[at(D0INIT)] = 1.0;
    w[at(EPSLON)] = 0.1;
    w[at(INCFAC)] = 2.0;
    w[at(LMAX0)] = 1.0;
    w[at(LMAXS)] = 1.0;
    w[at(PHMNFC)] = -0.1;
    w[at(PHMXFC)] = 0.1;
    w[at(RDFCMN)] = 0.1;
    w[at(RDFCMX)] = 4.0;
    w[at(RFCTOL)] = dmax1(1.0e-10, mepcrt * mepcrt);
    w[at(SCTOL)] = w[at(RFCTOL)];
    w[at(TUNER1)] = 0.1;
    w[at(TUNER2)] = 1.0e-4;
    w[at(TUNER3)] = 0.75;
    w[at(TUNER4)] = 0.5;
    w[at(TUNER5)] = 0.75;
    w[at(XCTOL)] = sqteps;
    w[at(XFTOL)] = 1.0e+2 * machep;
    if alg == 1 {
        // regression constants (nls); not used by nlminb
        w[at(47)] = dmax1(1.0e-6, 1.0e+2 * machep); // COSMIN
        w[at(DINIT)] = 0.0;
        w[at(44)] = sqteps; // DELTA0
        w[at(42)] = mepcrt; // DLTFDC
        w[at(43)] = sqteps; // DLTFDJ
        w[at(45)] = 1.5; // FUZZ
        w[at(46)] = dr7mdc(5); // RLIMIT
        w[at(49)] = 1.0e-3; // RSPTOL
        w[at(50)] = 1.0e-4; // SIGMIN
    } else {
        w[at(BIAS)] = 0.8;
        w[at(DINIT)] = -1.0;
        w[at(ETA0)] = 1.0e+3 * machep;
    }
}

/// R's `Rf_divset(alg, iv, liv, lv, v)` (port.c) for alg = 2 (nlminb).
fn divset(alg: i32, iw: &mut [i32], liv: usize, lv: usize, w: &mut [f64]) {
    const MINIV: [i32; 5] = [0, 82, 59, 103, 103];
    const MINV: [i32; 5] = [0, 98, 71, 101, 85];
    if PRUNIT <= liv {
        iw[PRUNIT] = 0;
    }
    if ALGSAV <= liv {
        iw[ALGSAV] = alg;
    }
    let miv = MINIV[alg as usize];
    if (liv as i32) < miv {
        iw[1] = 15;
        return;
    }
    let mv = MINV[alg as usize];
    if (lv as i32) < mv {
        iw[1] = 16;
        return;
    }
    let alg1 = (alg - 1) % 2 + 1;
    dv7dfl(alg1, w, 1);
    iw[1] = 12;
    iw[IVNEED] = 0;
    iw[LASTIV] = miv;
    iw[LASTV] = mv;
    iw[LMAT] = mv + 1;
    iw[MXFCAL] = 200;
    iw[MXITER] = 150;
    iw[OUTLEV] = 0;
    iw[PARPRT] = 1;
    iw[PERM] = miv + 1;
    iw[SOLPRT] = 0;
    iw[STATPR] = 0;
    iw[VNEED] = 0;
    iw[X0PRT] = 1;
    if alg1 >= 2 {
        iw[DTYPE] = 0;
        iw[INITS] = 1;
        iw[NFCOV] = 0;
        iw[NGCOV] = 0;
        iw[NVDFLT] = 25;
        iw[PARSAV] = if alg > 2 { 61 } else { 47 };
        w[vx::AFCTOL] = 0.0;
    }
}

/// `DPARCK(ALG, D, IV, LIV, LV, N, V)`: checks IV and V (no printing: R sets PRUNIT = 0).
fn dparck(alg: i32, d: usize, iw: &mut [i32], liv: usize, lv: usize, n: i32, w: &mut [f64]) {
    use vx::*;
    const VM: [f64; 35] = [0.0, 1.0e-3, -0.99, 1.0e-3, 1.0e-2, 1.2, 1.0e-2, 1.2, 0.0, 0.0, 1.0e-3, -1.0, 0.0 /*12*/, 0.0, 0.0 /*14*/, 0.0, 0.0, 0.0 /*17*/, 0.0 /*18*/, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0 /*24*/, 0.0 /*25*/, 0.0 /*26*/, 1.01, 1.0e10, 0.0 /*29*/, 0.0, 0.0, 0.0, 0.0 /*33*/, 0.0];
    const VX: [f64; 35] = [0.0, 0.9, -1.0e-3, 1.0e1, 0.8, 1.0e2, 0.8, 1.0e2, 0.5, 0.5, 1.0, 1.0, 0.0 /*12*/, 0.0 /*13*/, 0.1, 1.0, 1.0, 0.0 /*17*/, 0.0 /*18*/, 1.0, 0.0 /*20*/, 0.0 /*21*/, 0.0 /*22*/, 1.0, 1.0, 1.0, 1.0, 1.0e10, 0.0 /*28*/, 1.0, 0.0 /*30*/, 1.0, 1.0, 1.0, 1.0];
    const JLIM: [i32; 5] = [0, 0, 24, 0, 24];
    const NDFLT: [i32; 5] = [0, 32, 25, 32, 25];
    const MINIV: [i32; 5] = [0, 82, 59, 103, 103];
    const IJMP: usize = 33;
    let mut vm = VM;
    let mut vxx = VX;
    let pu = if PRUNIT <= liv { iw[PRUNIT] } else { 0 };
    if ALGSAV <= liv && alg != iw[ALGSAV] {
        iw[1] = 67;
        return;
    }
    if !(1..=4).contains(&alg) {
        iw[1] = 67;
        return;
    }
    let mut miv1 = MINIV[alg as usize];
    if iw[1] == 15 {
        if LASTIV <= liv {
            iw[LASTIV] = miv1;
        }
        if LASTV <= liv {
            iw[LASTV] = 0;
        }
        return;
    }
    let alg1 = (alg - 1) % 2 + 1;
    if iw[1] == 0 {
        divset(alg, iw, liv, lv, w);
    }
    let mut iv1 = iw[1];
    if iv1 == 13 || iv1 == 12 {
        if PERM <= liv {
            miv1 = miv1.max(iw[PERM] - 1);
        }
        let mut miv2 = 0;
        if IVNEED <= liv {
            miv2 = miv1 + iw[IVNEED].max(0);
        }
        if LASTIV <= liv {
            iw[LASTIV] = miv2;
        }
        if (liv as i32) < miv1 {
            iw[1] = 15;
            return;
        }
        iw[IVNEED] = 0;
        iw[LASTV] = iw[VNEED].max(0) + iw[LMAT] - 1;
        iw[VNEED] = 0;
        if (liv as i32) < miv2 {
            iw[1] = 15;
            return;
        }
        if (lv as i32) < iw[LASTV] {
            iw[1] = 16;
            return;
        }
    }
    if (12..=14).contains(&iv1) {
        if n < 1 {
            iw[1] = 81;
            return;
        }
        if iv1 != 14 {
            iw[NEXTIV] = iw[PERM];
            iw[NEXTV] = iw[LMAT];
        }
        if iv1 == 13 {
            return;
        }
        let k = (iw[PARSAV] - EPSLON as i32) as usize;
        dv7dfl(alg1, w, k + 1);
        iw[DTYPE0] = 2 - alg1;
        iw[OLDN] = n;
    } else {
        if n != iw[OLDN] {
            iw[1] = 17;
            return;
        }
        if !(1..=11).contains(&iv1) {
            iw[1] = 80;
            return;
        }
    }
    if iv1 == 14 {
        iv1 = 12;
    }
    // 110/120: machine-dependent limits
    let tiny = dr7mdc(1);
    let machep = dr7mdc(3);
    let big = dr7mdc(6);
    vm[12] = machep;
    vxx[12] = big;
    vxx[13] = big;
    vm[14] = machep;
    vm[17] = tiny;
    vxx[17] = big;
    vm[18] = tiny;
    vxx[18] = big;
    vxx[20] = big;
    vxx[21] = big;
    vxx[22] = big;
    vm[24] = machep;
    vm[25] = machep;
    vm[26] = machep;
    vxx[28] = dr7mdc(5);
    vm[29] = machep;
    vxx[30] = big;
    vm[33] = machep;
    let mut m = 0i32;
    let mut i = 1usize;
    let j = JLIM[alg1 as usize] as usize;
    let mut k = EPSLON;
    let ndfalt = NDFLT[alg1 as usize];
    for _ in 1..=ndfalt {
        let vk = w[k];
        if !(vk >= vm[i] && vk <= vxx[i]) {
            m = k as i32;
        }
        k += 1;
        i += 1;
        if i == j {
            i = IJMP;
        }
    }
    if iw[NVDFLT] != ndfalt {
        iw[1] = 51;
        return;
    }
    if !((iw[DTYPE] > 0 || w[DINIT] > 0.0) && iv1 == 12) {
        for i in 0..n as usize {
            if w[d + i] <= 0.0 {
                m = 18;
            }
        }
    }
    if m != 0 {
        iw[1] = m;
        return;
    }
    if pu == 0 || iw[PARPRT] == 0 {
        return;
    }
    unreachable!("PORT printing is disabled in R (PRUNIT = 0)");
}

// --- step assessment --------------------------------------------------------------------------

/// `DA7SST(IV, LIV, LV, V)`: assess the candidate step. Every label of the Fortran is a state.
fn da7sst(iw: &mut [i32], w: &mut [f64]) {
    use vx::*;
    #[derive(Clone, Copy)]
    enum L {
        L10, L20, L30, L40, L50, L60, L70, L80, L110, L140, L150, L160, L170, L180, L190, L200,
        L210, L220, L230, L240, L250, L260, L270, L280, L290,
    }
    let mut nfc = iw[NFCALL];
    iw[SWITCH] = 0;
    iw[RESTOR] = 0;
    let mut rfac1 = 1.0;
    let mut goodx = true;
    let mut i = iw[IRC];
    let mut emax: f64;
    let mut emaxs = 0.0f64;
    let mut lab = match i {
        1 => L::L20,
        2 => L::L30,
        3 | 4 => L::L10,
        5 => L::L40,
        6 => L::L280,
        7..=11 => L::L220,
        12 => L::L170,
        _ => {
            iw[IRC] = 13;
            return;
        }
    };
    loop {
        lab = match lab {
            L::L10 => {
                iw[STAGE] = 1;
                iw[RADINC] = 0;
                w[FLSTGD] = w[F0];
                if iw[TOOBIG] == 0 {
                    L::L110
                } else {
                    iw[STAGE] = -1;
                    iw[XIRC] = i;
                    L::L60
                }
            }
            L::L20 => {
                if iw[MODEL] != iw[MLSTGD] {
                    L::L30
                } else {
                    iw[STAGE] = iw[STGLIM];
                    iw[RADINC] = -1;
                    L::L110
                }
            }
            L::L30 => {
                iw[STAGE] += 1;
                L::L40
            }
            L::L40 => {
                if iw[STAGE] > 0 {
                    L::L50
                } else if iw[TOOBIG] != 0 {
                    L::L60
                } else {
                    iw[STAGE] = -iw[STAGE];
                    i = iw[XIRC];
                    match i {
                        1 => L::L20,
                        2 => L::L30,
                        3 | 4 => L::L110,
                        5 => L::L70,
                        // a computed GO TO out of range continues with the next statement
                        _ => L::L50,
                    }
                }
            }
            L::L50 => {
                if iw[TOOBIG] == 0 {
                    L::L70
                } else {
                    iw[TOOBIG] = 0;
                    if iw[RADINC] > 0 {
                        L::L80
                    } else {
                        iw[STAGE] = -iw[STAGE];
                        iw[XIRC] = iw[IRC];
                        L::L60
                    }
                }
            }
            L::L60 => {
                iw[TOOBIG] = 0;
                w[RADFAC] = w[DECFAC];
                iw[RADINC] -= 1;
                iw[IRC] = 5;
                iw[RESTOR] = 1;
                w[F] = w[FLSTGD];
                return;
            }
            L::L70 => {
                if w[F] < w[FLSTGD] {
                    L::L110
                } else {
                    if iw[MODEL] != iw[MLSTGD] {
                        iw[MODEL] = iw[MLSTGD];
                        iw[SWITCH] = 1;
                    }
                    L::L80
                }
            }
            L::L80 => {
                if w[FLSTGD] >= w[F0] {
                    L::L110
                } else {
                    if iw[STAGE] < iw[STGLIM] {
                        goodx = false;
                    } else if nfc < iw[NFGCAL] + iw[STGLIM] + 2 {
                        goodx = false;
                    } else if iw[SWITCH] != 0 {
                        goodx = false;
                    }
                    iw[RESTOR] = 3;
                    w[F] = w[FLSTGD];
                    w[PREDUC] = w[PLSTGD];
                    w[GTSTEP] = w[GTSLST];
                    if iw[SWITCH] == 0 {
                        rfac1 = w[DSTNRM] / w[DSTSAV];
                    }
                    w[DSTNRM] = w[DSTSAV];
                    if goodx {
                        w[FDIF] = w[F0] - w[F];
                        iw[IRC] = 4;
                        w[RADFAC] = rfac1;
                        return;
                    }
                    nfc = iw[NFGCAL];
                    L::L110
                }
            }
            L::L110 => {
                w[FDIF] = w[F0] - w[F];
                if w[FDIF] > w[TUNER2] * w[PREDUC] {
                    L::L140
                } else if iw[RADINC] > 0 {
                    L::L140
                } else {
                    if w[F] < w[F0] {
                        iw[NFGCAL] = nfc; // 120
                    } else {
                        iw[MLSTGD] = iw[MODEL];
                        w[FLSTGD] = w[F];
                        w[F] = w[F0];
                        iw[RESTOR] = 1;
                    }
                    iw[IRC] = 1; // 130
                    if iw[STAGE] >= iw[STGLIM] {
                        iw[IRC] = 5;
                        iw[RADINC] -= 1;
                    }
                    L::L160
                }
            }
            L::L140 => {
                iw[NFGCAL] = nfc;
                rfac1 = 1.0;
                w[DSTSAV] = w[DSTNRM];
                if w[FDIF] > w[PREDUC] * w[TUNER1] {
                    L::L190
                } else if iw[STAGE] >= iw[STGLIM] {
                    L::L150
                } else {
                    iw[IRC] = 2;
                    L::L160
                }
            }
            L::L150 => {
                iw[IRC] = 4;
                L::L160
            }
            L::L160 => {
                iw[XIRC] = iw[IRC];
                emax = w[GTSTEP] + w[FDIF];
                w[RADFAC] = 0.5 * rfac1;
                if emax < w[GTSTEP] {
                    w[RADFAC] = rfac1 * dmax1(w[RDFCMN], 0.5 * w[GTSTEP] / emax);
                }
                L::L170
            }
            L::L170 => {
                if w[RELDX] <= w[XFTOL] {
                    L::L180
                } else {
                    iw[IRC] = iw[XIRC];
                    if w[F] < w[F0] { L::L200 } else { L::L230 }
                }
            }
            L::L180 => {
                iw[IRC] = 12;
                L::L240
            }
            L::L190 => {
                if w[FDIF] < -w[TUNER3] * w[GTSTEP] {
                    L::L210
                } else if iw[RADINC] < 0 || iw[RESTOR] == 1 || iw[RESTOR] == 3 {
                    L::L210
                } else {
                    w[RADFAC] = w[RDFCMX];
                    let gts = w[GTSTEP];
                    if w[FDIF] < (0.5 / w[RADFAC] - 1.0) * gts {
                        w[RADFAC] = dmax1(w[INCFAC], 0.5 * gts / (gts + w[FDIF]));
                    }
                    iw[IRC] = 4;
                    if w[STPPAR] == 0.0 {
                        L::L230
                    } else if w[DST0] >= 0.0 && (w[DST0] < 2.0 * w[DSTNRM] || w[NREDUC] < 1.2 * w[FDIF]) {
                        L::L230
                    } else {
                        iw[IRC] = 5;
                        iw[RADINC] += 1;
                        L::L200
                    }
                }
            }
            L::L200 => {
                w[FLSTGD] = w[F];
                iw[MLSTGD] = iw[MODEL];
                if iw[RESTOR] == 0 {
                    iw[RESTOR] = 2;
                }
                w[DSTSAV] = w[DSTNRM];
                iw[NFGCAL] = nfc;
                w[PLSTGD] = w[PREDUC];
                w[GTSLST] = w[GTSTEP];
                L::L230
            }
            L::L210 => {
                w[RADFAC] = 1.0;
                iw[IRC] = 3;
                L::L230
            }
            L::L220 => {
                iw[IRC] = iw[XIRC];
                if w[DSTSAV] < 0.0 {
                    iw[IRC] = 12;
                }
                L::L240
            }
            L::L230 => {
                iw[XIRC] = iw[IRC];
                L::L240
            }
            L::L240 => {
                if iw[RESTOR] == 1 && w[FLSTGD] < w[F0] {
                    iw[RESTOR] = 3;
                }
                if w[F].abs() < w[AFCTOL] {
                    iw[IRC] = 10;
                }
                if 0.5 * w[FDIF] > w[PREDUC] {
                    return;
                }
                emax = w[RFCTOL] * w[F0].abs();
                emaxs = w[SCTOL] * w[F0].abs();
                if w[PREDUC] <= emaxs && (w[DSTNRM] > w[LMAXS] || w[STPPAR] == 0.0) {
                    iw[IRC] = 11;
                }
                if w[DST0] < 0.0 {
                    L::L250
                } else {
                    let mut i2 = 0;
                    if (w[NREDUC] > 0.0 && w[NREDUC] <= emax) || (w[NREDUC] == 0.0 && w[PREDUC] == 0.0) {
                        i2 = 2;
                    }
                    if w[STPPAR] == 0.0 && w[RELDX] <= w[XCTOL] && goodx {
                        i2 += 1;
                    }
                    if i2 > 0 {
                        iw[IRC] = i2 + 6;
                    }
                    L::L250
                }
            }
            L::L250 => {
                if iw[IRC] > 5 && iw[IRC] != 12 {
                    return;
                }
                if w[STPPAR] == 0.0 {
                    return;
                }
                if w[DSTNRM] > w[LMAXS] {
                    L::L260
                } else if w[PREDUC] >= emaxs {
                    return;
                } else if w[DST0] <= 0.0 {
                    L::L270
                } else if 0.5 * w[DST0] <= w[LMAXS] {
                    return;
                } else {
                    L::L270
                }
            }
            L::L260 => {
                if 0.5 * w[DSTNRM] <= w[LMAXS] {
                    return;
                }
                let xmax = w[LMAXS] / w[DSTNRM];
                if xmax * (2.0 - xmax) * w[PREDUC] >= emaxs {
                    return;
                }
                L::L270
            }
            L::L270 => {
                if w[NREDUC] < 0.0 {
                    L::L290
                } else {
                    w[GTSLST] = w[GTSTEP];
                    w[DSTSAV] = w[DSTNRM];
                    if iw[IRC] == 12 {
                        w[DSTSAV] = -w[DSTSAV];
                    }
                    w[PLSTGD] = w[PREDUC];
                    let ii = iw[RESTOR];
                    iw[RESTOR] = 2;
                    if ii == 3 {
                        iw[RESTOR] = 0;
                    }
                    iw[IRC] = 6;
                    return;
                }
            }
            L::L280 => {
                w[GTSTEP] = w[GTSLST];
                w[DSTNRM] = w[DSTSAV].abs();
                iw[IRC] = iw[XIRC];
                if w[DSTSAV] <= 0.0 {
                    iw[IRC] = 12;
                }
                w[NREDUC] = -w[PREDUC];
                w[PREDUC] = w[PLSTGD];
                iw[RESTOR] = 3;
                L::L290
            }
            L::L290 => {
                if -w[NREDUC] <= w[SCTOL] * w[F0].abs() {
                    iw[IRC] = 11;
                }
                return;
            }
        };
    }
}

// --- trust-region step ------------------------------------------------------------------------

/// `DG7QTS(D, DIG, DIHDI, KA, L, P, STEP, V, W)`: the optimally locally constrained step
/// (Moré-Hebden) for the scaled Hessian DIHDI and scaled gradient DIG. `V` is the global V
/// (position 1 of `w`).
#[allow(clippy::too_many_arguments)]
fn dg7qts(w: &mut [f64], d: usize, dig: usize, dihdi: usize, ka: &mut i32, l: usize, p: i32, step: usize, wv: usize) {
    use vx::*;
    #[derive(Clone, Copy)]
    enum L {
        L30, L50, L60, L70, L210, L240, L250, L260, L270, L290, L310, L320, L330, L340, L370, L380, L390, L410,
    }
    const EPSFAC: f64 = 50.0;
    const KAPPA: f64 = 2.0;
    const P001: f64 = 1.0e-3;
    let big = dr7mdc(6);
    let wa = |i: i32| wv + i as usize - 1;
    let la = |i: i32| l + i as usize - 1;
    let ha = |i: i32| dihdi + i as usize - 1;
    let dgg_dmx = p + 1;
    let emax_i = dgg_dmx + 1;
    let emin_i = emax_i + 1;
    let lk0 = emin_i + 1;
    let phipin = lk0 + 1;
    let uk0 = phipin + 1;
    let dstsav = uk0 + 1;
    let diag0 = dstsav;
    let diag = diag0 + 1;
    let q0 = diag0 + p;
    let q = q0 + 1;
    let xw = q + p;
    let rad = w[RADIUS];
    let radsq = rad * rad;
    let phimax = w[PHMXFC] * rad;
    let phimin = w[PHMNFC] * rad;
    let mut psifac = big;
    let mut t1 = 2.0 * w[EPSLON] / (3.0 * (4.0 * (w[PHMNFC] + 1.0) * (KAPPA + 1.0) + KAPPA + 2.0) * rad);
    if t1 < big * dmin1(rad, 1.0) {
        psifac = t1 / rad;
    }
    let mut oldphi = 0.0;
    let eps = w[EPSLON];
    let mut irc = 0;
    let mut restrt = false;
    let mut kalim = *ka + 50;
    let (mut k, mut uk, mut lk, mut alphak) = (0i32, 0.0f64, 0.0f64, 0.0f64);
    let (mut kamin, mut phi, mut dst, mut gtsta) = (0i32, 0.0f64, 0.0f64, 0.0f64);
    let (mut t, mut delta, mut twopsi, mut sw, mut si) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut t2: f64;
    let mut lab;
    if *ka >= 0 {
        lab = L::L290;
    } else {
        k = 0;
        uk = -1.0;
        *ka = 0;
        kalim = 50;
        w[DGNORM] = dv2nrm(p, w, dig);
        w[NREDUC] = 0.0;
        w[DST0] = 0.0;
        kamin = 3;
        if w[DGNORM] == 0.0 {
            kamin = 0;
        }
        let mut j = 0;
        for i in 1..=p {
            j += i;
            let k1 = diag0 + i;
            w[wa(k1)] = w[ha(j)];
        }
        t1 = 0.0;
        let j = p * (p + 1) / 2;
        for i in 1..=j {
            let t = w[ha(i)].abs();
            if t1 < t {
                t1 = t;
            }
        }
        w[wa(dgg_dmx)] = t1;
        lab = L::L30;
    }
    let _ = k;
    loop {
        lab = match lab {
            L::L30 => {
                irc = dl7srt(1, p, w, l, dihdi);
                if irc == 0 {
                    L::L50
                } else {
                    let j = irc * (irc + 1) / 2;
                    t = w[la(j)];
                    w[la(j)] = 1.0;
                    for i in 1..=irc {
                        w[wa(i)] = 0.0;
                    }
                    w[wa(irc)] = 1.0;
                    dl7itv(irc, w, wa(1), l, wa(1));
                    t1 = dv2nrm(irc, w, wa(1));
                    lk = -t / t1 / t1;
                    w[DST0] = -lk;
                    if restrt { L::L210 } else { L::L70 }
                }
            }
            L::L50 => {
                lk = 0.0;
                t = dl7svn(p, w, l, wa(q), wa(q));
                if t >= 1.0 {
                    L::L60
                } else if w[DGNORM] >= t * t * big {
                    L::L70
                } else {
                    L::L60
                }
            }
            L::L60 => {
                dl7ivm(p, w, wa(q), l, dig);
                gtsta = dd7tpr(p, w, wa(q), wa(q));
                w[NREDUC] = 0.5 * gtsta;
                dl7itv(p, w, wa(q), l, wa(q));
                dst = dv2nrm(p, w, wa(q));
                w[DST0] = dst;
                phi = dst - rad;
                if phi <= phimax {
                    L::L260
                } else if restrt {
                    L::L210
                } else {
                    L::L70
                }
            }
            L::L70 => {
                let mut kk = 0;
                for i in 1..=p {
                    let mut wi = 0.0;
                    if i != 1 {
                        let im1 = i - 1;
                        for j in 1..=im1 {
                            kk += 1;
                            let tt = w[ha(kk)].abs();
                            wi += tt;
                            w[wa(j)] += tt;
                        }
                    }
                    w[wa(i)] = wi;
                    kk += 1;
                }
                // smallest eigenvalue bound
                let mut kidx = 1;
                t1 = w[wa(diag)] - w[wa(1)];
                if p > 1 {
                    for i in 2..=p {
                        let j = diag0 + i;
                        let tt = w[wa(j)] - w[wa(i)];
                        if tt < t1 {
                            t1 = tt;
                            kidx = i;
                        }
                    }
                }
                let mut sk = w[wa(kidx)];
                let mut j = diag0 + kidx;
                let mut akk = w[wa(j)];
                let mut k1 = kidx * (kidx - 1) / 2 + 1;
                let mut inc = 1;
                t = 0.0;
                for i in 1..=p {
                    if i == kidx {
                        inc = i; // 130
                    } else {
                        let aki = w[ha(k1)].abs();
                        let si_ = w[wa(i)];
                        j = diag0 + i;
                        t1 = 0.5 * (akk - w[wa(j)] + si_ - aki);
                        t1 += (t1 * t1 + sk * aki).sqrt();
                        if t < t1 {
                            t = t1;
                        }
                        if i >= kidx {
                            inc = i; // 130
                        }
                    }
                    k1 += inc; // 140
                }
                w[wa(emin_i)] = akk - t;
                uk = w[DGNORM] / rad - w[wa(emin_i)];
                if w[DGNORM] == 0.0 {
                    uk = uk + P001 + P001 * uk;
                }
                if uk <= 0.0 {
                    uk = P001;
                }
                // largest eigenvalue bound
                kidx = 1;
                t1 = w[wa(diag)] + w[wa(1)];
                if p > 1 {
                    for i in 2..=p {
                        let j = diag0 + i;
                        let tt = w[wa(j)] + w[wa(i)];
                        if tt > t1 {
                            t1 = tt;
                            kidx = i;
                        }
                    }
                }
                sk = w[wa(kidx)];
                j = diag0 + kidx;
                akk = w[wa(j)];
                k1 = kidx * (kidx - 1) / 2 + 1;
                inc = 1;
                t = 0.0;
                for i in 1..=p {
                    if i == kidx {
                        inc = i; // 180
                    } else {
                        let aki = w[ha(k1)].abs();
                        let si_ = w[wa(i)];
                        j = diag0 + i;
                        t1 = 0.5 * (w[wa(j)] + si_ - aki - akk);
                        t1 += (t1 * t1 + sk * aki).sqrt();
                        if t < t1 {
                            t = t1;
                        }
                        if i >= kidx {
                            inc = i; // 180
                        }
                    }
                    k1 += inc; // 190
                }
                w[wa(emax_i)] = akk + t;
                lk = dmax1(lk, w[DGNORM] / rad - w[wa(emax_i)]);
                alphak = w[STPPAR].abs() * w[RAD0] / rad;
                alphak = dmin1(uk, dmax1(alphak, lk));
                if irc != 0 {
                    L::L210
                } else {
                    dl7ivm(p, w, wa(1), l, wa(q));
                    t = dv2nrm(p, w, wa(1));
                    w[wa(phipin)] = rad / t / t;
                    lk = dmax1(lk, phi * w[wa(phipin)]);
                    L::L210
                }
            }
            L::L210 => {
                *ka += 1;
                if -w[DST0] >= alphak || alphak < lk || alphak >= uk {
                    alphak = uk * dmax1(P001, (lk / uk).sqrt());
                }
                if alphak <= 0.0 {
                    alphak = 0.5 * uk;
                }
                if alphak <= 0.0 {
                    alphak = uk;
                }
                let mut kk = 0;
                for i in 1..=p {
                    kk += i;
                    let j = diag0 + i;
                    w[ha(kk)] = w[wa(j)] + alphak;
                }
                irc = dl7srt(1, p, w, l, dihdi);
                if irc == 0 {
                    L::L240
                } else {
                    let j = (irc * (irc + 1)) / 2;
                    t = w[la(j)];
                    w[la(j)] = 1.0;
                    for i in 1..=irc {
                        w[wa(i)] = 0.0;
                    }
                    w[wa(irc)] = 1.0;
                    dl7itv(irc, w, wa(1), l, wa(1));
                    t1 = dv2nrm(irc, w, wa(1));
                    lk = alphak - t / t1 / t1;
                    w[DST0] = -lk;
                    if uk < lk {
                        uk = lk;
                    }
                    if alphak < lk {
                        L::L210
                    } else {
                        t = P001 * alphak;
                        if t <= 0.0 {
                            t = P001;
                        }
                        lk = alphak + t;
                        if uk <= lk {
                            uk = lk + t;
                        }
                        L::L210
                    }
                }
            }
            L::L240 => {
                dl7ivm(p, w, wa(q), l, dig);
                gtsta = dd7tpr(p, w, wa(q), wa(q));
                dl7itv(p, w, wa(q), l, wa(q));
                dst = dv2nrm(p, w, wa(q));
                phi = dst - rad;
                if phi <= phimax && phi >= phimin {
                    L::L270
                } else if phi == oldphi {
                    L::L270
                } else {
                    oldphi = phi;
                    if phi < 0.0 { L::L330 } else { L::L250 }
                }
            }
            L::L250 => {
                if *ka >= kalim {
                    L::L270
                } else {
                    if phi < 0.0 {
                        uk = dmin1(uk, alphak);
                    }
                    if kamin == 0 {
                        L::L210
                    } else {
                        dl7ivm(p, w, wa(1), l, wa(q));
                        t1 = dv2nrm(p, w, wa(1));
                        alphak += (phi / t1) * (dst / t1) * (dst / rad);
                        lk = dmax1(lk, alphak);
                        alphak = lk;
                        L::L210
                    }
                }
            }
            L::L260 => {
                alphak = 0.0;
                L::L270
            }
            L::L270 => {
                for i in 1..=p {
                    let j = q0 + i;
                    w[step + i as usize - 1] = -w[wa(j)] / w[d + i as usize - 1];
                }
                w[GTSTEP] = -gtsta;
                w[PREDUC] = 0.5 * (alphak.abs() * dst * dst + gtsta);
                L::L410
            }
            L::L290 => {
                if w[DST0] <= 0.0 || w[DST0] - rad > phimax {
                    L::L310
                } else {
                    restrt = true;
                    *ka += 1;
                    let mut kk = 0;
                    for i in 1..=p {
                        kk += i;
                        let j = diag0 + i;
                        w[ha(kk)] = w[wa(j)];
                    }
                    uk = -1.0;
                    L::L30
                }
            }
            L::L310 => {
                kamin = *ka + 3;
                if w[DGNORM] == 0.0 {
                    kamin = 0;
                }
                if *ka == 0 {
                    L::L50
                } else {
                    dst = w[wa(dstsav)];
                    alphak = w[STPPAR].abs();
                    phi = dst - rad;
                    t = w[DGNORM] / rad;
                    uk = t - w[wa(emin_i)];
                    if w[DGNORM] == 0.0 {
                        uk = uk + P001 + P001 * uk;
                    }
                    if uk <= 0.0 {
                        uk = P001;
                    }
                    if rad > w[RAD0] {
                        L::L320
                    } else {
                        lk = 0.0;
                        if alphak > 0.0 {
                            lk = w[wa(lk0)];
                        }
                        lk = dmax1(lk, t - w[wa(emax_i)]);
                        if w[DST0] > 0.0 {
                            lk = dmax1(lk, (w[DST0] - rad) * w[wa(phipin)]);
                        }
                        L::L250
                    }
                }
            }
            L::L320 => {
                if alphak > 0.0 {
                    uk = dmin1(uk, w[wa(uk0)]);
                }
                lk = dmax1(dmax1(0.0, -w[DST0]), t - w[wa(emax_i)]);
                if w[DST0] > 0.0 {
                    lk = dmax1(lk, (w[DST0] - rad) * w[wa(phipin)]);
                }
                L::L250
            }
            L::L330 => {
                delta = alphak + dmin1(0.0, w[DST0]);
                twopsi = alphak * dst * dst + gtsta;
                if *ka >= kamin || psifac >= big {
                    L::L340
                } else if delta >= psifac * twopsi {
                    L::L370
                } else {
                    L::L340
                }
            }
            L::L340 => {
                t = dl7svn(p, w, l, wa(xw), wa(1));
                for i in 1..=p {
                    w[wa(i)] *= t;
                }
                dl7itv(p, w, wa(1), l, wa(1));
                t2 = 1.0 / dv2nrm(p, w, wa(1));
                for i in 1..=p {
                    w[wa(i)] *= t2;
                }
                t *= t2;
                sw = dd7tpr(p, w, wa(q), wa(1));
                t1 = (rad + dst) * (rad - dst);
                let mut root = (sw * sw + t1).sqrt();
                if sw < 0.0 {
                    root = -root;
                }
                si = t1 / (sw + root);
                if (t2 * si) * (t2 * si) <= eps * (dst * dst + alphak * radsq) {
                    L::L380
                } else {
                    if w[DST0] <= 0.0 {
                        w[DST0] = dmin1(w[DST0], t2 * t2 - alphak);
                    }
                    lk = dmax1(lk, -w[DST0]);
                    L::L370
                }
            }
            L::L370 => {
                let dgxfac = EPSFAC * dr7mdc(3);
                if delta > dgxfac * w[wa(dgg_dmx)] { L::L250 } else { L::L270 }
            }
            L::L380 => {
                alphak = -alphak;
                w[PREDUC] = 0.5 * twopsi;
                t1 = 0.0;
                let tt = si * (alphak * sw - 0.5 * si * (alphak + t * dd7tpr(p, w, wa(xw), wa(1))));
                if tt >= eps * twopsi / 6.0 {
                    w[PREDUC] += tt;
                    dst = rad;
                    t1 = -si;
                }
                L::L390
            }
            L::L390 => {
                for i in 1..=p {
                    let j = q0 + i;
                    w[wa(j)] = t1 * w[wa(i)] - w[wa(j)];
                    w[step + i as usize - 1] = w[wa(j)] / w[d + i as usize - 1];
                }
                w[GTSTEP] = dd7tpr(p, w, dig, wa(q));
                L::L410
            }
            L::L410 => {
                w[DSTNRM] = dst;
                w[STPPAR] = alphak;
                w[wa(lk0)] = lk;
                w[wa(uk0)] = uk;
                w[RAD0] = rad;
                w[wa(dstsav)] = dst;
                let mut j = 0;
                for i in 1..=p {
                    j += i;
                    let kk = diag0 + i;
                    w[ha(j)] = w[wa(kk)];
                }
                return;
            }
        };
    }
}

/// `DS7BQN(B, D, DST, IPIV, IPIV1, IPIV2, KB, L, LV, NS, P, P1, STEP, TD, TG, V, W, X, X0)`:
/// bounded modified Newton step. `b` is the position of `B(1,1)`; the IPIV arrays are in `iw`.
#[allow(clippy::too_many_arguments)]
fn ds7bqn(w: &mut [f64], iw: &mut [i32], b: usize, d: usize, dst: usize, ipiv: usize, ipiv1: usize, ipiv2: usize, kb: &mut i32, l: usize, ns: &mut i32, p: i32, p1: &mut i32, step: usize, td: usize, tg: usize, wv: usize, x: usize, x0: usize) {
    use vx::*;
    let b1 = |j: i32| b + 2 * (j as usize - 1);
    let b2 = |j: i32| b + 2 * (j as usize - 1) + 1;
    let at = |base: usize, i: i32| base + i as usize - 1;
    let dstmax = 1.0001 * (1.0 + w[PHMXFC]) * w[RADIUS];
    let dstmin = (1.0 + w[PHMNFC]) * w[RADIUS];
    let mut dst1 = 0.0;
    let meps2 = 2.0 * dr7mdc(3);
    let p0 = *p1;
    *ns = 0;
    for i in 1..=p {
        iw[at(ipiv1, i)] = i;
        iw[at(ipiv2, i)] = i;
    }
    for i in 1..=*p1 {
        w[at(wv, i)] = -w[at(step, i)] * w[at(td, i)];
    }
    let alpha = w[STPPAR].abs();
    w[PREDUC] = 0.0;
    let mut gts = -w[GTSTEP];
    if *kb < 0 {
        dv7scp(p, w, dst, 0.0);
    }
    *kb = 1;
    loop {
        // 30
        let mut t = 1.0;
        let mut k = 0i32;
        for i in 1..=*p1 {
            let j = iw[at(ipiv, i)];
            let dx = w[at(wv, i)] / w[at(d, j)];
            let xi = w[at(x, j)] - dx;
            let ti;
            if xi < w[b1(j)] {
                ti = (w[at(x, j)] - w[b1(j)]) / dx;
                k = -i;
            } else if xi <= w[b2(j)] {
                continue;
            } else {
                ti = (w[at(x, j)] - w[b2(j)]) / dx;
                k = i;
            }
            // 50
            if t > ti {
                t = ti;
            }
        }
        if p > *p1 {
            dv7cpy(p - *p1, w, at(step, *p1 + 1), at(dst, *p1 + 1));
        }
        dv2axy(*p1, w, step, -t, wv, dst);
        let dst0 = dst1;
        dst1 = dv2nrm(p, w, step);
        // check for oversize step
        if dst1 > dstmax {
            if *p1 >= p0 {
                k = 0; // 70
            } else {
                if dst0 < dstmin {
                    *kb = 0;
                }
                break; // GO TO 110
            }
        }
        // 80
        w[DSTNRM] = dst1;
        dv7cpy(*p1, w, dst, step);
        let t1 = 1.0 - t;
        for i in 1..=*p1 {
            w[at(tg, i)] *= t1;
        }
        if alpha > 0.0 {
            dv2axy(*p1, w, tg, t * alpha, wv, tg);
        }
        w[PREDUC] += t * ((1.0 - 0.5 * t) * gts + 0.5 * alpha * t * dd7tpr(*p1, w, wv, wv));
        if k == 0 {
            break;
        }
        // permute L, etc. if necessary
        let p1m1 = *p1 - 1;
        let j = k.abs();
        if j != *p1 {
            *ns += 1;
            iw[at(ipiv2, *p1)] = j;
            dq7rsh(j, *p1, w, l, wv);
            i7shft(*p1, j, iw, ipiv);
            i7shft(*p1, j, iw, ipiv1);
            dv7shf(*p1, j, w, tg);
            dv7shf(*p1, j, w, dst);
        }
        // 100
        if k < 0 {
            iw[at(ipiv, *p1)] = -iw[at(ipiv, *p1)];
        }
        *p1 = p1m1;
        if *p1 <= 0 {
            break;
        }
        dl7ivm(*p1, w, wv, l, tg);
        gts = dd7tpr(*p1, w, wv, wv);
        dl7itv(*p1, w, wv, l, wv);
    }
    // 110: unscale step
    for i in 1..=p {
        let j = iw[at(ipiv, i)].abs();
        w[at(step, j)] = w[at(dst, i)] / w[at(d, j)];
    }
    // fudge step to force components to their bounds
    if *p1 < p0 {
        let k = *p1 + 1;
        for i in k..=p0 {
            let mut j = iw[at(ipiv, i)];
            let mut t = meps2;
            if j <= 0 {
                t = -t;
                j = -j;
                iw[at(ipiv, i)] = j;
            }
            t *= dmax1(w[at(x, j)].abs(), w[at(x0, j)].abs());
            w[at(step, j)] += t;
        }
    }
    // 150
    dv2axy(p, w, x, 1.0, step, x0);
    if *ns > 0 {
        dv7ipr(p0, iw, ipiv1, w, td);
    }
}

/// `DG7QSB(B, D, DIHDI, G, IPIV, IPIV1, IPIV2, KA, L, LV, P, P0, PC, STEP, TD, TG, V, W, X, X0)`:
/// heuristic bounded Newton step.
#[allow(clippy::too_many_arguments)]
fn dg7qsb(w: &mut [f64], iw: &mut [i32], b: usize, d: usize, dihdi: usize, g: usize, ipiv: usize, ipiv1: usize, ipiv2: usize, ka: &mut i32, l: usize, p: i32, p0: &mut i32, pc: i32, step: usize, td: usize, tg: usize, wv: usize, x: usize, x0: usize) {
    use vx::*;
    let mut p1 = pc;
    let (mut nred, mut ds0) = (0.0, 0.0);
    if *ka >= 0 {
        nred = w[NREDUC];
        ds0 = w[DST0];
    } else {
        *p0 = 0;
        *ka = -1;
    }
    // 20
    let mut kinit = -1;
    if *p0 == p1 {
        kinit = *ka;
    }
    dv7cpy(p, w, x, x0);
    let mut pred = 0.0;
    let rad = w[RADIUS];
    let mut kb = -1;
    w[DSTNRM] = 0.0;
    if p1 <= 0 {
        nred = 0.0;
        ds0 = 0.0;
        dv7scp(p, w, step, 0.0);
    } else {
        dv7cpy(p, w, td, d);
        dv7ipr(p, iw, ipiv, w, td);
        dv7vmp(p, w, tg, g, d, -1);
        dv7ipr(p, iw, ipiv, w, tg);
        loop {
            // 40
            let mut k = kinit;
            kinit = -1;
            w[RADIUS] = rad - w[DSTNRM];
            dg7qts(w, td, tg, dihdi, &mut k, l, p1, step, wv);
            *p0 = p1;
            if *ka < 0 {
                nred = w[NREDUC];
                ds0 = w[DST0];
            }
            // 50
            *ka = k;
            w[RADIUS] = rad;
            let p10 = p1;
            let mut ns = 0;
            ds7bqn(w, iw, b, d, step + p as usize, ipiv, ipiv1, ipiv2, &mut kb, l, &mut ns, p, &mut p1, step, td, tg, wv, x, x0);
            if ns > 0 {
                ds7ipr(p10, iw, ipiv1, w, dihdi);
            }
            pred += w[PREDUC];
            if ns != 0 {
                *p0 = 0;
            }
            if kb > 0 {
                break;
            }
        }
    }
    // 60
    w[DST0] = ds0;
    w[NREDUC] = nred;
    w[PREDUC] = pred;
    w[GTSTEP] = dd7tpr(p, w, g, step);
}

// --- routines of DRMNGB (no Hessian: BFGS updates of a Cholesky factor) -------------------------

/// `DL7TVM(N, X, L, Y)`: `X = L**T * Y`.
///
/// DRMNGB calls it with X and Y the same array (`CALL DL7TVM(N, V(X01), V(L), V(X01))`), which
/// Fortran does not allow, and the compiled code in R's `stats.so` does not do what the source
/// says then: gfortran stores `X(I) = 0` before it loads `Y(I)` (`dl7tvm_`: `movq $0x0, X(I)`,
/// then `movsd Y(I)`), so with aliased arguments every `YI` reads a zero and the result is all
/// zeros. The port keeps that order, which is the same as the source's when X and Y differ.
fn dl7tvm(n: i32, w: &mut [f64], x: usize, l: usize, y: usize) {
    let mut i0 = 0;
    for i in 1..=n {
        w[x + i as usize - 1] = 0.0;
        let yi = w[y + i as usize - 1];
        for j in 1..=i {
            let ij = i0 + j;
            w[x + j as usize - 1] += yi * w[l + ij as usize - 1];
        }
        i0 += i;
    }
}

/// `DL7VML(N, X, L, Y)`: `X = L * Y` (X and Y may be the same array).
fn dl7vml(n: i32, w: &mut [f64], x: usize, l: usize, y: usize) {
    let np1 = n + 1;
    let mut i0 = n * (n + 1) / 2;
    for ii in 1..=n {
        let i = np1 - ii;
        i0 -= i;
        let mut t = 0.0;
        for j in 1..=i {
            let ij = i0 + j;
            t += w[l + ij as usize - 1] * w[y + j as usize - 1];
        }
        w[x + i as usize - 1] = t;
    }
}

/// `DW7ZBF(L, N, S, W, Y, Z)`: the W and Z of the BFGS update of the Cholesky factor L.
fn dw7zbf(w: &mut [f64], l: usize, n: i32, s: usize, wv: usize, y: usize, z: usize) {
    const EPS: f64 = 0.1;
    dl7tvm(n, w, wv, l, s);
    let shs = dd7tpr(n, w, wv, wv);
    let ys = dd7tpr(n, w, y, s);
    let (cy, cs);
    if ys >= EPS * shs {
        cy = 1.0 / (ys.sqrt() * shs.sqrt());
        cs = 1.0 / shs;
    } else {
        let theta = (1.0 - EPS) * shs / (shs - ys);
        let epsrt = EPS.sqrt();
        cy = theta / (shs * epsrt);
        cs = (1.0 + (theta - 1.0) / epsrt) / shs;
    }
    dl7ivm(n, w, z, l, y);
    for i in 0..n as usize {
        w[z + i] = cy * w[z + i] - cs * w[wv + i];
    }
}

/// `DL7UPD(BETA, GAMMA, L, LAMBDA, LPLUS, N, W, Z)`: `LPLUS = L * (I + Z*W**T)` (LPLUS may be L).
#[allow(clippy::too_many_arguments)]
fn dl7upd(w: &mut [f64], beta: usize, gamma: usize, l: usize, lambda: usize, lplus: usize, n: i32, wv: usize, z: usize) {
    let at = |base: usize, i: i32| base + i as usize - 1;
    let mut nu = 1.0;
    let mut eta = 0.0;
    if n > 1 {
        let nm1 = n - 1;
        let mut s = 0.0;
        for i in 1..=nm1 {
            let j = n - i;
            let wj1 = w[at(wv, j + 1)];
            s += wj1 * wj1;
            w[at(lambda, j)] = s;
        }
        for j in 1..=nm1 {
            let wj = w[at(wv, j)];
            let a = nu * w[at(z, j)] - eta * wj;
            let theta = 1.0 + a * wj;
            let s = a * w[at(lambda, j)];
            let mut lj = (theta * theta + a * s).sqrt();
            if theta > 0.0 {
                lj = -lj;
            }
            w[at(lambda, j)] = lj;
            let b = theta * wj + s;
            w[at(gamma, j)] = b * nu / lj;
            w[at(beta, j)] = (a - b * eta) / lj;
            nu = -nu / lj;
            eta = -(eta + (a * a) / (theta - lj)) / lj;
        }
    }
    w[at(lambda, n)] = 1.0 + (nu * w[at(z, n)] - eta * w[at(wv, n)]) * w[at(wv, n)];
    let np1 = n + 1;
    let mut jj = n * (n + 1) / 2;
    for k in 1..=n {
        let j = np1 - k;
        let lj = w[at(lambda, j)];
        let ljj = w[at(l, jj)];
        w[at(lplus, jj)] = lj * ljj;
        let wj = w[at(wv, j)];
        w[at(wv, j)] = ljj * wj;
        let zj = w[at(z, j)];
        w[at(z, j)] = ljj * zj;
        if k != 1 {
            let bj = w[at(beta, j)];
            let gj = w[at(gamma, j)];
            let mut ij = jj + j;
            let jp1 = j + 1;
            for i in jp1..=n {
                let lij = w[at(l, ij)];
                w[at(lplus, ij)] = lj * lij + bj * w[at(wv, i)] + gj * w[at(z, i)];
                w[at(wv, i)] += lij * wj;
                w[at(z, i)] += lij * zj;
                ij += i;
            }
        }
        jj -= j;
    }
}

/// `DD7DOG(DIG, LV, N, NWTSTP, STEP, V)`: the double-dogleg step.
fn dd7dog(w: &mut [f64], dig: usize, n: i32, nwtstp: usize, step: usize) {
    use vx::*;
    const GTHG: usize = 44;
    const GRDFAC: usize = 45;
    const NWTFAC: usize = 46;
    let nwtnrm = w[DST0];
    let mut rlambd = 1.0;
    if nwtnrm > 0.0 {
        rlambd = w[RADIUS] / nwtnrm;
    }
    let gnorm = w[DGNORM];
    let ghinvg = 2.0 * w[NREDUC];
    w[GRDFAC] = 0.0;
    w[NWTFAC] = 0.0;
    if rlambd >= 1.0 {
        w[STPPAR] = 0.0;
        w[DSTNRM] = nwtnrm;
        w[GTSTEP] = -ghinvg;
        w[PREDUC] = w[NREDUC];
        w[NWTFAC] = -1.0;
        for i in 0..n as usize {
            w[step + i] = -w[nwtstp + i];
        }
        return;
    }
    w[DSTNRM] = w[RADIUS];
    let cf = gnorm / w[GTHG];
    let cfact = cf * cf;
    let cnorm = gnorm * cfact;
    let relax = 1.0 - w[BIAS] * (1.0 - gnorm * cnorm / ghinvg);
    if rlambd >= relax {
        w[STPPAR] = 1.0 - (rlambd - relax) / (1.0 - relax);
        let t = -rlambd;
        w[GTSTEP] = t * ghinvg;
        w[PREDUC] = rlambd * (1.0 - 0.5 * rlambd) * ghinvg;
        w[NWTFAC] = t;
        for i in 0..n as usize {
            w[step + i] = t * w[nwtstp + i];
        }
        return;
    }
    if cnorm >= w[RADIUS] {
        let t = -w[RADIUS] / gnorm;
        w[GRDFAC] = t;
        w[STPPAR] = 1.0 + cnorm / w[RADIUS];
        w[GTSTEP] = -w[RADIUS] * gnorm;
        let q = w[GTHG] / gnorm;
        w[PREDUC] = w[RADIUS] * (gnorm - 0.5 * w[RADIUS] * (q * q));
        for i in 0..n as usize {
            w[step + i] = t * w[dig + i];
        }
        return;
    }
    let ctrnwt = cfact * relax * ghinvg / gnorm;
    let mut t1 = ctrnwt - gnorm * (cfact * cfact);
    let mut t2 = w[RADIUS] * (w[RADIUS] / gnorm) - gnorm * (cfact * cfact);
    let mut t = relax * nwtnrm;
    let femnsq = (t / gnorm) * t - ctrnwt - t1;
    t = t2 / (t1 + (t1 * t1 + femnsq * t2).sqrt());
    t1 = (t - 1.0) * cfact;
    w[GRDFAC] = t1;
    t2 = -t * relax;
    w[NWTFAC] = t2;
    w[STPPAR] = 2.0 - t;
    w[GTSTEP] = t1 * (gnorm * gnorm) + t2 * ghinvg;
    let q = w[GTHG] * t1;
    w[PREDUC] = -t1 * gnorm * ((t2 + 1.0) * gnorm) - t2 * (1.0 + 0.5 * t2) * ghinvg - 0.5 * (q * q);
    for i in 0..n as usize {
        w[step + i] = t1 * w[dig + i] + t2 * w[nwtstp + i];
    }
}

/// `DD7DGB(B, D, DIG, DST, G, IPIV, KA, L, LV, P, PC, NWTST, STEP, TD, TG, V, W, X0)`: the
/// bounded double-dogleg step.
#[allow(clippy::too_many_arguments)]
fn dd7dgb(w: &mut [f64], iw: &mut [i32], b: usize, d: usize, dig: usize, dst: usize, g: usize, ipiv: usize, ka: &mut i32, l: usize, p: i32, pc: i32, nwtst: usize, step: usize, td: usize, tg: usize, wv: usize, x0: usize) {
    use vx::*;
    const GTHG: usize = 44;
    const GRDFAC: usize = 45;
    const NWTFAC: usize = 46;
    let at = |base: usize, i: i32| base + i as usize - 1;
    let b1 = |j: i32| b + 2 * (j as usize - 1);
    let b2 = |j: i32| b + 2 * (j as usize - 1) + 1;
    let meps2 = 2.0 * dr7mdc(3);
    let gnorm0 = w[DGNORM];
    w[DSTNRM] = 0.0;
    let (mut dnwtst, mut nred) = (0.0, 0.0);
    if *ka >= 0 {
        dnwtst = w[DST0];
        nred = w[NREDUC];
    }
    let mut pred = 0.0;
    w[STPPAR] = 0.0;
    let rad = w[RADIUS];
    let mut p1 = pc;
    if pc <= 0 {
        dnwtst = 0.0;
        dv7scp(p, w, step, 0.0);
    } else {
        dv7cpy(p, w, td, d);
        dv7ipr(p, iw, ipiv, w, td);
        dv7scp(pc, w, dst, 0.0);
        dv7cpy(p, w, tg, g);
        dv7ipr(p, iw, ipiv, w, tg);
        loop {
            // 30
            dl7ivm(p1, w, nwtst, l, tg);
            let ghinvg = dd7tpr(p1, w, nwtst, nwtst);
            w[NREDUC] = 0.5 * ghinvg;
            dl7itv(p1, w, nwtst, l, nwtst);
            dv7vmp(p1, w, step, nwtst, td, 1);
            w[DST0] = dv2nrm(pc, w, step);
            if *ka < 0 {
                *ka = 0;
                dnwtst = w[DST0];
                nred = w[NREDUC];
            }
            // 40
            w[RADIUS] = rad - w[DSTNRM];
            if w[RADIUS] <= 0.0 {
                break;
            }
            dv7vmp(p1, w, dig, tg, td, -1);
            let gnorm = dv2nrm(p1, w, dig);
            if gnorm <= 0.0 {
                break;
            }
            w[DGNORM] = gnorm;
            dv7vmp(p1, w, dig, dig, td, -1);
            dl7tvm(p1, w, wv, l, dig);
            w[GTHG] = dv2nrm(p1, w, wv);
            *ka += 1;
            dd7dog(w, dig, p1, nwtst, step);
            let mut t = 1.0;
            let mut k = 0i32;
            for i in 1..=p1 {
                let mut j = iw[at(ipiv, i)];
                let x0i = w[at(x0, j)] + w[at(dst, i)] / w[at(td, i)];
                let xi = x0i + w[at(step, i)];
                let ti;
                if xi < w[b1(j)] {
                    ti = (w[b1(j)] - x0i) / w[at(step, i)];
                    j = -i;
                } else if xi <= w[b2(j)] {
                    continue;
                } else {
                    ti = (w[b2(j)] - x0i) / w[at(step, i)];
                    j = i;
                }
                // 60
                if t > ti {
                    k = j;
                    t = ti;
                }
            }
            dv7vmp(p1, w, step, step, td, 1);
            dv2axy(p1, w, dst, t, step, dst);
            w[DSTNRM] = dv2nrm(pc, w, dst);
            let t1 = t * w[GRDFAC];
            let t2 = t * w[NWTFAC];
            let q = w[GTHG] * t1;
            pred = pred - t1 * gnorm * ((t2 + 1.0) * gnorm) - t2 * (1.0 + 0.5 * t2) * ghinvg - 0.5 * (q * q);
            if k == 0 {
                break;
            }
            dl7vml(p1, w, wv, l, wv);
            let t2 = 1.0 - t2;
            for i in 1..=p1 {
                w[at(tg, i)] = t2 * w[at(tg, i)] - t1 * w[at(wv, i)];
            }
            let p1m1 = p1 - 1;
            let j = k.abs();
            if j != p1 {
                dq7rsh(j, p1, w, l, wv);
                i7shft(p1, j, iw, ipiv);
                dv7shf(p1, j, w, tg);
                dv7shf(p1, j, w, td);
                dv7shf(p1, j, w, dst);
            }
            // 90
            if k < 0 {
                iw[at(ipiv, p1)] = -iw[at(ipiv, p1)];
            }
            p1 = p1m1;
            if p1 <= 0 {
                break;
            }
        }
        // 100
        dv7scp(p, w, step, 0.0);
        for i in 1..=pc {
            let j = iw[at(ipiv, i)].abs();
            w[at(step, j)] = w[at(dst, i)] / w[at(td, i)];
        }
        if p1 < pc {
            dv2axy(p, w, td, 1.0, step, x0);
            let k = p1 + 1;
            for i in k..=pc {
                let mut j = iw[at(ipiv, i)];
                let mut t = meps2;
                if j <= 0 {
                    t = -t;
                    j = -j;
                    iw[at(ipiv, i)] = j;
                }
                t *= dmax1(w[at(td, j)].abs(), w[at(x0, j)].abs());
                w[at(step, j)] += t;
            }
        }
    }
    // 140
    w[DGNORM] = gnorm0;
    w[NREDUC] = nred;
    w[PREDUC] = pred;
    w[RADIUS] = rad;
    w[DST0] = dnwtst;
    w[GTSTEP] = dd7tpr(p, w, step, g);
}

/// `DRMNGB(B, D, FX, G, IV, LIV, LV, N, V, X)`: one reverse-communication call.
fn drmngb(w: &mut [f64], iw: &mut [i32], lay: &Layout, fx: f64) {
    use vx::*;
    #[derive(Clone, Copy)]
    enum L {
        L10, L20, L30, L70, L80, L140, L150, L160, L170, L180, L190, L200, L210, L220, L230, L240, L250, L260, L330, L340, L360, L370, L390, L400, L410, L420, L430, L440,
    }
    const STLSTG: usize = 41;
    const NWTSTP: usize = 34;
    const INITH: usize = 25;
    let n = lay.n as i32;
    let (liv, lv) = (lay.liv, lay.lv);
    let b1 = |i: i32| lay.b + 2 * (i as usize - 1);
    let b2 = |i: i32| lay.b + 2 * (i as usize - 1) + 1;
    let dat = |i: i32| lay.d + i as usize - 1;
    let gat = |i: i32| lay.g + i as usize - 1;
    let xat = |i: i32| lay.x + i as usize - 1;
    let mut rstrst = 3;
    let mut lab = match iw[1] {
        1 => L::L70,
        2 => L::L80,
        _ => {
            if iw[1] == 0 {
                divset(2, iw, liv, lv, w);
            }
            if (12..=13).contains(&iw[1]) {
                iw[VNEED] += n * (n + 19) / 2;
                iw[IVNEED] += n;
            }
            L::L10
        }
    };
    loop {
        lab = match lab {
            L::L10 => {
                dparck(2, lay.d, iw, liv, lv, n, w);
                let i = iw[1] - 2;
                if i > 12 {
                    return;
                }
                match i {
                    1..=6 => L::L250,
                    7 | 9 => L::L190,
                    8 => L::L150,
                    12 => L::L30,
                    _ => L::L20,
                }
            }
            L::L20 => {
                let l = iw[LMAT];
                iw[X0] = l + n * (n + 1) / 2;
                iw[STEP] = iw[X0] + 2 * n;
                iw[STLSTG] = iw[STEP] + 2 * n;
                iw[NWTSTP] = iw[STLSTG] + n;
                iw[DG] = iw[NWTSTP] + 2 * n;
                iw[NEXTV] = iw[DG] + 2 * n;
                iw[NEXTIV] = iw[PERM] + n;
                if iw[1] != 13 {
                    L::L30
                } else {
                    iw[1] = 14;
                    return;
                }
            }
            L::L30 => {
                iw[NITER] = 0;
                iw[NFCALL] = 1;
                iw[NGCALL] = 1;
                iw[NFGCAL] = 1;
                iw[MODE] = -1;
                iw[MODEL] = 1;
                iw[STGLIM] = 1;
                iw[TOOBIG] = 0;
                iw[CNVCOD] = 0;
                iw[RADINC] = 0;
                iw[NC] = n;
                w[RAD0] = 0.0;
                let mut ipi = iw[PERM] as usize;
                let mut bad = false;
                for i in 1..=n {
                    iw[ipi] = i;
                    ipi += 1;
                    if w[b1(i)] > w[b2(i)] {
                        bad = true;
                        break;
                    }
                }
                if bad {
                    L::L410
                } else {
                    if w[DINIT] >= 0.0 {
                        let c = w[DINIT];
                        dv7scp(n, w, lay.d, c);
                    }
                    if iw[INITH] == 1 {
                        let l = iw[LMAT] as usize;
                        dv7scp(n * (n + 1) / 2, w, l, 0.0);
                        let mut k = l - 1;
                        for i in 1..=n {
                            k += i as usize;
                            let mut t = w[dat(i)];
                            if t <= 0.0 {
                                t = 1.0;
                            }
                            w[k] = t;
                        }
                    }
                    iw[1] = 1;
                    L::L440
                }
            }
            L::L70 => {
                w[F] = fx;
                if iw[MODE] >= 0 {
                    L::L250
                } else {
                    w[F0] = fx;
                    iw[1] = 2;
                    if iw[TOOBIG] == 0 {
                        return;
                    }
                    iw[1] = 63;
                    L::L430
                }
            }
            L::L80 => {
                if iw[TOOBIG] != 0 {
                    iw[1] = 65;
                    L::L430
                } else {
                    // 90
                    let ipi = iw[PERM] as usize;
                    let mut ipn = ipi + n as usize;
                    let mut n1 = n;
                    let np1 = n + 1;
                    let l = iw[LMAT] as usize;
                    let w1 = (iw[NWTSTP] + n) as usize;
                    let k = n - iw[NC];
                    for i in 1..=n {
                        ipn -= 1;
                        let j = iw[ipn];
                        let free = !(w[b1(j)] >= w[b2(j)]) && {
                            let xi = w[xat(j)];
                            let gi = w[gat(j)];
                            !(xi <= w[b1(j)] && gi > 0.0) && !(xi >= w[b2(j)] && gi < 0.0)
                        };
                        if free {
                            if i <= k {
                                iw[CNVCOD] = 0;
                            }
                            continue;
                        }
                        // 100
                        let i1 = np1 - i;
                        if i1 < n1 {
                            i7shft(n1, i1, iw, ipi);
                            dq7rsh(i1, n1, w, l, w1);
                        }
                        n1 -= 1;
                    }
                    iw[NC] = n1;
                    w[DGNORM] = 0.0;
                    if n1 > 0 {
                        let dg1 = iw[DG] as usize;
                        dv7vmp(n, w, dg1, lay.g, lay.d, -1);
                        dv7ipr(n, iw, ipi, w, dg1);
                        w[DGNORM] = dv2nrm(n1, w, dg1);
                    }
                    if iw[CNVCOD] != 0 {
                        L::L420
                    } else if iw[MODE] == 0 {
                        L::L370
                    } else {
                        w[RADIUS] = w[LMAX0];
                        iw[MODE] = 0;
                        L::L140
                    }
                }
            }
            L::L140 => L::L150,
            L::L150 => {
                let k = iw[NITER];
                if k < iw[MXITER] {
                    iw[NITER] = k + 1;
                    L::L160
                } else {
                    iw[1] = 10;
                    L::L430
                }
            }
            L::L160 => {
                let k = iw[NITER] - 1;
                if k != 0 {
                    let t = w[RADFAC] * w[DSTNRM];
                    if w[RADFAC] < 1.0 || t > w[RADIUS] {
                        w[RADIUS] = t;
                    }
                }
                L::L170
            }
            L::L170 => {
                let x01 = iw[X0] as usize;
                w[F0] = w[F];
                iw[IRC] = 4;
                iw[KAGQT] = -1;
                dv7cpy(n, w, x01, lay.x);
                L::L180
            }
            L::L180 => L::L200,
            L::L190 => {
                if w[F] >= w[F0] {
                    L::L200
                } else {
                    w[RADFAC] = 1.0;
                    let k = iw[NITER];
                    iw[NITER] = k + 1;
                    L::L160
                }
            }
            L::L200 => {
                if iw[NFCALL] < iw[MXFCAL] {
                    L::L220
                } else {
                    iw[1] = 9;
                    L::L210
                }
            }
            L::L210 => {
                if w[F] >= w[F0] {
                    L::L430
                } else {
                    iw[CNVCOD] = iw[1];
                    L::L360
                }
            }
            L::L220 => {
                let step1 = iw[STEP] as usize;
                let dg1 = iw[DG] as usize;
                let nwtst1 = iw[NWTSTP] as usize;
                let w1 = nwtst1 + n as usize;
                let dstep1 = step1 + n as usize;
                let ipi = iw[PERM] as usize;
                let l = iw[LMAT] as usize;
                let tg1 = dg1 + n as usize;
                let x01 = iw[X0] as usize;
                let td1 = x01 + n as usize;
                let mut ka = iw[KAGQT];
                let nc = iw[NC];
                dd7dgb(w, iw, lay.b, lay.d, dg1, dstep1, lay.g, ipi, &mut ka, l, n, nc, nwtst1, step1, td1, tg1, w1, x01);
                iw[KAGQT] = ka;
                if iw[IRC] != 6 {
                    L::L230
                } else if iw[RESTOR] != 2 {
                    L::L250
                } else {
                    rstrst = 2;
                    L::L260
                }
            }
            L::L230 => {
                iw[TOOBIG] = 0;
                if w[DSTNRM] <= 0.0 {
                    L::L250
                } else if iw[IRC] != 5 {
                    L::L240
                } else if w[RADFAC] <= 1.0 {
                    L::L240
                } else if w[PREDUC] > 1.2 * w[FDIF] {
                    L::L240
                } else if iw[RESTOR] != 2 {
                    L::L250
                } else {
                    rstrst = 0;
                    L::L260
                }
            }
            L::L240 => {
                let step1 = iw[STEP] as usize;
                let x01 = iw[X0] as usize;
                dv2axy(n, w, lay.x, 1.0, step1, x01);
                iw[NFCALL] += 1;
                iw[1] = 1;
                L::L440
            }
            L::L250 => {
                rstrst = 3;
                L::L260
            }
            L::L260 => {
                let x01 = iw[X0] as usize;
                w[RELDX] = drldst(n, w, lay.d, lay.x, x01);
                da7sst(iw, w);
                let step1 = iw[STEP] as usize;
                let lstgst = iw[STLSTG] as usize;
                match iw[RESTOR] + 1 {
                    2 => dv7cpy(n, w, lay.x, x01),
                    3 => dv7cpy(n, w, lstgst, lay.x),
                    4 => {
                        dv7cpy(n, w, lay.x, lstgst);
                        dv2axy(n, w, step1, -1.0, x01, lay.x);
                        w[RELDX] = drldst(n, w, lay.d, lay.x, x01);
                        iw[RESTOR] = rstrst;
                    }
                    _ => {}
                }
                let k = iw[IRC];
                match k {
                    1 | 5 => {
                        w[RADIUS] = w[RADFAC] * w[DSTNRM];
                        L::L180
                    }
                    2..=4 => L::L340,
                    6 => {
                        w[RADIUS] = w[LMAXS];
                        L::L220
                    }
                    7..=12 => L::L330,
                    13 => L::L400,
                    14 => L::L370,
                    _ => L::L340,
                }
            }
            L::L330 => {
                let k = iw[IRC];
                iw[CNVCOD] = k - 4;
                if w[F] >= w[F0] || iw[XIRC] == 14 {
                    L::L420
                } else {
                    iw[XIRC] = 14;
                    L::L340
                }
            }
            L::L340 => {
                let mut x01 = iw[X0] as usize;
                let step1 = iw[STEP] as usize;
                dv2axy(n, w, step1, -1.0, x01, lay.x);
                if iw[IRC] == 3 {
                    let mut ipi = iw[PERM] as usize;
                    dv7cpy(n, w, x01, step1);
                    dv7ipr(n, iw, ipi, w, x01);
                    let l = iw[LMAT] as usize;
                    dl7tvm(n, w, x01, l, x01);
                    dl7vml(n, w, x01, l, x01);
                    let temp1 = iw[STLSTG] as usize;
                    let temp0 = temp1 - 1;
                    for _ in 1..=n {
                        let j = iw[ipi];
                        ipi += 1;
                        let k = temp0 + j as usize;
                        w[k] = w[x01];
                        x01 += 1;
                    }
                }
                L::L360
            }
            L::L360 => {
                let g01 = (iw[NWTSTP] + n) as usize;
                dv7cpy(n, w, g01, lay.g);
                iw[NGCALL] += 1;
                iw[TOOBIG] = 0;
                iw[1] = 2;
                return;
            }
            L::L370 => {
                let g01 = (iw[NWTSTP] + n) as usize;
                dv2axy(n, w, g01, -1.0, g01, lay.g);
                let step1 = iw[STEP] as usize;
                let temp1 = iw[STLSTG] as usize;
                if iw[IRC] == 3 {
                    dv2axy(n, w, temp1, -1.0, g01, temp1);
                    dv7vmp(n, w, temp1, temp1, lay.d, -1);
                    if dv2nrm(n, w, temp1) <= w[DGNORM] * w[TUNER4] || dd7tpr(n, w, lay.g, step1) < w[GTSTEP] * w[TUNER5] {
                        w[RADFAC] = w[INCFAC];
                    }
                }
                L::L390
            }
            L::L390 => {
                let g01 = (iw[NWTSTP] + n) as usize;
                let step1 = iw[STEP] as usize;
                let temp1 = iw[STLSTG] as usize;
                let w1 = iw[NWTSTP] as usize;
                let z = iw[X0] as usize;
                let l = iw[LMAT] as usize;
                let ipi = iw[PERM] as usize;
                dv7ipr(n, iw, ipi, w, step1);
                dv7ipr(n, iw, ipi, w, g01);
                dw7zbf(w, l, n, step1, w1, g01, z);
                dl7upd(w, temp1, step1, l, g01, l, n, w1, z);
                iw[1] = 2;
                L::L140
            }
            L::L400 => {
                iw[1] = 64;
                L::L430
            }
            L::L410 => {
                iw[1] = 82;
                L::L430
            }
            L::L420 => {
                iw[1] = iw[CNVCOD];
                iw[CNVCOD] = 0;
                L::L430
            }
            L::L430 => return,
            L::L440 => {
                for i in 1..=n {
                    if w[xat(i)] < w[b1(i)] {
                        w[xat(i)] = w[b1(i)];
                    }
                    if w[xat(i)] > w[b2(i)] {
                        w[xat(i)] = w[b2(i)];
                    }
                }
                return;
            }
        };
    }
}

// --- the minimiser ----------------------------------------------------------------------------

/// The caller's arrays of one `nlminb` problem, placed after V in the workspace.
struct Layout {
    liv: usize,
    lv: usize,
    n: usize,
    b: usize,
    d: usize,
    g: usize,
    h: usize,
    x: usize,
}

/// `DRMNHB(B, D, FX, G, H, IV, LH, LIV, LV, N, V, X)`: one reverse-communication call.
fn drmnhb(w: &mut [f64], iw: &mut [i32], lay: &Layout, fx: f64) {
    use vx::*;
    #[derive(Clone, Copy)]
    enum L {
        L10, L20, L30, L50, L60, L70, L140, L150, L160, L180, L190, L200, L210, L220, L230, L240, L250, L260, L340, L370, L380, L400, L410, L420, L430, L440, L450,
    }
    let n = lay.n as i32;
    let (liv, lv) = (lay.liv, lay.lv);
    let b1 = |i: i32| lay.b + 2 * (i as usize - 1);
    let b2 = |i: i32| lay.b + 2 * (i as usize - 1) + 1;
    let dat = |i: i32| lay.d + i as usize - 1;
    let gat = |i: i32| lay.g + i as usize - 1;
    let hat = |i: i32| lay.h + i as usize - 1;
    let xat = |i: i32| lay.x + i as usize - 1;
    let vat = |i: i32| i as usize; // V(i) is w[i]
    let mut i = iw[1];
    let mut rstrst = 3;
    let mut lab = if i == 1 {
        L::L50
    } else if i == 2 {
        L::L60
    } else {
        if iw[1] == 0 {
            divset(2, iw, liv, lv, w);
        }
        if (12..=13).contains(&iw[1]) {
            iw[VNEED] += n * (n + 27) / 2 + 7;
            iw[IVNEED] += 3 * n;
        }
        L::L10
    };
    let nn1o2 = n * (n + 1) / 2;
    loop {
        lab = match lab {
            L::L10 => {
                dparck(2, lay.d, iw, liv, lv, n, w);
                i = iw[1] - 2;
                if i > 12 {
                    return;
                }
                // LH >= NN1O2 always holds here
                match i {
                    1..=6 => L::L250,
                    7 | 9 => L::L190,
                    8 => L::L150,
                    10 | 11 => L::L20,
                    12 => L::L30,
                    _ => {
                        // the computed select case falls through to the next statement
                        iw[1] = 81;
                        L::L440
                    }
                }
            }
            L::L20 => {
                iw[DTOL] = iw[LMAT] + nn1o2;
                iw[X0] = iw[DTOL] + 2 * n;
                iw[STEP] = iw[X0] + 2 * n;
                iw[DG] = iw[STEP] + 3 * n;
                iw[W] = iw[DG] + 2 * n;
                iw[NEXTV] = iw[W] + 4 * n + 7;
                iw[NEXTIV] = iw[PERM] + 3 * n;
                if iw[1] != 13 {
                    L::L30
                } else {
                    iw[1] = 14;
                    return;
                }
            }
            L::L30 => {
                iw[NITER] = 0;
                iw[NFCALL] = 1;
                iw[NGCALL] = 1;
                iw[NFGCAL] = 1;
                iw[MODE] = -1;
                iw[MODEL] = 1;
                iw[STGLIM] = 1;
                iw[TOOBIG] = 0;
                iw[CNVCOD] = 0;
                iw[RADINC] = 0;
                iw[NC] = n;
                w[RAD0] = 0.0;
                w[STPPAR] = 0.0;
                if w[DINIT] >= 0.0 {
                    let c = w[DINIT];
                    dv7scp(n, w, lay.d, c);
                }
                let mut k = iw[DTOL];
                if w[DTINIT] > 0.0 {
                    let c = w[DTINIT];
                    dv7scp(n, w, vat(k), c);
                }
                k += n;
                if w[D0INIT] > 0.0 {
                    let c = w[D0INIT];
                    dv7scp(n, w, vat(k), c);
                }
                let mut ipi = iw[PERM] as usize;
                let mut bad = false;
                for i in 1..=n {
                    iw[ipi] = i;
                    ipi += 1;
                    if w[b1(i)] > w[b2(i)] {
                        bad = true;
                        break;
                    }
                }
                if bad {
                    L::L420
                } else {
                    iw[1] = 1;
                    L::L450
                }
            }
            L::L50 => {
                w[F] = fx;
                if iw[MODE] >= 0 {
                    L::L250
                } else {
                    w[F0] = fx;
                    iw[1] = 2;
                    if iw[TOOBIG] == 0 {
                        return;
                    }
                    iw[1] = 63;
                    L::L440
                }
            }
            L::L60 => {
                if iw[TOOBIG] != 0 {
                    iw[1] = 65;
                    L::L440
                } else {
                    L::L70
                }
            }
            L::L70 => {
                let dg1 = iw[DG] as usize;
                if iw[DTYPE] > 0 {
                    // DD7DUP: nlminb uses DTYPE = 0
                    unreachable!("DTYPE > 0 is not used by nlminb");
                }
                // 90: scaled gradient and its norm
                dv7vmp(n, w, dg1, lay.g, lay.d, -1);
                // scaled Hessian
                let mut k = 1;
                for i in 1..=n {
                    let t = 1.0 / w[dat(i)];
                    for j in 1..=i {
                        w[hat(k)] = t * w[hat(k)] / w[dat(j)];
                        k += 1;
                    }
                }
                // initial permutation
                let mut ipi = iw[PERM] as usize;
                let mut ipn = ipi + n as usize;
                let ipiv2 = ipn - 1;
                i7pnvr(n, iw, ipn, ipi);
                let k = iw[NC];
                for i in 1..=n {
                    let free = !(w[b1(i)] >= w[b2(i)]) && {
                        let xi = w[xat(i)];
                        let gi = w[gat(i)];
                        !(xi <= w[b1(i)] && gi > 0.0) && !(xi >= w[b2(i)] && gi < 0.0)
                    };
                    if free {
                        iw[ipi] = i;
                        ipi += 1;
                        let j = ipiv2 + i as usize;
                        if iw[j] > k {
                            iw[CNVCOD] = 0;
                        }
                    } else {
                        ipn -= 1;
                        iw[ipn] = i;
                    }
                }
                iw[NC] = ipn as i32 - iw[PERM];
                // permute scaled gradient and Hessian
                let ipi = iw[PERM] as usize;
                ds7ipr(n, iw, ipi, w, lay.h);
                dv7ipr(n, iw, ipi, w, dg1);
                w[DGNORM] = 0.0;
                if iw[NC] > 0 {
                    w[DGNORM] = dv2nrm(iw[NC], w, dg1);
                }
                if iw[CNVCOD] != 0 {
                    L::L430
                } else if iw[MODE] == 0 {
                    L::L380
                } else {
                    w[RADIUS] = w[LMAX0] / (1.0 + w[PHMXFC]);
                    iw[MODE] = 0;
                    L::L140
                }
            }
            L::L140 => L::L150, // DITSUM: no output (OUTLEV = 0)
            L::L150 => {
                let k = iw[NITER];
                if k < iw[MXITER] {
                    iw[NITER] = k + 1;
                    L::L160
                } else {
                    iw[1] = 10;
                    L::L440
                }
            }
            L::L160 => {
                let k = iw[NITER] - 1;
                let x01 = iw[X0] as usize;
                w[F0] = w[F];
                iw[IRC] = 4;
                iw[KAGQT] = -1;
                dv7cpy(n, w, x01, lay.x);
                if k != 0 {
                    let step1 = iw[STEP] as usize;
                    let mut kk = step1;
                    for i in 1..=n {
                        w[kk] *= w[dat(i)];
                        kk += 1;
                    }
                    let t = w[RADFAC] * dv2nrm(n, w, step1);
                    if w[RADFAC] < 1.0 || t > w[RADIUS] {
                        w[RADIUS] = t;
                    }
                }
                L::L180
            }
            L::L180 => L::L200, // STOPX() is false
            L::L190 => {
                if w[F] >= w[F0] {
                    L::L200
                } else {
                    w[RADFAC] = 1.0;
                    // K = IV(NITER); GO TO 160 (which sets NITER = K + 1)
                    let k = iw[NITER];
                    iw[NITER] = k + 1;
                    L::L160
                }
            }
            L::L200 => {
                if iw[NFCALL] < iw[MXFCAL] {
                    L::L220
                } else {
                    iw[1] = 9;
                    L::L210
                }
            }
            L::L210 => {
                if w[F] >= w[F0] {
                    L::L440
                } else {
                    iw[CNVCOD] = iw[1];
                    L::L370
                }
            }
            L::L220 => {
                let step1 = iw[STEP] as usize;
                let l = iw[LMAT] as usize;
                let w1 = iw[W] as usize;
                let ipi = iw[PERM] as usize;
                let ipn = ipi + n as usize;
                let ipiv2 = ipn + n as usize;
                let tg1 = iw[DG] as usize;
                let td1 = tg1 + n as usize;
                let x01 = iw[X0] as usize;
                let x11 = x01 + n as usize;
                let mut ka = iw[KAGQT];
                let mut p0 = iw[N0];
                let nc = iw[NC];
                dg7qsb(w, iw, lay.b, lay.d, lay.h, lay.g, ipi, ipn, ipiv2, &mut ka, l, n, &mut p0, nc, step1, td1, tg1, w1, x11, x01);
                iw[KAGQT] = ka;
                iw[N0] = p0;
                if iw[IRC] != 6 {
                    L::L230
                } else if iw[RESTOR] != 2 {
                    L::L250
                } else {
                    rstrst = 2;
                    L::L260
                }
            }
            L::L230 => {
                iw[TOOBIG] = 0;
                if w[DSTNRM] <= 0.0 {
                    L::L250
                } else if iw[IRC] != 5 {
                    L::L240
                } else if w[RADFAC] <= 1.0 {
                    L::L240
                } else if w[PREDUC] > 1.2 * w[FDIF] {
                    L::L240
                } else if iw[RESTOR] != 2 {
                    L::L250
                } else {
                    rstrst = 0;
                    L::L260
                }
            }
            L::L240 => {
                let step1 = iw[STEP] as usize;
                let x01 = iw[X0] as usize;
                dv2axy(n, w, lay.x, 1.0, step1, x01);
                iw[NFCALL] += 1;
                iw[1] = 1;
                L::L450
            }
            L::L250 => {
                rstrst = 3;
                L::L260
            }
            L::L260 => {
                let x01 = iw[X0] as usize;
                w[RELDX] = drldst(n, w, lay.d, lay.x, x01);
                da7sst(iw, w);
                let step1 = iw[STEP] as usize;
                let lstgst = step1 + 2 * n as usize;
                match iw[RESTOR] + 1 {
                    2 => dv7cpy(n, w, lay.x, x01),
                    3 => dv7cpy(n, w, lstgst, lay.x),
                    4 => {
                        dv7cpy(n, w, lay.x, lstgst);
                        dv2axy(n, w, step1, -1.0, x01, lay.x);
                        w[RELDX] = drldst(n, w, lay.d, lay.x, x01);
                        iw[RESTOR] = rstrst;
                    }
                    _ => {}
                }
                // 300
                let k = iw[IRC];
                match k {
                    1 | 5 => {
                        w[RADIUS] = w[RADFAC] * w[DSTNRM];
                        L::L180
                    }
                    2..=4 => L::L340,
                    6 => {
                        w[RADIUS] = w[LMAXS];
                        L::L220
                    }
                    7..=12 => {
                        iw[CNVCOD] = k - 4;
                        if w[F] >= w[F0] || iw[XIRC] == 14 {
                            L::L430
                        } else {
                            iw[XIRC] = 14;
                            L::L340
                        }
                    }
                    13 => L::L410,
                    14 => L::L380,
                    _ => L::L340, // out-of-range computed GO TO continues at 310... not reached
                }
            }
            L::L340 => {
                if iw[IRC] == 3 {
                    let step1 = iw[STEP] as usize;
                    let temp1 = step1 + 2 * n as usize;
                    let mut x01 = iw[X0] as usize;
                    // TEMP1 = HESSIAN * STEP + G(X0)
                    let mut k = temp1;
                    let step0 = step1 - 1;
                    let mut ipi = iw[PERM] as usize;
                    for _ in 1..=n {
                        let j = iw[ipi];
                        ipi += 1;
                        let s1 = step0 + j as usize;
                        w[k] = w[dat(j)] * w[s1];
                        k += 1;
                    }
                    ds7lvm(n, w, x01, lay.h, temp1);
                    let temp0 = temp1 - 1;
                    let mut ipi = iw[PERM] as usize;
                    for _ in 1..=n {
                        let j = iw[ipi];
                        ipi += 1;
                        let t1 = temp0 + j as usize;
                        w[t1] = w[dat(j)] * w[x01] + w[gat(j)];
                        x01 += 1;
                    }
                }
                L::L370
            }
            L::L370 => {
                iw[NGCALL] += 1;
                iw[TOOBIG] = 0;
                iw[1] = 2;
                L::L450
            }
            L::L380 => {
                iw[1] = 2;
                if iw[IRC] != 3 {
                    L::L140
                } else {
                    let step1 = iw[STEP] as usize;
                    let temp1 = step1 + 2 * n as usize;
                    let mut k = temp1;
                    for i in 1..=n {
                        w[k] = (w[k] - w[gat(i)]) / w[dat(i)];
                        k += 1;
                    }
                    if dv2nrm(n, w, temp1) <= w[DGNORM] * w[TUNER4] {
                        L::L400
                    } else if dd7tpr(n, w, lay.g, step1) >= w[GTSTEP] * w[TUNER5] {
                        L::L140
                    } else {
                        L::L400
                    }
                }
            }
            L::L400 => {
                w[RADFAC] = w[INCFAC];
                L::L140
            }
            L::L410 => {
                iw[1] = 64;
                L::L440
            }
            L::L420 => {
                iw[1] = 82;
                L::L440
            }
            L::L430 => {
                iw[1] = iw[CNVCOD];
                iw[CNVCOD] = 0;
                L::L440
            }
            L::L440 => return, // DITSUM: no output
            L::L450 => {
                for i in 1..=n {
                    if w[xat(i)] < w[b1(i)] {
                        w[xat(i)] = w[b1(i)];
                    }
                    if w[xat(i)] > w[b2(i)] {
                        w[xat(i)] = w[b2(i)];
                    }
                }
                return;
            }
        };
    }
}

/// The result of [`nlminb`]: R's `list(par, objective, convergence, iterations, evaluations,
/// message)`; `code` is `iv[1]`, from which R's message is made.
#[derive(Debug, Clone, PartialEq)]
pub struct NlminbResult {
    pub par: Vec<f64>,
    pub objective: f64,
    /// 0 when `code` is 3 to 6 (convergence), else 1
    pub convergence: i32,
    pub iterations: i32,
    pub evaluations_function: i32,
    pub evaluations_gradient: i32,
    pub code: i32,
}

/// The objective, gradient and Hessian of an `nlminb` problem.
pub trait Problem {
    fn objective(&mut self, x: &[f64]) -> f64;
    fn gradient(&mut self, x: &[f64], g: &mut [f64]);
    /// The full symmetric Hessian, column-major `n x n` (not called when `nlminb` runs without
    /// one).
    fn hessian(&mut self, x: &[f64], h: &mut [f64]);
}

/// `nlminb(start, objective, gradient, hessian, lower, upper, control = list(iter.max =
/// iter_max))` with scale 1 and the other controls at their defaults, as R 4.5.3 runs it
/// (`nlminb.R`, `port_nlminb` in `port.c`): with the Hessian (`DRMNHB`) when `use_hessian`,
/// else with the gradient only (`DRMNGB`). Bounds must be given (R uses other routines without).
pub fn nlminb(start: &[f64], prob: &mut dyn Problem, lower: &[f64], upper: &[f64], iter_max: Option<i32>, use_hessian: bool) -> Result<NlminbResult, PortError> {
    let n = start.len();
    assert!(n >= 1 && lower.len() == n && upper.len() == n);
    let liv = 78 + 3 * n;
    let lv = 130 + (n * (n + 27)) / 2;
    let lh = n * (n + 1) / 2;
    let lay = Layout { liv, lv, n, b: lv + 1, d: lv + 1 + 2 * n, g: lv + 1 + 3 * n, h: lv + 1 + 4 * n, x: lv + 1 + 4 * n + lh };
    let mut w = vec![0.0f64; lay.x + n];
    let mut iw = vec![0i32; liv + 1];
    divset(2, &mut iw, liv, lv, &mut w);
    if let Some(m) = iter_max {
        iw[MXITER] = m;
    }
    for i in 0..n {
        w[lay.b + 2 * i] = lower[i];
        w[lay.b + 2 * i + 1] = upper[i];
        w[lay.d + i] = 1.0;
        w[lay.x + i] = start[i];
    }
    let mut fx = f64::INFINITY;
    let mut g = vec![0.0; n];
    let mut hfull = vec![0.0; n * n];
    loop {
        if use_hessian {
            drmnhb(&mut w, &mut iw, &lay, fx);
        } else {
            drmngb(&mut w, &mut iw, &lay, fx);
        }
        let x: Vec<f64> = w[lay.x..lay.x + n].to_vec();
        if iw[1] == 2 {
            prob.gradient(&x, &mut g);
            if g.iter().any(|v| v.is_nan()) {
                return Err(PortError::NanDerivative("gradient"));
            }
            w[lay.g..lay.g + n].copy_from_slice(&g);
            if !use_hessian {
                if iw[1] >= 3 {
                    break;
                }
                continue;
            }
            prob.hessian(&x, &mut hfull);
            let mut pos = 0;
            for i in 0..n {
                for j in 0..=i {
                    let v = hfull[i + j * n];
                    if v.is_nan() {
                        return Err(PortError::NanDerivative("Hessian"));
                    }
                    w[lay.h + pos] = v;
                    pos += 1;
                }
            }
        } else {
            fx = prob.objective(&x);
            if fx.is_nan() {
                fx = f64::INFINITY;
            }
        }
        if iw[1] >= 3 {
            break;
        }
    }
    let code = iw[1];
    Ok(NlminbResult {
        par: w[lay.x..lay.x + n].to_vec(),
        objective: w[vx::F],
        convergence: if (3..=6).contains(&code) { 0 } else { 1 },
        iterations: iw[NITER],
        evaluations_function: iw[NFCALL],
        evaluations_gradient: iw[NGCALL],
        code,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Quad;
    impl Problem for Quad {
        fn objective(&mut self, x: &[f64]) -> f64 {
            (x[0] - 3.0) * (x[0] - 3.0)
        }
        fn gradient(&mut self, x: &[f64], g: &mut [f64]) {
            g[0] = 2.0 * (x[0] - 3.0);
        }
        fn hessian(&mut self, _x: &[f64], h: &mut [f64]) {
            h[0] = 2.0;
        }
    }

    #[test]
    fn minimises_a_quadratic() {
        for h in [true, false] {
            let r = nlminb(&[0.0], &mut Quad, &[-10.0], &[10.0], None, h).unwrap();
            assert!((r.par[0] - 3.0).abs() < 1e-8, "{r:?}");
            assert_eq!(r.convergence, 0);
            let r = nlminb(&[0.0], &mut Quad, &[-10.0], &[1.0], None, h).unwrap();
            assert_eq!(r.par[0], 1.0, "{r:?}");
        }
    }
}

#[cfg(test)]
mod bfgs_tests {
    use super::*;

    /// L+ L+^T from DW7ZBF/DL7UPD must be the BFGS update of H = L L^T.
    #[test]
    fn bfgs_update_is_bfgs() {
        // w layout (1-based): L at 1..=3, s at 4..=5, y at 6..=7, w at 8..=9, z at 10..=11,
        // lambda at 12..=13, beta at 14..=15, gamma at 16..=17
        let mut w = vec![0.0; 18];
        let (l, s, y, wv, z, lam, beta, gamma) = (1, 4, 6, 8, 10, 12, 14, 16);
        w[1] = 2.0; w[2] = 0.5; w[3] = 1.5; // L = [[2, 0], [0.5, 1.5]]
        w[4] = 0.3; w[5] = -0.2; // s
        w[6] = 1.1; w[7] = 0.4; // y
        let h = |w: &[f64]| {
            let (a, b, c) = (w[1], w[2], w[3]);
            [a * a, a * b, a * b, b * b + c * c]
        };
        let h0 = h(&w);
        let (s0, y0) = ([w[4], w[5]], [w[6], w[7]]);
        dw7zbf(&mut w, l, 2, s, wv, y, z);
        dl7upd(&mut w, beta, gamma, l, lam, l, 2, wv, z);
        let h1 = h(&w);
        // BFGS: H+ = H - H s s^T H / (s^T H s) + y y^T / (y^T s)
        let hs = [h0[0] * s0[0] + h0[2] * s0[1], h0[1] * s0[0] + h0[3] * s0[1]];
        let shs = s0[0] * hs[0] + s0[1] * hs[1];
        let ys = y0[0] * s0[0] + y0[1] * s0[1];
        let want = [
            h0[0] - hs[0] * hs[0] / shs + y0[0] * y0[0] / ys,
            h0[1] - hs[1] * hs[0] / shs + y0[1] * y0[0] / ys,
            h0[2] - hs[0] * hs[1] / shs + y0[0] * y0[1] / ys,
            h0[3] - hs[1] * hs[1] / shs + y0[1] * y0[1] / ys,
        ];
        for k in 0..4 {
            assert!((h1[k] - want[k]).abs() < 1e-12, "{h1:?} vs {want:?}");
        }
    }
}
