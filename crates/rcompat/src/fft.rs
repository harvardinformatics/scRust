//! R's `fft()`: Singleton's mixed-radix fast Fourier transform as R 4.5.3 implements it
//! (`fft_factor`, `fft_work` and `fftmx` in `src/library/stats/src/fft.c`, called from `fft()` in
//! `fourier.c`), ported statement by statement so that every result is bit-identical to R's.
//!
//! R transforms a complex vector in place, with the real and imaginary parts interleaved and the
//! index increment `|isn| = 2` (`isn = -2` forward, `+2` inverse). The port keeps that layout and
//! that increment, so all index arithmetic (including its integer divisions) is R's. The inverse
//! transform is unnormalised, as in R. `sin` and `cos` are glibc's, as R calls them.
//!
//! The C code is a translation of Fortran with `goto`s; here each label is a state of a loop
//! over [`Label`], and the statements between labels are copied in order.

use crate::libm::{cos, sin};

/// A complex number, as R stores one (`Rcomplex`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const fn new(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }

}

impl std::ops::Mul for Complex {
    type Output = Complex;

    /// C99 complex multiplication as GCC compiles it without FMA (R's `*` on complex vectors,
    /// `src/main/complex.c`): `(ac - bd) + (ad + bc)i`. When both parts come out NaN, C's
    /// Annex G recovers infinities through `__muldc3`; that case is not reproduced.
    fn mul(self, other: Complex) -> Complex {
        let (a, b, c, d) = (self.re, self.im, other.re, other.im);
        let (ac, bd, ad, bc) = (a * c, b * d, a * d, b * c);
        Complex { re: ac - bd, im: ad + bc }
    }
}

/// `fft_factor` failed: n has more than 20 factors (impossible for n < 2^31).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FftError;

impl std::fmt::Display for FftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("fft factorization error")
    }
}

impl std::error::Error for FftError {}

/// The factorization `fft_factor` computes, which `fftmx` then consumes.
#[derive(Debug, Clone)]
struct Factors {
    nfac: [i32; 20],
    m_fac: i32,
    kt: i32,
    maxf: i32,
    maxp: i32,
}

/// `fft_factor(n, &maxf, &maxp)` for n >= 2: square factors first (16 = 4^2 separately, then
/// odd squares), then the square-free rest; the square factors are repeated at the end.
fn fft_factor(n: i32) -> Result<Factors, FftError> {
    debug_assert!(n >= 2);
    let mut f = Factors { nfac: [0; 20], m_fac: 0, kt: 0, maxf: 0, maxp: 0 };
    let mut k = n;
    let push = |f: &mut Factors, v: i32| -> Result<(), FftError> {
        if f.m_fac as usize >= f.nfac.len() {
            return Err(FftError);
        }
        f.nfac[f.m_fac as usize] = v;
        f.m_fac += 1;
        Ok(())
    };
    // extract 4^2 = 16 separately: at most one remaining factor 2^2 = 4, done below
    while k % 16 == 0 {
        push(&mut f, 4)?;
        k /= 16;
    }
    // extract 3^2, 5^2, ...
    let mut kchanged = false;
    let mut sqrtk = (k as f64).sqrt() as i32;
    let mut j = 3;
    while j <= sqrtk {
        let jj = j * j;
        while k % jj == 0 {
            push(&mut f, j)?;
            k /= jj;
            kchanged = true;
        }
        if kchanged {
            kchanged = false;
            sqrtk = (k as f64).sqrt() as i32;
        }
        j += 2;
    }
    if k <= 4 {
        f.kt = f.m_fac;
        f.nfac[f.m_fac as usize] = k;
        if k != 1 {
            f.m_fac += 1;
        }
    } else {
        if k % 4 == 0 {
            push(&mut f, 2)?;
            k /= 4;
        }
        // all square factors out now, but k >= 5 still
        f.kt = f.m_fac;
        f.maxp = (f.kt + f.kt + 2).max(k - 1);
        let mut j = 2;
        loop {
            if k % j == 0 {
                push(&mut f, j)?;
                k /= j;
            }
            if j > i32::MAX - 2 {
                break;
            }
            j = ((j + 1) / 2) * 2 + 1;
            if j > k {
                break;
            }
        }
    }
    if f.m_fac <= f.kt + 1 {
        f.maxp = f.m_fac + f.kt + 1;
    }
    if f.m_fac + f.kt > 20 {
        return Err(FftError);
    }
    if f.kt != 0 {
        let mut j = f.kt;
        while j != 0 {
            j -= 1;
            let v = f.nfac[j as usize];
            push(&mut f, v)?;
        }
    }
    let nf = |i: i32| f.nfac[i as usize];
    let mut maxf = nf(f.m_fac - f.kt - 1);
    // the last squared factor is not necessarily the largest (PR#1429)
    if f.kt > 0 {
        maxf = maxf.max(nf(f.kt - 1));
    }
    if f.kt > 1 {
        maxf = maxf.max(nf(f.kt - 2));
    }
    if f.kt > 2 {
        maxf = maxf.max(nf(f.kt - 3));
    }
    f.maxf = maxf;
    Ok(f)
}

/// The labels of `fftmx` that are jumped to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Label {
    Start,
    L60,
    L70,
    L80,
    L100,
    L110,
    L120,
    L130,
    L140,
    L150,
    L160,
    L170,
    L180,
    L190,
    L200,
    F5,
    L220,
    FOdd,
    L250,
    L260,
    L270,
    L290,
    L300,
    Fin,
    L370,
    L380,
    L390,
    L440,
    L460,
    L470,
    L480,
    L490,
    L500,
    Ord,
    L520,
    L570,
    Return,
}

/// `fftmx(a, b, ntot, n, nspan, isn, m, kt, at, ck, bt, sk, np, nfac)` for one variate
/// (`nseg = nspn = 1`, so `ntot = nspan = n`), on the interleaved array `z` with `a = z` and
/// `b = z + 1`.
#[allow(unused_assignments, clippy::cognitive_complexity, clippy::assign_op_pattern)]
fn fftmx(z: &mut [f64], n: i32, isn: i32, fac: &Factors) {
    let ntot = n;
    let nspan = n;
    let m = fac.m_fac;
    let mut kt = fac.kt;
    let mut nfac = fac.nfac;
    let mf = fac.maxf as usize;
    // 1-based work arrays: at, ck, bt, sk of maxf entries, np of maxp entries
    let mut at = vec![0.0f64; mf + 1];
    let mut ck = vec![0.0f64; mf + 1];
    let mut bt = vec![0.0f64; mf + 1];
    let mut sk = vec![0.0f64; mf + 1];
    let mut np = vec![0i32; fac.maxp.max(1) as usize + 1];

    // a[k] = z[k - 1], b[k] = z[k] (1-based, after `a--; b--;`)
    macro_rules! a {
        ($k:expr) => {
            z[($k - 1) as usize]
        };
    }
    macro_rules! b {
        ($k:expr) => {
            z[($k) as usize]
        };
    }
    macro_rules! nf {
        ($k:expr) => {
            nfac[($k) as usize]
        };
    }
    macro_rules! npx {
        ($k:expr) => {
            np[($k) as usize]
        };
    }

    // the C function's locals, all zero-initialised here
    let (mut aa, mut aj, mut ajm, mut ajp, mut ak, mut akm, mut akp) = (0.0f64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut bb, mut bj, mut bjm, mut bjp, mut bk, mut bkm, mut bkp) = (0.0f64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut c1, mut c2, mut c3, mut cd) = (0.0f64, 0.0, 0.0, 0.0);
    let (mut s1, mut s2, mut s3, mut sd) = (0.0f64, 0.0, 0.0, 0.0);
    let mut dr = 0.0f64;
    let (mut i, mut j, mut jc, mut jf, mut jj) = (0i32, 0i32, 0i32, 0i32, 0i32);
    let (mut k, mut k1, mut k2, mut k3, mut k4) = (0i32, 0i32, 0i32, 0i32, 0i32);
    let (mut kk, mut kspan, mut kspnn) = (0i32, 0i32, 0i32);
    let (mut maxf, mut mm, mut nn, mut nt) = (0i32, 0i32, 0i32, 0i32);

    let inc = isn.abs();
    nt = inc * ntot;
    let ks = inc * nspan;
    let mut rad = std::f64::consts::FRAC_PI_4; // M_PI_4: 45 degrees
    let mut s72 = rad / 0.625; // 72 = 45 / .625 degrees
    let c72 = cos(s72);
    s72 = sin(s72);
    let mut s120 = 0.5 * M_SQRT_3; // sin(120) = sqrt(3)/2
    if isn <= 0 {
        s72 = -s72;
        s120 = -s120;
        rad = -rad;
    }
    // (SCALING is not defined in R: no scaling by 1/n for isn > 0)

    kspan = ks;
    nn = nt - inc;
    jc = ks / n;

    // sin, cos values are re-initialized each lim steps
    let lim = 32;
    let klim = lim * jc;
    i = 0;
    jf = 0;
    maxf = nf!(m - kt - 1);
    if kt > 0 {
        maxf = maxf.max(nf!(kt - 1));
    }

    let mut label = Label::Start;
    loop {
        match label {
            Label::Start => {
                dr = (8.0 * jc as f64) / kspan as f64;
                cd = sin(0.5 * dr * rad);
                cd = 2.0 * cd * cd;
                sd = sin(dr * rad);
                kk = 1;
                i += 1;
                if nf!(i - 1) != 2 {
                    label = Label::L110;
                    continue;
                }
                // transform for factor of 2 (including rotation factor)
                kspan /= 2;
                k1 = kspan + 2;
                loop {
                    loop {
                        k2 = kk + kspan;
                        ak = a!(k2);
                        bk = b!(k2);
                        a!(k2) = a!(kk) - ak;
                        b!(k2) = b!(kk) - bk;
                        a!(kk) += ak;
                        b!(kk) += bk;
                        kk = k2 + kspan;
                        if kk > nn {
                            break;
                        }
                    }
                    kk -= nn;
                    if kk > jc {
                        break;
                    }
                }
                label = if kk > kspan { Label::Fin } else { Label::L60 };
            }
            Label::L60 => {
                c1 = 1.0 - cd;
                s1 = sd;
                mm = (k1 / 2).min(klim);
                label = Label::L80;
            }
            Label::L70 => {
                ak = c1 - (cd * c1 + sd * s1);
                s1 = (sd * c1 - cd * s1) + s1;
                // (TRUNCATED_ARITHMETIC is not defined: c1 = ak)
                c1 = ak;
                label = Label::L80;
            }
            Label::L80 => {
                loop {
                    k2 = kk + kspan;
                    ak = a!(kk) - a!(k2);
                    bk = b!(kk) - b!(k2);
                    a!(kk) += a!(k2);
                    b!(kk) += b!(k2);
                    a!(k2) = c1 * ak - s1 * bk;
                    b!(k2) = s1 * ak + c1 * bk;
                    kk = k2 + kspan;
                    if kk >= nt {
                        break;
                    }
                }
                k2 = kk - nt;
                c1 = -c1;
                kk = k1 - k2;
                if kk > k2 {
                    label = Label::L80;
                    continue;
                }
                kk += jc;
                if kk <= mm {
                    label = Label::L70;
                    continue;
                }
                if kk >= k2 {
                    k1 += inc + inc;
                    kk = (k1 - kspan) / 2 + jc;
                    label = if kk <= jc + jc { Label::L60 } else { Label::Start };
                    continue;
                }
                s1 = ((kk - 1) / jc) as f64 * dr * rad;
                c1 = cos(s1);
                s1 = sin(s1);
                mm = (k1 / 2).min(mm + klim);
                label = Label::L80;
            }
            Label::L100 => {
                // transform for factor of 3 (optional code)
                k1 = kk + kspan;
                k2 = k1 + kspan;
                ak = a!(kk);
                bk = b!(kk);
                aj = a!(k1) + a!(k2);
                bj = b!(k1) + b!(k2);
                a!(kk) = ak + aj;
                b!(kk) = bk + bj;
                ak = -0.5 * aj + ak;
                bk = -0.5 * bj + bk;
                aj = (a!(k1) - a!(k2)) * s120;
                bj = (b!(k1) - b!(k2)) * s120;
                a!(k1) = ak - bj;
                b!(k1) = bk + aj;
                a!(k2) = ak + bj;
                b!(k2) = bk - aj;
                kk = k2 + kspan;
                if kk < nn {
                    label = Label::L100;
                    continue;
                }
                kk -= nn;
                label = if kk <= kspan { Label::L100 } else { Label::L290 };
            }
            Label::L110 => {
                // transform for factor of 4
                if nf!(i - 1) != 4 {
                    label = Label::FOdd;
                    continue;
                }
                kspnn = kspan;
                kspan /= 4;
                label = Label::L120;
            }
            Label::L120 => {
                c1 = 1.0;
                s1 = 0.0;
                mm = kspan.min(klim);
                label = Label::L150;
            }
            Label::L130 => {
                c2 = c1 - (cd * c1 + sd * s1);
                s1 = (sd * c1 - cd * s1) + s1;
                c1 = c2;
                label = Label::L140;
            }
            Label::L140 => {
                c2 = c1 * c1 - s1 * s1;
                s2 = c1 * s1 * 2.0;
                c3 = c2 * c1 - s2 * s1;
                s3 = c2 * s1 + s2 * c1;
                label = Label::L150;
            }
            Label::L150 => {
                k1 = kk + kspan;
                k2 = k1 + kspan;
                k3 = k2 + kspan;
                akp = a!(kk) + a!(k2);
                akm = a!(kk) - a!(k2);
                ajp = a!(k1) + a!(k3);
                ajm = a!(k1) - a!(k3);
                a!(kk) = akp + ajp;
                ajp = akp - ajp;
                bkp = b!(kk) + b!(k2);
                bkm = b!(kk) - b!(k2);
                bjp = b!(k1) + b!(k3);
                bjm = b!(k1) - b!(k3);
                b!(kk) = bkp + bjp;
                bjp = bkp - bjp;
                if isn < 0 {
                    label = Label::L180;
                    continue;
                }
                akp = akm - bjm;
                akm += bjm;
                bkp = bkm + ajm;
                bkm -= ajm;
                label = if s1 == 0.0 { Label::L190 } else { Label::L160 };
            }
            Label::L160 => {
                a!(k1) = akp * c1 - bkp * s1;
                b!(k1) = akp * s1 + bkp * c1;
                a!(k2) = ajp * c2 - bjp * s2;
                b!(k2) = ajp * s2 + bjp * c2;
                a!(k3) = akm * c3 - bkm * s3;
                b!(k3) = akm * s3 + bkm * c3;
                kk = k3 + kspan;
                label = if kk <= nt { Label::L150 } else { Label::L170 };
            }
            Label::L170 => {
                kk = kk - nt + jc;
                if kk <= mm {
                    label = Label::L130;
                    continue;
                }
                if kk < kspan {
                    label = Label::L200;
                    continue;
                }
                kk = kk - kspan + inc;
                if kk <= jc {
                    label = Label::L120;
                    continue;
                }
                label = if kspan == jc { Label::Fin } else { Label::Start };
            }
            Label::L180 => {
                akp = akm + bjm;
                akm -= bjm;
                bkp = bkm - ajm;
                bkm += ajm;
                label = if s1 != 0.0 { Label::L160 } else { Label::L190 };
            }
            Label::L190 => {
                a!(k1) = akp;
                b!(k1) = bkp;
                a!(k2) = ajp;
                b!(k2) = bjp;
                a!(k3) = akm;
                b!(k3) = bkm;
                kk = k3 + kspan;
                label = if kk <= nt { Label::L150 } else { Label::L170 };
            }
            Label::L200 => {
                s1 = ((kk - 1) / jc) as f64 * dr * rad;
                c1 = cos(s1);
                s1 = sin(s1);
                mm = kspan.min(mm + klim);
                label = Label::L140;
            }
            Label::F5 => {
                // transform for factor of 5 (optional code)
                c2 = c72 * c72 - s72 * s72;
                s2 = 2.0 * c72 * s72;
                label = Label::L220;
            }
            Label::L220 => {
                k1 = kk + kspan;
                k2 = k1 + kspan;
                k3 = k2 + kspan;
                k4 = k3 + kspan;
                akp = a!(k1) + a!(k4);
                akm = a!(k1) - a!(k4);
                bkp = b!(k1) + b!(k4);
                bkm = b!(k1) - b!(k4);
                ajp = a!(k2) + a!(k3);
                ajm = a!(k2) - a!(k3);
                bjp = b!(k2) + b!(k3);
                bjm = b!(k2) - b!(k3);
                aa = a!(kk);
                bb = b!(kk);
                a!(kk) = aa + akp + ajp;
                b!(kk) = bb + bkp + bjp;
                ak = akp * c72 + ajp * c2 + aa;
                bk = bkp * c72 + bjp * c2 + bb;
                aj = akm * s72 + ajm * s2;
                bj = bkm * s72 + bjm * s2;
                a!(k1) = ak - bj;
                a!(k4) = ak + bj;
                b!(k1) = bk + aj;
                b!(k4) = bk - aj;
                ak = akp * c2 + ajp * c72 + aa;
                bk = bkp * c2 + bjp * c72 + bb;
                aj = akm * s2 - ajm * s72;
                bj = bkm * s2 - bjm * s72;
                a!(k2) = ak - bj;
                a!(k3) = ak + bj;
                b!(k2) = bk + aj;
                b!(k3) = bk - aj;
                kk = k4 + kspan;
                if kk < nn {
                    label = Label::L220;
                    continue;
                }
                kk -= nn;
                label = if kk <= kspan { Label::L220 } else { Label::L290 };
            }
            Label::FOdd => {
                // transform for odd factors
                k = nf!(i - 1);
                kspnn = kspan;
                kspan /= k;
                if k == 3 {
                    label = Label::L100;
                    continue;
                }
                if k == 5 {
                    label = Label::F5;
                    continue;
                }
                if k == jf {
                    label = Label::L250;
                    continue;
                }
                jf = k;
                s1 = rad / (k as f64 / 8.0);
                c1 = cos(s1);
                s1 = sin(s1);
                ck[jf as usize] = 1.0;
                sk[jf as usize] = 0.0;
                j = 1;
                while j < k {
                    // k is changing as well
                    ck[j as usize] = ck[k as usize] * c1 + sk[k as usize] * s1;
                    sk[j as usize] = ck[k as usize] * s1 - sk[k as usize] * c1;
                    k -= 1;
                    ck[k as usize] = ck[j as usize];
                    sk[k as usize] = -sk[j as usize];
                    j += 1;
                }
                label = Label::L250;
            }
            Label::L250 => {
                k1 = kk;
                k2 = kk + kspnn;
                aa = a!(kk);
                bb = b!(kk);
                ak = aa;
                bk = bb;
                j = 1;
                k1 += kspan;
                label = Label::L260;
            }
            Label::L260 => {
                k2 -= kspan;
                j += 1;
                at[j as usize] = a!(k1) + a!(k2);
                ak = at[j as usize] + ak;
                bt[j as usize] = b!(k1) + b!(k2);
                bk = bt[j as usize] + bk;
                j += 1;
                at[j as usize] = a!(k1) - a!(k2);
                bt[j as usize] = b!(k1) - b!(k2);
                k1 += kspan;
                if k1 < k2 {
                    label = Label::L260;
                    continue;
                }
                a!(kk) = ak;
                b!(kk) = bk;
                k1 = kk;
                k2 = kk + kspnn;
                j = 1;
                label = Label::L270;
            }
            Label::L270 => {
                k1 += kspan;
                k2 -= kspan;
                jj = j;
                ak = aa;
                bk = bb;
                aj = 0.0;
                bj = 0.0;
                k = 2;
                while k < jf {
                    ak += at[k as usize] * ck[jj as usize];
                    bk += bt[k as usize] * ck[jj as usize];
                    k += 1;
                    aj += at[k as usize] * sk[jj as usize];
                    bj += bt[k as usize] * sk[jj as usize];
                    jj += j;
                    if jj > jf {
                        jj -= jf;
                    }
                    k += 1;
                }
                k = jf - j;
                a!(k1) = ak - bj;
                b!(k1) = bk + aj;
                a!(k2) = ak + bj;
                b!(k2) = bk - aj;
                j += 1;
                if j < k {
                    label = Label::L270;
                    continue;
                }
                kk += kspnn;
                if kk <= nn {
                    label = Label::L250;
                    continue;
                }
                kk -= nn;
                label = if kk <= kspan { Label::L250 } else { Label::L290 };
            }
            Label::L290 => {
                // multiply by rotation factor (except for factors of 2 and 4)
                if i == m {
                    label = Label::Fin;
                    continue;
                }
                kk = jc + 1;
                label = Label::L300;
            }
            Label::L300 => {
                c2 = 1.0 - cd;
                s1 = sd;
                mm = kspan.min(klim);
                label = loop {
                    // L320
                    c1 = c2;
                    s2 = s1;
                    kk += kspan;
                    loop {
                        // L330
                        loop {
                            ak = a!(kk);
                            a!(kk) = c2 * ak - s2 * b!(kk);
                            b!(kk) = s2 * ak + c2 * b!(kk);
                            kk += kspnn;
                            if kk > nt {
                                break;
                            }
                        }
                        ak = s1 * s2;
                        s2 = s1 * c2 + c1 * s2;
                        c2 = c1 * c2 - ak;
                        kk += -nt + kspan;
                        if kk > kspnn {
                            break;
                        }
                    }
                    kk += -kspnn + jc;
                    if kk <= mm {
                        // L310
                        c2 = c1 - (cd * c1 + sd * s1);
                        s1 += sd * c1 - cd * s1;
                        continue;
                    }
                    if kk >= kspan {
                        kk = kk - kspan + jc + inc;
                        break if kk <= jc + jc { Label::L300 } else { Label::Start };
                    }
                    s1 = ((kk - 1) / jc) as f64 * dr * rad;
                    c2 = cos(s1);
                    s1 = sin(s1);
                    mm = kspan.min(mm + klim);
                };
            }
            Label::Fin => {
                // permute the results to normal order, done in two stages;
                // permutation for square factors of n
                npx!(1) = ks;
                if kt == 0 {
                    label = Label::L440;
                    continue;
                }
                k = kt + kt + 1;
                if m < k {
                    k -= 1;
                }
                npx!(k + 1) = jc;
                j = 1;
                while j < k {
                    npx!(j + 1) = npx!(j) / nf!(j - 1);
                    npx!(k) = npx!(k + 1) * nf!(j - 1);
                    j += 1;
                    k -= 1;
                }
                k3 = npx!(k + 1);
                kspan = npx!(2);
                kk = jc + 1;
                k2 = kspan + 1;
                j = 1;
                // n == ntot: permutation for a single-variate transform
                label = Label::L370;
            }
            Label::L370 => {
                loop {
                    ak = a!(kk);
                    a!(kk) = a!(k2);
                    a!(k2) = ak;
                    bk = b!(kk);
                    b!(kk) = b!(k2);
                    b!(k2) = bk;
                    kk += inc;
                    k2 += kspan;
                    if k2 >= ks {
                        break;
                    }
                }
                label = Label::L380;
            }
            Label::L380 => {
                loop {
                    k2 -= npx!(j);
                    j += 1;
                    k2 += npx!(j + 1);
                    if k2 <= npx!(j) {
                        break;
                    }
                }
                j = 1;
                label = Label::L390;
            }
            Label::L390 => {
                let mut to_l370 = false;
                loop {
                    if kk < k2 {
                        to_l370 = true;
                        break;
                    }
                    kk += inc;
                    k2 += kspan;
                    if k2 >= ks {
                        break;
                    }
                }
                if to_l370 {
                    label = Label::L370;
                    continue;
                }
                if kk < ks {
                    label = Label::L380;
                    continue;
                }
                jc = k3;
                label = Label::L440;
            }
            Label::L440 => {
                if 2 * kt + 1 >= m {
                    label = Label::Return;
                    continue;
                }
                kspnn = npx!(kt + 1);
                // permutation for square-free factors of n; nfac becomes cumulative products
                nn = m - kt;
                nf!(nn) = 1;
                j = nn;
                while j > kt {
                    nf!(j - 1) *= nf!(j);
                    j -= 1;
                }
                kt += 1;
                nn = nf!(kt - 1) - 1;
                jj = 0;
                j = 0;
                label = Label::L480;
            }
            Label::L460 => {
                jj -= k2;
                k2 = kk;
                k += 1;
                kk = nf!(k - 1);
                label = Label::L470;
            }
            Label::L470 => {
                jj += kk;
                if jj >= k2 {
                    label = Label::L460;
                    continue;
                }
                npx!(j) = jj;
                label = Label::L480;
            }
            Label::L480 => {
                k2 = nf!(kt - 1);
                k = kt + 1;
                kk = nf!(k - 1);
                j += 1;
                if j <= nn {
                    label = Label::L470;
                    continue;
                }
                // determine the permutation cycles of length greater than 1
                j = 0;
                label = Label::L500;
            }
            Label::L490 => {
                loop {
                    k = kk;
                    kk = npx!(k);
                    npx!(k) = -kk;
                    if kk == j {
                        break;
                    }
                }
                k3 = kk;
                label = Label::L500;
            }
            Label::L500 => {
                loop {
                    j += 1;
                    kk = npx!(j);
                    if kk >= 0 {
                        break;
                    }
                }
                if kk != j {
                    label = Label::L490;
                    continue;
                }
                npx!(j) = -j;
                if j != nn {
                    label = Label::L500;
                    continue;
                }
                maxf *= inc;
                label = Label::L570;
            }
            Label::Ord => {
                // reorder a and b, following the permutation cycles
                loop {
                    j -= 1;
                    if npx!(j) >= 0 {
                        break;
                    }
                }
                jj = jc;
                label = Label::L520;
            }
            Label::L520 => {
                kspan = jj.min(maxf);
                jj -= kspan;
                k = npx!(j);
                kk = jc * k + i + jj;
                k1 = kk + kspan;
                k2 = 1;
                while k1 != kk {
                    at[k2 as usize] = a!(k1);
                    bt[k2 as usize] = b!(k1);
                    k1 -= inc;
                    k2 += 1;
                }
                loop {
                    k1 = kk + kspan;
                    k2 = k1 - jc * (k + npx!(k));
                    k = -npx!(k);
                    loop {
                        a!(k1) = a!(k2);
                        b!(k1) = b!(k2);
                        k1 -= inc;
                        k2 -= inc;
                        if k1 == kk {
                            break;
                        }
                    }
                    kk = k2;
                    if k == j {
                        break;
                    }
                }
                k1 = kk + kspan;
                k2 = 1;
                while k1 > kk {
                    a!(k1) = at[k2 as usize];
                    b!(k1) = bt[k2 as usize];
                    k1 -= inc;
                    k2 += 1;
                }
                if jj != 0 {
                    label = Label::L520;
                    continue;
                }
                label = if j != 1 { Label::Ord } else { Label::L570 };
            }
            Label::L570 => {
                j = k3 + 1;
                nt -= kspnn;
                i = nt - inc + 1;
                label = if nt >= 0 { Label::Ord } else { Label::Return };
            }
            Label::Return => return,
        }
    }
}

/// `M_SQRT_3` from Rmath.h.
#[allow(clippy::excessive_precision)]
const M_SQRT_3: f64 = 1.732050807568877293527446341505872367;

/// `fft(z, inverse)` for a complex vector: the unnormalised discrete Fourier transform,
/// `sum_k z[k] exp(-2 pi i j k / n)` (forward) or with `+i` (inverse), bit-identical to R.
/// Vectors of length 0 or 1 come back unchanged, as in R.
pub fn fft(z: &[Complex], inverse: bool) -> Result<Vec<Complex>, FftError> {
    let n = z.len();
    if n <= 1 {
        return Ok(z.to_vec());
    }
    let n32 = i32::try_from(n).map_err(|_| FftError)?;
    let fac = fft_factor(n32)?;
    let mut buf: Vec<f64> = Vec::with_capacity(2 * n);
    for c in z {
        buf.push(c.re);
        buf.push(c.im);
    }
    fftmx(&mut buf, n32, if inverse { 2 } else { -2 }, &fac);
    Ok(buf.as_chunks::<2>().0.iter().map(|&[re, im]| Complex::new(re, im)).collect())
}

/// `fft(x, inverse)` of a real vector (R coerces it to complex with zero imaginary parts).
pub fn fft_real(x: &[f64], inverse: bool) -> Result<Vec<Complex>, FftError> {
    let z: Vec<Complex> = x.iter().map(|&re| Complex::new(re, 0.0)).collect();
    fft(&z, inverse)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The direct O(n^2) transform, for a tolerance check of the port's algebra (the R vectors
    /// in tests/fixtures/rcompat check the bits).
    fn dft(z: &[Complex], inverse: bool) -> Vec<Complex> {
        let n = z.len();
        let sign = if inverse { 1.0 } else { -1.0 };
        (0..n)
            .map(|j| {
                let mut acc = Complex::default();
                for (k, c) in z.iter().enumerate() {
                    let t = sign * 2.0 * std::f64::consts::PI * ((j * k) % n) as f64 / n as f64;
                    let w = Complex::new(t.cos(), t.sin());
                    let p = *c * w;
                    acc.re += p.re;
                    acc.im += p.im;
                }
                acc
            })
            .collect()
    }

    #[test]
    fn factorizations() {
        let f = fft_factor(1024).unwrap();
        assert_eq!(&f.nfac[..f.m_fac as usize], &[4, 4, 4, 4, 4]);
        assert_eq!((f.kt, f.maxf), (2, 4));
        let f = fft_factor(512).unwrap();
        assert_eq!(&f.nfac[..f.m_fac as usize], &[4, 4, 2, 4, 4]);
        let f = fft_factor(2 * 3 * 3 * 5 * 7 * 11).unwrap();
        assert_eq!(&f.nfac[..f.m_fac as usize], &[3, 2, 5, 7, 11, 3]);
    }

    #[test]
    fn matches_the_direct_transform() {
        let mut state = 12345u64;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        for n in [2usize, 3, 4, 5, 6, 7, 8, 9, 12, 16, 25, 27, 30, 32, 49, 60, 64, 97, 100, 121, 210, 256, 360, 512, 1000, 1024] {
            let z: Vec<Complex> = (0..n).map(|_| Complex::new(next(), next())).collect();
            for inverse in [false, true] {
                let got = fft(&z, inverse).unwrap();
                let want = dft(&z, inverse);
                for (g, w) in got.iter().zip(&want) {
                    assert!((g.re - w.re).abs() < 1e-9 && (g.im - w.im).abs() < 1e-9,
                            "n = {n}, inverse = {inverse}: {g:?} vs {w:?}");
                }
            }
        }
    }
}
