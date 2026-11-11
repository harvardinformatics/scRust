//! Doubles to strings as R's `as.character()` and `paste()` write them (`StringFromReal` in
//! `src/main/printutils.c`, `formatReal` and `scientific` in `src/main/format.c`, R 4.5.3).
//!
//! R sets 15 significant digits (`DBL_DIG`) and, per element: finds the number of significant
//! digits needed (at most 15) by scaling |x| by a power of ten in long double and rounding to a
//! 15-digit integer; chooses fixed notation if it is no wider than scientific (plus
//! `options(scipen)`, 0 by default); prints with C's `%.*f` or `%.*e`; and drops trailing zeros.
//! So `5e-04`, `0.001`, `1e+05`, `123456` and `0.3` (for 0.1 + 0.2).
//!
//! The digits printed are the exact decimal expansion of the double rounded half to even, as
//! glibc's printf produces them; Rust's float formatting is exact in the same way.

use crate::ldouble::F80;
use crate::libm;
use crate::powl10::{POWL10, POWL10_KMIN};

/// `DBL_DIG`: the precision `as.character` uses.
pub const DBL_DIG: i32 = 15;
/// Exponents up to which 10^k is exact in the 64-bit significand (`KP_MAX` in format.c).
const KP_MAX: i32 = 27;
/// `R_print.na_width`: the width of "NA".
const NA_WIDTH: i32 = 2;

/// `tbl[k]`: 10^k, exact in long double for k <= 27 (5^27 < 2^64).
fn tbl(k: i32) -> F80 {
    debug_assert!((0..=KP_MAX).contains(&k));
    let mut p = F80::from_u64(1);
    let ten = F80::from_u64(10);
    for _ in 0..k {
        p = p * ten; // exact
    }
    p
}

/// glibc's `powl(10, k)`.
fn powl10(k: i32) -> F80 {
    let (se, m) = POWL10[(k - POWL10_KMIN) as usize];
    F80::from_bits(se, m)
}

/// The result of `scientific()` for one number.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sci {
    pub neg: bool,
    /// Exponent of 10: |x| = alpha * 10^kpower with 1 <= alpha < 10.
    pub kpower: i32,
    /// Significant digits needed, at most `digits`.
    pub nsig: i32,
    /// Rounding to `digits` makes the number one digit wider (e.g. 9996 to 1e+04).
    pub roundingwidens: bool,
}

/// `scientific(x)` with `R_print.digits = digits`.
pub fn scientific(x: f64, digits: i32) -> Sci {
    if x == 0.0 {
        return Sci { neg: false, kpower: 0, nsig: 1, roundingwidens: false };
    }
    let neg = x < 0.0;
    let r = x.abs();
    if digits > DBL_DIG {
        // format_via_sprintf: "%#.*e" with digits - 1 decimals
        let s = format!("{:.*e}", (digits - 1) as usize, r);
        let (mant, exp) = s.split_once('e').unwrap();
        let mant = mant.as_bytes();
        let mut i = digits as usize;
        while i >= 2 && mant[i] == b'0' {
            i -= 1;
        }
        return Sci { neg, kpower: exp.parse().unwrap(), nsig: i as i32, roundingwidens: false };
    }
    let mut kp = libm::log10(r).floor() as i32 - digits + 1;
    let mut r_prec = F80::from_f64(r);
    if kp.abs() <= KP_MAX {
        if kp > 0 {
            r_prec = r_prec / tbl(kp);
        } else if kp < 0 {
            r_prec = r_prec * tbl(-kp);
        }
    } else {
        r_prec = r_prec / powl10(kp);
    }
    if r_prec < tbl(digits - 1) {
        r_prec = r_prec * F80::from_u64(10);
        kp -= 1;
    }
    // alpha: |x| rounded to a `digits`-digit integer
    let mut alpha = r_prec.round_ties_even().to_f64();
    let mut nsig = digits;
    for _ in 1..=digits {
        alpha /= 10.0;
        if alpha == alpha.floor() {
            nsig -= 1;
        } else {
            break;
        }
    }
    if nsig == 0 && digits > 0 {
        nsig = 1;
        kp += 1;
    }
    let kpower = kp + digits - 1;
    // Scientific format may round more than fixed: e.g. 9996 with 3 digits is 1e+04 but 9996.
    let rgt = (digits - kpower).clamp(0, KP_MAX);
    let fuzz = 0.5 / tbl(rgt).to_f64();
    let roundingwidens =
        kpower > 0 && kpower <= KP_MAX && F80::from_f64(r) < tbl(kpower) - F80::from_f64(fuzz);
    Sci { neg, kpower, nsig, roundingwidens }
}

/// The format `formatReal` chooses for a vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealFormat {
    /// Field width.
    pub w: i32,
    /// Digits after the decimal point (fixed) or in the mantissa (scientific).
    pub d: i32,
    /// 0 for fixed notation; else scientific, with `e + 1` exponent digits.
    pub e: i32,
}

/// `formatReal(x, n, &w, &d, &e, nsmall)` with `R_print.digits = digits` and
/// `R_print.scipen = scipen`: one common format for all of `x`.
pub fn format_real(x: &[f64], digits: i32, scipen: i32, nsmall: i32) -> RealFormat {
    let (mut naflag, mut nanflag, mut posinf, mut neginf) = (false, false, false, false);
    let mut neg = 0;
    let mut mnl = i32::MAX;
    let (mut mxl, mut rgt, mut mxsl, mut mxns) = (i32::MIN, i32::MIN, i32::MIN, i32::MIN);
    for &xi in x {
        if !xi.is_finite() {
            if crate::summary::is_na(xi) {
                naflag = true;
            } else if xi.is_nan() {
                nanflag = true;
            } else if xi > 0.0 {
                posinf = true;
            } else {
                neginf = true;
            }
        } else {
            let s = scientific(xi, digits);
            let mut left = s.kpower + 1;
            if s.roundingwidens {
                left -= 1;
            }
            let sleft = s.neg as i32 + if left <= 0 { 1 } else { left };
            let right = s.nsig - left;
            if s.neg {
                neg = 1;
            }
            rgt = rgt.max(right);
            mxl = mxl.max(left);
            mnl = mnl.min(left);
            mxsl = mxsl.max(sleft);
            mxns = mxns.max(s.nsig);
        }
    }
    if digits == 0 {
        rgt = 0;
    }
    if mxl < 0 {
        mxsl = 1 + neg;
    }
    if rgt < 0 {
        rgt = 0;
    }
    let mut w_f = mxsl + rgt + (rgt != 0) as i32;
    let mut e = if mxl > 100 || mnl <= -99 { 2 } else { 1 };
    let (mut w, mut d);
    if mxns != i32::MIN {
        d = mxns - 1;
        w = neg + (d > 0) as i32 + d + 4 + e;
        if w_f <= w + scipen {
            e = 0;
            if nsmall > rgt {
                rgt = nsmall;
                w_f = mxsl + rgt + (rgt != 0) as i32;
            }
            d = rgt;
            w = w_f;
        }
    } else {
        w = 0;
        d = 0;
        e = 0;
    }
    if naflag && w < NA_WIDTH {
        w = NA_WIDTH;
    }
    if (nanflag || posinf) && w < 3 {
        w = 3;
    }
    if neginf && w < 4 {
        w = 4;
    }
    RealFormat { w, d, e }
}

/// C's `%.*e`: a mantissa with `d` decimals and an exponent of at least two digits.
fn c_sci(x: f64, d: usize) -> String {
    let s = format!("{:.*e}", d, x);
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    format!("{}e{}{:02}", mant, if exp < 0 { '-' } else { '+' }, exp.abs())
}

/// `EncodeRealDrop0(x, w, d, e, ".")`: prints with `%w.de` or `%w.df` (right-justified in `w`
/// characters), then drops the trailing zeros of the decimal part (and a bare decimal point).
pub fn encode_real_drop0(x: f64, f: RealFormat) -> String {
    let w = f.w.max(0) as usize;
    let x = if x == 0.0 { 0.0 } else { x }; // no "-0"
    let mut s = if !x.is_finite() {
        let t = if crate::summary::is_na(x) {
            "NA"
        } else if x.is_nan() {
            "NaN"
        } else if x > 0.0 {
            "Inf"
        } else {
            "-Inf"
        };
        format!("{:>w$}", t)
    } else if f.e != 0 {
        format!("{:>w$}", c_sci(x, f.d as usize))
    } else {
        format!("{:>w$}", format!("{:.*}", f.d as usize, x))
    };
    if let Some(dot) = s.find('.') {
        let b = s.as_bytes();
        let mut end = dot + 1;
        let mut keep = dot; // drop everything from here if all decimals are zero
        while end < b.len() && b[end].is_ascii_digit() {
            if b[end] != b'0' {
                keep = end + 1;
            }
            end += 1;
        }
        s.replace_range(keep..end, "");
    }
    s
}

/// `as.character(x)` (and `paste(x)`) for one double: R's 15-significant-digit representation,
/// with `options(scipen = 0)`. `NA` gives "NA", as `paste()` writes it.
pub fn as_character(x: f64) -> String {
    as_character_scipen(x, 0)
}

/// `as.character(x)` under `options(scipen = scipen)`.
pub fn as_character_scipen(x: f64, scipen: i32) -> String {
    encode_real_drop0(x, format_real(&[x], DBL_DIG, scipen, 0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn familiar_strings() {
        for (x, s) in [
            (5e-4, "5e-04"),
            (0.001, "0.001"),
            (0.1 + 0.2, "0.3"),
            (1.0 / 3.0, "0.333333333333333"),
            (1e5, "1e+05"),
            (123456.0, "123456"),
            (1e15, "1e+15"),
            (0.25, "0.25"),
            (-0.0, "0"),
            (-1.5, "-1.5"),
            (1e-300, "1e-300"),
            (f64::NAN, "NaN"),
            (f64::NEG_INFINITY, "-Inf"),
            (crate::summary::NA_REAL, "NA"),
        ] {
            assert_eq!(as_character(x), s, "{x:e}");
        }
    }

    #[test]
    fn tables() {
        assert_eq!(tbl(27).to_bits(), (0x4058, 0xcecb_8f27_f420_0f3a));
        assert_eq!(powl10(0).to_f64(), 1.0);
        assert_eq!(powl10(-300).to_f64(), 1e-300);
    }

    /// Rust's exact float formatting against glibc's printf, which R calls.
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    #[test]
    fn rust_formatting_matches_glibc_printf() {
        use std::ffi::{c_char, c_int, CStr, CString};
        extern "C" {
            fn snprintf(buf: *mut c_char, n: usize, fmt: *const c_char, ...) -> c_int;
        }
        fn c_format(fmt: &str, d: c_int, x: f64) -> String {
            let mut buf = vec![0 as c_char; 512];
            let f = CString::new(fmt).unwrap();
            // SAFETY: the buffer is large enough for these formats and its size is passed.
            unsafe {
                snprintf(buf.as_mut_ptr(), buf.len(), f.as_ptr(), d, x);
                CStr::from_ptr(buf.as_ptr()).to_string_lossy().into_owned()
            }
        }
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for i in 0..200_000 {
            let x = match i % 3 {
                0 => f64::from_bits(next() & !(1 << 62)), // any magnitude below 2
                1 => f64::from_bits(next()),
                // exact ties: a few significant decimal digits ending in 5
                _ => ((next() % 2_000_000) as f64 + 0.5) / [1.0, 2.0, 4.0, 8.0, 16.0][(next() % 5) as usize],
            };
            if !x.is_finite() {
                continue;
            }
            let d = (next() % 18) as c_int;
            assert_eq!(c_sci(x, d as usize), c_format("%.*e", d, x), "%.{d}e of {x:e}");
            if x.abs() < 1e30 {
                assert_eq!(format!("{:.*}", d as usize, x), c_format("%.*f", d, x), "%.{d}f of {x:e}");
            }
        }
    }
}
