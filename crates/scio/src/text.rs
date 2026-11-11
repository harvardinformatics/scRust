//! Numbers as text, exactly as the R reference scripts write them (`fmt17` in
//! `validation/R/lib/io.R`): C's `printf("%.17g")`, which glibc rounds correctly, so 17
//! significant digits identify every double and Rust's correctly rounded parser gives back the
//! same bits. Integer-valued doubles come out as integers; trailing zeros are dropped.

/// R's `NA_real_` test (a NaN whose low word is 1954).
pub fn is_na(x: f64) -> bool {
    x.is_nan() && (x.to_bits() & 0xffff_ffff) == 1954
}

/// C's `%.{p}g` for a finite double, p >= 1: scientific notation when the decimal exponent X of
/// the value rounded to p digits is below -4 or at least p, fixed otherwise; trailing zeros
/// (and a bare decimal point) removed; the exponent has at least two digits.
pub fn fmt_g(x: f64, p: usize) -> String {
    debug_assert!(x.is_finite() && p >= 1);
    // the %e conversion decides the exponent (Rust's float formatting is exact, as glibc's)
    let e = format!("{:.*e}", p - 1, x);
    let (mant, exp) = e.split_once('e').unwrap();
    let xexp: i32 = exp.parse().unwrap();
    let s = if xexp < -4 || xexp >= p as i32 {
        let mant = strip_zeros(mant);
        format!("{mant}e{}{:02}", if xexp < 0 { '-' } else { '+' }, xexp.abs())
    } else {
        strip_zeros(&format!("{:.*}", (p as i32 - 1 - xexp) as usize, x)).to_string()
    };
    if x == 0.0 && x.is_sign_negative() && !s.starts_with('-') {
        format!("-{s}")
    } else {
        s
    }
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// `fmt17(x)`: `%.17g`, with R's spellings `NA`, `NaN`, `Inf` and `-Inf`.
pub fn fmt17(x: f64) -> String {
    if is_na(x) {
        "NA".into()
    } else if x.is_nan() {
        "NaN".into()
    } else if x == f64::INFINITY {
        "Inf".into()
    } else if x == f64::NEG_INFINITY {
        "-Inf".into()
    } else {
        fmt_g(x, 17)
    }
}

/// Parses what [`fmt17`] writes back into the exact double (`NA` becomes R's `NA_real_`).
pub fn parse17(s: &str) -> Option<f64> {
    match s {
        "NA" => Some(f64::from_bits(0x7ff0_0000_0000_07a2)),
        "NaN" => Some(f64::from_bits(0x7ff8_0000_0000_0000)),
        "Inf" => Some(f64::INFINITY),
        "-Inf" => Some(f64::NEG_INFINITY),
        _ => s.parse().ok(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn familiar_values() {
        for (x, s) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (1.0, "1"),
            (24.0, "24"),
            (0.25, "0.25"),
            (0.29, "0.28999999999999998"),
            (0.2413793103448276, "0.2413793103448276"),
            (1e-5, "1.0000000000000001e-05"),
            (0.0001, "0.0001"),
            (1e17, "1e+17"),
            (123456789012345678.0, "1.2345678901234568e+17"),
            (12345678901234567.0, "12345678901234568"),
            (202.115_264_261_982_8, "202.11526426198279"),
            (f64::INFINITY, "Inf"),
            (f64::from_bits(0x7ff0_0000_0000_07a2), "NA"),
        ] {
            assert_eq!(fmt17(x), s, "{x:e}");
            assert_eq!(parse17(s).unwrap().to_bits(), x.to_bits(), "{s}");
        }
    }

    /// `%.17g` against glibc's printf, which R calls.
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    #[test]
    fn matches_glibc_printf() {
        use std::ffi::{c_char, c_int, CStr, CString};
        extern "C" {
            fn snprintf(buf: *mut c_char, n: usize, fmt: *const c_char, ...) -> c_int;
        }
        let fmt = CString::new("%.*g").unwrap();
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for i in 0..300_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let x = match i % 3 {
                0 => f64::from_bits(state),
                1 => (state >> 11) as f64 / (1u64 << 53) as f64,
                _ => ((state % 100_000) as f64) / [1.0, 7.0, 64.0, 1000.0][(i / 3) % 4],
            };
            if !x.is_finite() {
                continue;
            }
            let p = [17, 15, 1, 6, 20][i % 5] as c_int;
            let mut buf = vec![0 as c_char; 512];
            // SAFETY: the buffer is large enough for %.20g and its size is passed.
            let c = unsafe {
                snprintf(buf.as_mut_ptr(), buf.len(), fmt.as_ptr(), p, x);
                CStr::from_ptr(buf.as_ptr()).to_string_lossy().into_owned()
            };
            assert_eq!(fmt_g(x, p as usize), c, "%.{p}g of {x:e}");
        }
    }
}
