//! R's decimal parser and `as.numeric()` of strings: `R_strtod5` in `src/main/util.c` and
//! `RealFromString` in `src/main/coerce.c`, R 4.5.3.
//!
//! R's parser is not correctly rounded. It accumulates the digits in an x87 long double
//! (`ans = 10*ans + digit`, each step rounded to a 64-bit significand), scales by a power of ten
//! built by repeated squaring in long double, divides or multiplies once, and finally rounds the
//! long double to a double. The two roundings can leave the result one unit in the last place
//! away from the correctly rounded value, so it is emulated here operation by operation with
//! [`F80`]. DoubletFinder passes pK through this parser (`summarizeSweep` and
//! `as.numeric(as.character(...))`, plan B1), and the parsed value sets k and the column names.

use crate::ldouble::F80;
use crate::summary::NA_REAL;

/// R's `NaN` (`R_NaN`).
const R_NAN: f64 = f64::from_bits(0x7ff8_0000_0000_0000);
/// `MAX_EXPONENT_PREFIX`: an exponent stops growing once it reaches this many.
const MAX_EXPONENT_PREFIX: i32 = 9999;

/// C's `isspace` in the C locale.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn starts_with_ignore_case(s: &[u8], prefix: &[u8]) -> bool {
    s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

fn ld(x: f64) -> F80 {
    F80::from_f64(x)
}

fn ld_int(n: u64) -> F80 {
    F80::from_u64(n)
}

/// `R_strtod5(str, &endptr, dec, NA, exact = FALSE)`: the value and the number of bytes
/// consumed (`endptr - str`; 0 when nothing could be parsed, in which case the value is `NA`).
///
/// `na` accepts the string "NA" (as `scan()` does); `as.numeric()` passes false.
pub fn r_strtod5(s: &[u8], dec: u8, na: bool) -> (f64, usize) {
    let at = |i: usize| -> u8 { s.get(i).copied().unwrap_or(0) };
    let mut ans = F80::ZERO;
    let mut sign = 1.0f64;
    let mut p = 0usize;

    // optional whitespace
    while is_space(at(p)) {
        p += 1;
    }
    if na && s[p..].starts_with(b"NA") {
        return done(F80::from_f64(NA_REAL), sign, p + 2);
    }
    // optional sign
    match at(p) {
        b'-' => {
            sign = -1.0;
            p += 1;
        }
        b'+' => p += 1,
        _ => {}
    }
    if starts_with_ignore_case(&s[p..], b"nan") {
        return done(ld(R_NAN), sign, p + 3);
    } else if starts_with_ignore_case(&s[p..], b"infinity") {
        return done(F80::INFINITY, sign, p + 8);
    } else if starts_with_ignore_case(&s[p..], b"inf") {
        return done(F80::INFINITY, sign, p + 3);
    }

    let mut expn: i32 = 0;
    // strlen(p) > 2: the hexadecimal form needs at least one character after "0x"
    let rest = &s[p..];
    let len = rest.iter().position(|&c| c == 0).unwrap_or(rest.len());
    if len > 2 && rest[0] == b'0' && (rest[1] == b'x' || rest[1] == b'X') {
        let sixteen = ld_int(16);
        let mut exph: i32 = -1;
        p += 2;
        loop {
            let c = at(p);
            let digit = match c {
                b'0'..=b'9' => Some(c - b'0'),
                b'a'..=b'f' => Some(c - b'a' + 10),
                b'A'..=b'F' => Some(c - b'A' + 10),
                _ => None,
            };
            if let Some(d) = digit {
                ans = sixteen * ans + ld_int(d as u64);
            } else if c == dec && c != 0 {
                exph = 0;
                p += 1;
                continue;
            } else {
                break;
            }
            if exph >= 0 {
                exph += 4;
            }
            p += 1;
        }
        // binary exponent, if any
        if at(p) == b'p' || at(p) == b'P' {
            let mut expsign = 1;
            p += 1;
            match at(p) {
                b'-' => {
                    expsign = -1;
                    p += 1;
                }
                b'+' => p += 1,
                _ => {}
            }
            let (mut n, mut ndig) = (0i32, 0);
            while at(p).is_ascii_digit() {
                n = if n < MAX_EXPONENT_PREFIX { n * 10 + (at(p) - b'0') as i32 } else { n };
                p += 1;
                ndig += 1;
            }
            if ndig == 0 {
                return done(F80::from_f64(NA_REAL), sign, 0);
            }
            expn += expsign * n;
        }
        if ans != F80::ZERO {
            // fac is a long double; p2 a double, squared in double
            let mut fac;
            let mut p2 = 2.0f64;
            if exph > 0 {
                if expn - exph < -122 {
                    let mut n = exph;
                    fac = ld(1.0);
                    while n != 0 {
                        if n & 1 != 0 {
                            fac = fac * ld(p2);
                        }
                        n >>= 1;
                        p2 *= p2;
                    }
                    ans = ans / fac;
                    p2 = 2.0;
                } else {
                    expn -= exph;
                }
            }
            let mut n = expn.unsigned_abs();
            fac = ld(1.0);
            while n != 0 {
                if n & 1 != 0 {
                    fac = fac * ld(p2);
                }
                n >>= 1;
                p2 *= p2;
            }
            ans = if expn < 0 { ans / fac } else { ans * fac };
        }
        return done(ans, sign, p);
    }

    let ten = ld_int(10);
    let mut ndigits: i32 = 0;
    while at(p).is_ascii_digit() {
        ans = ten * ans + ld_int((at(p) - b'0') as u64);
        p += 1;
        ndigits += 1;
    }
    if at(p) == dec && dec != 0 {
        p += 1;
        while at(p).is_ascii_digit() {
            ans = ten * ans + ld_int((at(p) - b'0') as u64);
            p += 1;
            ndigits += 1;
            expn -= 1;
        }
    }
    if ndigits == 0 {
        // back out
        return done(F80::from_f64(NA_REAL), sign, 0);
    }

    if at(p) == b'e' || at(p) == b'E' {
        let mut expsign = 1;
        p += 1;
        match at(p) {
            b'-' => {
                expsign = -1;
                p += 1;
            }
            b'+' => p += 1,
            _ => {}
        }
        let (mut n, mut ndig) = (0i32, 0);
        while at(p).is_ascii_digit() {
            n = if n < MAX_EXPONENT_PREFIX { n * 10 + (at(p) - b'0') as i32 } else { n };
            p += 1;
            ndig += 1;
        }
        if ndig == 0 {
            return done(F80::from_f64(NA_REAL), sign, 0);
        }
        expn += expsign * n;
    }

    // avoid unnecessary underflow for large negative exponents
    if expn + ndigits < -300 {
        for _ in 0..ndigits {
            ans = ans / ld(10.0);
        }
        expn += ndigits;
    }
    let mut p10 = ld(10.0);
    let mut fac = ld(1.0);
    if expn < -307 {
        // use underflow, not overflow
        let mut n = -expn;
        while n != 0 {
            if n & 1 != 0 {
                fac = fac / p10;
            }
            n >>= 1;
            p10 = p10 * p10;
        }
        ans = ans * fac;
    } else if expn < 0 {
        // positive powers are exact
        let mut n = -expn;
        while n != 0 {
            if n & 1 != 0 {
                fac = fac * p10;
            }
            n >>= 1;
            p10 = p10 * p10;
        }
        ans = ans / fac;
    } else if ans != F80::ZERO {
        let mut n = expn;
        while n != 0 {
            if n & 1 != 0 {
                fac = fac * p10;
            }
            n >>= 1;
            p10 = p10 * p10;
        }
        ans = ans * fac;
    }

    // explicit overflow to infinity
    if ans > ld(f64::MAX) {
        return (if sign > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY }, p);
    }
    done(ans, sign, p)
}

/// `done:` in `R_strtod5`: `sign * (double) ans`. A NaN keeps its bits whatever the sign (the
/// product of -1 and a NaN is that NaN on x86-64); `NA` comes back quietened, as loading it into
/// the FPU leaves it, which R still reads as `NA` (its low word is kept).
fn done(ans: F80, sign: f64, end: usize) -> (f64, usize) {
    let d = ans.to_f64();
    if d.is_nan() {
        return (d, end);
    }
    (if sign < 0.0 { -d } else { d }, end)
}

/// `R_strtod(str, &endptr)`: [`r_strtod5`] with `.` as the decimal point and no "NA".
pub fn r_strtod(s: &str) -> (f64, usize) {
    r_strtod5(s.as_bytes(), b'.', false)
}

/// `isBlankString` for ASCII text: true if every byte is C whitespace.
fn is_blank(s: &[u8]) -> bool {
    s.iter().all(|&c| is_space(c))
}

/// `as.numeric(x)` for one string (`RealFromString`): `NA` for a blank string or when anything
/// but whitespace follows the number (R also warns, "NAs introduced by coercion").
pub fn as_numeric(s: &str) -> f64 {
    let b = s.as_bytes();
    if is_blank(b) {
        return NA_REAL;
    }
    let (x, end) = r_strtod5(b, b'.', false);
    if is_blank(&b[end..]) {
        x
    } else {
        NA_REAL
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::summary::is_na;

    #[test]
    fn plain_numbers() {
        assert_eq!(as_numeric("0.25"), 0.25);
        assert_eq!(as_numeric(" 12 "), 12.0);
        assert_eq!(as_numeric("-1.5e3"), -1500.0);
        assert_eq!(as_numeric("1e-2"), 0.01);
        assert_eq!(as_numeric("5e-04"), 0.0005);
        assert!(as_numeric("-0").is_sign_negative());
        assert_eq!(as_numeric("0x1p3"), 8.0);
        assert_eq!(as_numeric("0x1.8"), 1.5);
        assert_eq!(as_numeric("Inf"), f64::INFINITY);
        assert_eq!(as_numeric("-infinity"), f64::NEG_INFINITY);
        assert_eq!(as_numeric("1e999"), f64::INFINITY);
        assert!(as_numeric("NaN").is_nan() && !is_na(as_numeric("NaN")));
        for bad in ["", "  ", "abc", "1.2.3", "1e", "e5", ".", "-", "0x", "NA", "1 2"] {
            assert!(is_na(as_numeric(bad)), "{bad:?}");
        }
    }

    #[test]
    fn not_correctly_rounded() {
        // "0.29": 29 / 100 in long double, then rounded to double. Correctly rounded here, but
        // other strings come out one ulp off (the R vectors in tests/fixtures/rcompat check).
        assert_eq!(as_numeric("0.29"), 0.29);
        assert_eq!(r_strtod("0.29x"), (0.29, 4));
        // an exponent without digits backs out of the whole number
        let (x, end) = r_strtod("0.29e");
        assert!(is_na(x) && end == 0);
    }
}
