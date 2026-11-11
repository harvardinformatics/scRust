//! Software emulation of the x87 80-bit extended-precision `long double`, which R uses through
//! its `LDOUBLE` type on x86-64: in `sum()`, `mean()`, `var()` and number formatting.
//!
//! Every arithmetic operation rounds its exact result once to a 64-bit significand, ties to even,
//! as the x87 FPU does with its default precision control (extended) and rounding mode (nearest).
//! Conversion to `f64` rounds again, ties to even, like a store to a C `double`. Gradual underflow,
//! overflow to infinity, signed zeros, infinities and NaN propagation follow the FPU's rules, so
//! results are identical to R's on x86-64 while the emulation itself gives the same answer on any
//! machine (plan A5b). On x86-64 the tests check it against the FPU instruction by instruction.

use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};

const BIAS: i32 = 16383;
/// Exponent of the smallest normal number.
const EMIN: i32 = 1 - BIAS;
const EXP_MASK: u16 = 0x7fff;
const SIGN_BIT: u16 = 0x8000;
/// The explicit integer bit of the significand.
const INT_BIT: u64 = 1 << 63;
const QUIET_BIT: u64 = 1 << 62;

/// An x87 extended-precision value, stored as the FPU stores it: 16 bits of sign and biased
/// exponent, and a 64-bit significand with an explicit integer bit.
///
/// `==` and `<` compare numerically, like C's operators on `long double` (NaN is unordered,
/// `-0 == +0`); [`F80::to_bits`] gives the exact representation.
#[derive(Clone, Copy, Debug)]
pub struct F80 {
    se: u16,
    m: u64,
}

/// A decoded value. Finite non-zero values are normalised: `m` has its top bit set and the value
/// is `m * 2^(exp - 63)`, even for denormals.
enum Kind {
    Zero,
    Finite { exp: i32, m: u64 },
    Inf,
    /// A NaN, or an encoding the FPU rejects as an invalid operand (unnormal, pseudo-NaN,
    /// pseudo-infinity).
    Nan,
}

impl F80 {
    pub const ZERO: F80 = F80 { se: 0, m: 0 };
    pub const NEG_ZERO: F80 = F80 { se: SIGN_BIT, m: 0 };
    pub const INFINITY: F80 = F80 { se: EXP_MASK, m: INT_BIT };
    pub const NEG_INFINITY: F80 = F80 { se: SIGN_BIT | EXP_MASK, m: INT_BIT };
    /// The x87 "real indefinite" NaN, the result of invalid operations such as `Inf - Inf`.
    pub const INDEFINITE: F80 = F80 { se: SIGN_BIT | EXP_MASK, m: INT_BIT | QUIET_BIT };

    /// From the FPU's representation: sign and biased exponent, and the 64-bit significand.
    pub const fn from_bits(se: u16, m: u64) -> F80 {
        F80 { se, m }
    }

    /// The FPU's representation: sign and biased exponent, and the 64-bit significand.
    pub const fn to_bits(self) -> (u16, u64) {
        (self.se, self.m)
    }

    fn negative(self) -> bool {
        self.se & SIGN_BIT != 0
    }

    fn sign_bits(neg: bool) -> u16 {
        if neg {
            SIGN_BIT
        } else {
            0
        }
    }

    fn kind(self) -> Kind {
        let e = self.se & EXP_MASK;
        if e == EXP_MASK {
            if self.m == INT_BIT {
                Kind::Inf
            } else {
                // NaNs, and pseudo-infinities and pseudo-NaNs, which the FPU rejects as invalid
                Kind::Nan
            }
        } else if e == 0 {
            if self.m == 0 {
                Kind::Zero
            } else {
                // denormal (or pseudo-denormal): m * 2^(EMIN - 63)
                let lz = self.m.leading_zeros();
                Kind::Finite { exp: EMIN - lz as i32, m: self.m << lz }
            }
        } else if self.m & INT_BIT == 0 {
            Kind::Nan // unnormal: an invalid operand since the 80387
        } else {
            Kind::Finite { exp: e as i32 - BIAS, m: self.m }
        }
    }

    pub fn is_nan(self) -> bool {
        matches!(self.kind(), Kind::Nan)
    }

    /// A proper NaN: maximal exponent, integer bit set, non-zero fraction.
    fn is_proper_nan(self) -> bool {
        self.se & EXP_MASK == EXP_MASK && self.m & INT_BIT != 0 && self.m != INT_BIT
    }

    /// The same value in the FPU's canonical encoding (a pseudo-denormal becomes normal).
    fn canonical(self) -> F80 {
        match self.kind() {
            Kind::Finite { exp, m } => round_pack(self.negative(), exp - 63, m as u128),
            _ => self,
        }
    }

    /// True if neither infinite nor NaN.
    pub fn is_finite(self) -> bool {
        matches!(self.kind(), Kind::Zero | Kind::Finite { .. })
    }

    /// Exact conversion, as the FPU loads a `double` (`fld`); signalling NaNs become quiet.
    pub fn from_f64(x: f64) -> F80 {
        let bits = x.to_bits();
        let sign = Self::sign_bits(bits >> 63 != 0);
        let e = ((bits >> 52) & 0x7ff) as i32;
        let f = bits & ((1u64 << 52) - 1);
        if e == 0x7ff {
            if f == 0 {
                F80 { se: sign | EXP_MASK, m: INT_BIT }
            } else {
                F80 { se: sign | EXP_MASK, m: INT_BIT | QUIET_BIT | (f << 11) }
            }
        } else if e == 0 {
            if f == 0 {
                F80 { se: sign, m: 0 }
            } else {
                // a double denormal is a normal extended number: f * 2^-1074
                let lz = f.leading_zeros();
                let exp = 63 - lz as i32 - 1074;
                F80 { se: sign | (exp + BIAS) as u16, m: f << lz }
            }
        } else {
            F80 { se: sign | (e - 1023 + BIAS) as u16, m: INT_BIT | (f << 11) }
        }
    }

    /// `self + x` for a double operand in memory (`fadd m64fp`), the way R's loops add doubles
    /// to a long-double accumulator (`s += x[i]`). It equals `self + F80::from_f64(x)` except
    /// when both are NaN: the double keeps its signalling status, so a signalling NaN (R's `NA`
    /// is one) loses to a quiet NaN already in the accumulator.
    pub fn add_f64(self, x: f64) -> F80 {
        let bits = x.to_bits();
        if x.is_nan() && bits & (1 << 51) == 0 {
            let snan = F80 {
                se: Self::sign_bits(bits >> 63 != 0) | EXP_MASK,
                m: INT_BIT | ((bits & ((1u64 << 52) - 1)) << 11),
            };
            return self.propagate_nan(snan).unwrap();
        }
        self + F80::from_f64(x)
    }

    /// Exact conversion of an integer, as `(long double) n`.
    pub fn from_u64(n: u64) -> F80 {
        if n == 0 {
            F80::ZERO
        } else {
            round_pack(false, 0, n as u128)
        }
    }

    /// Rounds to `f64`, ties to even, as a store to a C `double` does (`fstp`): overflow gives
    /// infinity and tiny values become denormal or zero. NaN payloads are kept, so R's `NA`
    /// stays `NA`.
    pub fn to_f64(self) -> f64 {
        let sign = ((self.se >> 15) as u64) << 63;
        match self.kind() {
            Kind::Zero => f64::from_bits(sign),
            Kind::Inf => f64::from_bits(sign | 0x7ff0_0000_0000_0000),
            Kind::Nan => f64::from_bits(
                sign | 0x7ff0_0000_0000_0000 | (1 << 51) | ((self.m >> 11) & ((1 << 52) - 1)),
            ),
            Kind::Finite { exp, m } => {
                if exp > 1023 {
                    return f64::from_bits(sign | 0x7ff0_0000_0000_0000);
                }
                if exp >= -1022 {
                    let k = shift_right_round(m as u128, 11) as u64;
                    let (k, exp) = if k == 1 << 53 { (1u64 << 52, exp + 1) } else { (k, exp) };
                    if exp > 1023 {
                        return f64::from_bits(sign | 0x7ff0_0000_0000_0000);
                    }
                    f64::from_bits(sign | (((exp + 1023) as u64) << 52) | (k & ((1 << 52) - 1)))
                } else {
                    // denormal result k * 2^-1074; k == 2^52 is the smallest normal number
                    let shift = (11 + (-1022 - exp)) as u32;
                    f64::from_bits(sign | shift_right_round(m as u128, shift) as u64)
                }
            }
        }
    }

    /// Rounds to an integral value, ties to even, as C's `nearbyintl` does in the default
    /// rounding mode.
    pub fn round_ties_even(self) -> F80 {
        let (exp, m) = match self.kind() {
            Kind::Finite { exp, m } => (exp, m),
            _ => return self,
        };
        if exp >= 63 {
            return self; // already an integer
        }
        let neg = self.negative();
        if exp < -1 {
            return F80 { se: Self::sign_bits(neg), m: 0 }; // |x| < 0.5
        }
        let k = shift_right_round(m as u128, (63 - exp) as u32);
        if k == 0 {
            F80 { se: Self::sign_bits(neg), m: 0 }
        } else {
            round_pack(neg, 0, k)
        }
    }

    /// x87 NaN propagation for a two-operand instruction (Intel SDM vol. 1, table 4-7): a single
    /// NaN operand is returned, quietened; of a quiet and a signalling NaN, the quiet one; of two
    /// NaNs of the same kind, the one with the larger significand, and on a tie the positive one
    /// (as measured on the FPU). Invalid encodings give the default NaN.
    fn propagate_nan(self, rhs: F80) -> Option<F80> {
        if !self.is_nan() && !rhs.is_nan() {
            return None;
        }
        if (self.is_nan() && !self.is_proper_nan()) || (rhs.is_nan() && !rhs.is_proper_nan()) {
            return Some(F80::INDEFINITE);
        }
        let quiet = |x: F80| F80 { se: x.se, m: x.m | QUIET_BIT };
        Some(match (self.is_nan(), rhs.is_nan()) {
            (true, false) => quiet(self),
            (false, true) => quiet(rhs),
            _ => {
                let (qa, qb) = (self.m & QUIET_BIT != 0, rhs.m & QUIET_BIT != 0);
                if qa != qb {
                    if qa { self } else { rhs }
                } else if rhs.m > self.m || (rhs.m == self.m && !rhs.negative()) {
                    quiet(rhs)
                } else {
                    quiet(self)
                }
            }
        })
    }
}

/// `w >> s`, rounded to nearest with ties to even.
fn shift_right_round(w: u128, s: u32) -> u128 {
    if s == 0 {
        return w;
    }
    if s > 128 {
        return 0; // w < 2^128 <= half
    }
    let (kept, rem) = if s == 128 { (0, w) } else { (w >> s, w & ((1u128 << s) - 1)) };
    let half = 1u128 << (s - 1);
    if rem > half || (rem == half && kept & 1 == 1) {
        kept + 1
    } else {
        kept
    }
}

/// Rounds the exact value `w * 2^q` (w > 0) to extended precision, ties to even, with gradual
/// underflow and overflow to infinity. When `w` carries a sticky bit for discarded low-order bits,
/// that bit must lie at least two places below the rounding position, which every caller ensures.
fn round_pack(neg: bool, q: i32, w: u128) -> F80 {
    debug_assert!(w != 0);
    let sign = F80::sign_bits(neg);
    let p = 127 - w.leading_zeros() as i32; // position of the top bit
    let e = q + p; // exponent of the top bit
    // Bits to drop: down to 64 significant bits, and further for a denormal result.
    let shift = (p - 63) + (EMIN - e).max(0);
    let (mut m, mut q) = if shift <= 0 {
        (w << (-shift) as u32, q + shift)
    } else {
        (shift_right_round(w, shift as u32), q + shift)
    };
    if m == 1u128 << 64 {
        m >>= 1; // rounding carried into a new bit
        q += 1;
    }
    let m = m as u64;
    if m & INT_BIT != 0 {
        let field = q + 63 + BIAS;
        if field >= EXP_MASK as i32 {
            return F80 { se: sign | EXP_MASK, m: INT_BIT };
        }
        F80 { se: sign | field as u16, m }
    } else {
        // denormal or zero: m * 2^(EMIN - 63)
        F80 { se: sign, m }
    }
}

fn add_finite(neg_a: bool, ea: i32, ma: u64, neg_b: bool, eb: i32, mb: u64) -> F80 {
    // Order by magnitude, then align the smaller operand to the larger one, keeping 62 extra
    // low-order bits and a sticky bit for anything shifted out beyond them.
    let ((neg_a, ea, ma), (neg_b, eb, mb)) = if (ea, ma) >= (eb, mb) {
        ((neg_a, ea, ma), (neg_b, eb, mb))
    } else {
        ((neg_b, eb, mb), (neg_a, ea, ma))
    };
    let d = (ea - eb) as u32;
    let a = (ma as u128) << 62;
    let b_full = (mb as u128) << 62;
    let b = if d == 0 {
        b_full
    } else if d < 128 {
        (b_full >> d) | ((b_full & ((1u128 << d) - 1) != 0) as u128)
    } else {
        1
    };
    let s = if neg_a == neg_b { a + b } else { a - b };
    if s == 0 {
        return F80::ZERO; // exact cancellation gives +0 when rounding to nearest
    }
    round_pack(neg_a, ea - 125, s)
}

impl Add for F80 {
    type Output = F80;
    fn add(self, rhs: F80) -> F80 {
        if let Some(nan) = self.propagate_nan(rhs) {
            return nan;
        }
        match (self.kind(), rhs.kind()) {
            (Kind::Inf, Kind::Inf) => {
                if self.negative() == rhs.negative() {
                    self
                } else {
                    F80::INDEFINITE
                }
            }
            (Kind::Inf, _) => self,
            (_, Kind::Inf) => rhs,
            (Kind::Zero, Kind::Zero) => {
                if self.negative() && rhs.negative() {
                    F80::NEG_ZERO
                } else {
                    F80::ZERO
                }
            }
            (Kind::Zero, _) => rhs.canonical(),
            (_, Kind::Zero) => self.canonical(),
            (Kind::Finite { exp: ea, m: ma }, Kind::Finite { exp: eb, m: mb }) => {
                add_finite(self.negative(), ea, ma, rhs.negative(), eb, mb)
            }
            (Kind::Nan, _) | (_, Kind::Nan) => unreachable!(),
        }
    }
}

impl Neg for F80 {
    type Output = F80;
    fn neg(self) -> F80 {
        F80 { se: self.se ^ SIGN_BIT, m: self.m }
    }
}

impl Sub for F80 {
    type Output = F80;
    fn sub(self, rhs: F80) -> F80 {
        if let Some(nan) = self.propagate_nan(rhs) {
            return nan;
        }
        self + (-rhs)
    }
}

impl Mul for F80 {
    type Output = F80;
    fn mul(self, rhs: F80) -> F80 {
        if let Some(nan) = self.propagate_nan(rhs) {
            return nan;
        }
        let neg = self.negative() != rhs.negative();
        match (self.kind(), rhs.kind()) {
            (Kind::Inf, Kind::Zero) | (Kind::Zero, Kind::Inf) => F80::INDEFINITE,
            (Kind::Inf, _) | (_, Kind::Inf) => {
                F80 { se: F80::sign_bits(neg) | EXP_MASK, m: INT_BIT }
            }
            (Kind::Zero, _) | (_, Kind::Zero) => F80 { se: F80::sign_bits(neg), m: 0 },
            (Kind::Finite { exp: ea, m: ma }, Kind::Finite { exp: eb, m: mb }) => {
                round_pack(neg, ea + eb - 126, ma as u128 * mb as u128)
            }
            (Kind::Nan, _) | (_, Kind::Nan) => unreachable!(),
        }
    }
}

impl Div for F80 {
    type Output = F80;
    fn div(self, rhs: F80) -> F80 {
        if let Some(nan) = self.propagate_nan(rhs) {
            return nan;
        }
        let neg = self.negative() != rhs.negative();
        match (self.kind(), rhs.kind()) {
            (Kind::Inf, Kind::Inf) | (Kind::Zero, Kind::Zero) => F80::INDEFINITE,
            (Kind::Inf, _) | (_, Kind::Zero) => {
                F80 { se: F80::sign_bits(neg) | EXP_MASK, m: INT_BIT }
            }
            (_, Kind::Inf) | (Kind::Zero, _) => F80 { se: F80::sign_bits(neg), m: 0 },
            (Kind::Finite { exp: ea, m: ma }, Kind::Finite { exp: eb, m: mb }) => {
                // 66 quotient bits (at least 2 beyond the 64 kept), then a sticky bit
                let (ma, mb) = (ma as u128, mb as u128);
                let n = ma << 64;
                let (q, r) = (n / mb, n % mb);
                let q2 = (q << 2) | ((r << 2) / mb);
                let r2 = (r << 2) % mb;
                round_pack(neg, ea - eb - 66, q2 | (r2 != 0) as u128)
            }
            (Kind::Nan, _) | (_, Kind::Nan) => unreachable!(),
        }
    }
}

impl PartialEq for F80 {
    fn eq(&self, other: &F80) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

impl PartialOrd for F80 {
    fn partial_cmp(&self, other: &F80) -> Option<Ordering> {
        // (class, exponent, significand) orders magnitudes; zeros compare equal whatever the sign
        fn magnitude(x: F80) -> Option<(u8, i32, u64)> {
            match x.kind() {
                Kind::Zero => Some((0, 0, 0)),
                Kind::Finite { exp, m } => Some((1, exp, m)),
                Kind::Inf => Some((2, 0, 0)),
                Kind::Nan => None,
            }
        }
        let (ma, mb) = (magnitude(*self)?, magnitude(*other)?);
        let neg_a = self.negative() && ma.0 != 0;
        let neg_b = other.negative() && mb.0 != 0;
        Some(match (neg_a, neg_b) {
            (false, true) => Ordering::Greater,
            (true, false) => Ordering::Less,
            (false, false) => ma.cmp(&mb),
            (true, true) => mb.cmp(&ma),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_values() {
        assert_eq!(F80::from_f64(1.0).to_bits(), (0x3fff, INT_BIT));
        assert_eq!(F80::from_f64(-2.0).to_bits(), (0xc000, INT_BIT));
        assert_eq!(F80::from_u64(10).to_bits(), (0x4002, 0xa000_0000_0000_0000));
        assert_eq!(F80::from_f64(f64::from_bits(1)).to_f64().to_bits(), 1);
        assert_eq!((F80::from_f64(0.1) + F80::from_f64(0.2)).to_f64(), 0.30000000000000004);
        // (1 + 2^-60) is exact in extended precision but rounds to 1 as a double
        let x = F80::from_f64(1.0) + F80::from_f64(2f64.powi(-60));
        assert!(x > F80::from_f64(1.0));
        assert_eq!(x.to_f64(), 1.0);
        assert_eq!((F80::from_f64(1.0) / F80::from_u64(3)).to_bits(), (0x3ffd, 0xaaaa_aaaa_aaaa_aaab));
        assert_eq!(F80::from_f64(2.5).round_ties_even().to_f64(), 2.0);
        assert_eq!(F80::from_f64(-3.5).round_ties_even().to_f64(), -4.0);
        assert!(F80::from_f64(-0.25).round_ties_even().to_f64().is_sign_negative());
        // R's NA (a NaN with low word 1954) survives a round trip
        let na = f64::from_bits(0x7ff0_0000_0000_07a2);
        assert_eq!((F80::from_f64(na) + F80::from_f64(1.0)).to_f64().to_bits() & 0xffff_ffff, 1954);
        assert!(F80::INFINITY > F80::from_f64(f64::MAX));
        assert!(F80::from_f64(-0.0) == F80::ZERO);
    }

    #[test]
    fn overflow_and_underflow_of_doubles() {
        let big = F80::from_f64(f64::MAX) + F80::from_f64(f64::MAX);
        assert!(big.is_finite());
        assert_eq!(big.to_f64(), f64::INFINITY);
        assert_eq!((big / F80::from_u64(2)).to_f64(), f64::MAX);
        let tiny = F80::from_f64(f64::from_bits(1)) / F80::from_u64(3);
        assert!(tiny > F80::ZERO);
        assert_eq!(tiny.to_f64(), 0.0);
        let half_min = F80::from_f64(f64::from_bits(1)) / F80::from_u64(2);
        assert_eq!(half_min.to_f64(), 0.0); // tie: rounds to even (zero)
        let three_halves = F80::from_f64(f64::from_bits(3)) / F80::from_u64(2);
        assert_eq!(three_halves.to_f64().to_bits(), 2); // 1.5 ulp rounds to 2
    }

    /// The emulation against the x87 FPU itself, on random and adversarial operands.
    #[cfg(target_arch = "x86_64")]
    mod against_x87 {
        use super::super::*;
        use std::arch::asm;

        #[derive(Clone, Copy)]
        enum Op {
            Add,
            Sub,
            Mul,
            Div,
        }

        fn to_mem(x: F80) -> [u8; 16] {
            let mut b = [0u8; 16];
            b[..8].copy_from_slice(&x.m.to_le_bytes());
            b[8..10].copy_from_slice(&x.se.to_le_bytes());
            b
        }

        fn from_mem(b: &[u8; 16]) -> F80 {
            F80 {
                m: u64::from_le_bytes(b[..8].try_into().unwrap()),
                se: u16::from_le_bytes(b[8..10].try_into().unwrap()),
            }
        }

        fn control_word() -> u16 {
            let mut cw = 0u16;
            unsafe { asm!("fnstcw word ptr [{}]", in(reg) &mut cw as *mut u16, options(nostack)) };
            cw
        }

        /// `a op b` on the FPU: st(0) = a, st(1) = b, then `fop st(0), st(1)`.
        fn x87(op: Op, a: F80, b: F80) -> F80 {
            let (ma, mb, mut out) = (to_mem(a), to_mem(b), [0u8; 16]);
            unsafe {
                macro_rules! run {
                    ($inst:literal) => {
                        asm!(
                            "fld tbyte ptr [{b}]",
                            "fld tbyte ptr [{a}]",
                            $inst,
                            "fstp tbyte ptr [{o}]",
                            "fstp st(0)",
                            a = in(reg) ma.as_ptr(),
                            b = in(reg) mb.as_ptr(),
                            o = in(reg) out.as_mut_ptr(),
                            out("st(0)") _, out("st(1)") _,
                            options(nostack),
                        )
                    };
                }
                match op {
                    Op::Add => run!("fadd st(0), st(1)"),
                    Op::Sub => run!("fsub st(0), st(1)"),
                    Op::Mul => run!("fmul st(0), st(1)"),
                    Op::Div => run!("fdiv st(0), st(1)"),
                }
            }
            from_mem(&out)
        }

        /// `a + x` with x a double memory operand: st(0) = a, then `fadd qword ptr [x]`.
        fn x87_add_m64(a: F80, x: f64) -> F80 {
            let (ma, mut out) = (to_mem(a), [0u8; 16]);
            unsafe {
                asm!("fld tbyte ptr [{a}]", "fadd qword ptr [{x}]", "fstp tbyte ptr [{o}]",
                     a = in(reg) ma.as_ptr(), x = in(reg) &x as *const f64,
                     o = in(reg) out.as_mut_ptr(), out("st(0)") _, options(nostack));
            }
            from_mem(&out)
        }

        #[test]
        fn add_f64_matches_fadd_m64() {
            let mut g = SplitMix(1954);
            // R's NA (a signalling NaN), R's NaN, the default NaN, quietened NA, and others
            let nans = [0x7ff0_0000_0000_07a2u64, 0x7ff8_0000_0000_0000, 0xfff8_0000_0000_0000,
                        0x7ff8_0000_0000_07a2, 0x7ff0_0000_0000_0001, 0xfff4_0000_0000_0000];
            for i in 0..200_000 {
                let a = match i % 4 {
                    0 => F80::from_f64(f64::from_bits(nans[(g.next() % 6) as usize])),
                    _ => random_f80(&mut g, 16383, 70),
                };
                let x = match i % 3 {
                    0 => f64::from_bits(nans[(g.next() % 6) as usize]),
                    1 => f64::from_bits(g.next()),
                    _ => random_f80(&mut g, 16383, 70).to_f64(),
                };
                let (got, want) = (a.add_f64(x), x87_add_m64(a, x));
                assert_eq!(got.to_bits(), want.to_bits(), "{:?} + {:016x}", a.to_bits(), x.to_bits());
            }
        }

        fn x87_from_f64(x: f64) -> F80 {
            let mut out = [0u8; 16];
            unsafe {
                asm!("fld qword ptr [{src}]", "fstp tbyte ptr [{o}]",
                     src = in(reg) &x as *const f64, o = in(reg) out.as_mut_ptr(),
                     out("st(0)") _, options(nostack));
            }
            from_mem(&out)
        }

        fn x87_to_f64(a: F80) -> f64 {
            let (ma, mut out) = (to_mem(a), 0f64);
            unsafe {
                asm!("fld tbyte ptr [{a}]", "fstp qword ptr [{o}]",
                     a = in(reg) ma.as_ptr(), o = in(reg) &mut out as *mut f64,
                     out("st(0)") _, options(nostack));
            }
            out
        }

        fn x87_rint(a: F80) -> F80 {
            let (ma, mut out) = (to_mem(a), [0u8; 16]);
            unsafe {
                asm!("fld tbyte ptr [{a}]", "frndint", "fstp tbyte ptr [{o}]",
                     a = in(reg) ma.as_ptr(), o = in(reg) out.as_mut_ptr(),
                     out("st(0)") _, options(nostack));
            }
            from_mem(&out)
        }

        /// Bitwise equality, except that any two NaNs match only if both are NaN.
        fn same(a: F80, b: F80) -> bool {
            if a.is_nan() || b.is_nan() {
                a.is_nan() && b.is_nan()
            } else {
                a.to_bits() == b.to_bits()
            }
        }

        struct SplitMix(u64);
        impl SplitMix {
            fn next(&mut self) -> u64 {
                self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut z = self.0;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                z ^ (z >> 31)
            }
        }

        /// A random extended value: mostly normal numbers with exponents in a band around
        /// `centre`, some with short significands (to hit exact ties), and a few special values.
        fn random_f80(g: &mut SplitMix, centre: i32, spread: u32) -> F80 {
            let r = g.next();
            let sign = if r & 1 == 1 { SIGN_BIT } else { 0 };
            match (r >> 1) % 64 {
                0 => F80 { se: sign, m: 0 },
                1 => F80 { se: sign | EXP_MASK, m: INT_BIT },
                4 => F80 { se: sign | EXP_MASK, m: INT_BIT | QUIET_BIT | (g.next() % 3) << 11 },
                2 => F80 { se: sign, m: g.next() >> (g.next() % 64) }, // denormal
                3 => F80::from_f64(f64::from_bits(g.next())),
                _ => {
                    let e = centre + (g.next() % (2 * spread as u64 + 1)) as i32 - spread as i32;
                    let e = e.clamp(1, 0x7ffe);
                    let mut m = g.next() | INT_BIT;
                    if (r >> 8).is_multiple_of(4) {
                        m &= !((1u64 << (g.next() % 63)) - 1); // short significand
                    }
                    F80 { se: sign | e as u16, m }
                }
            }
        }

        #[test]
        fn fpu_is_in_default_mode() {
            // precision control = extended (bits 8-9 = 11), rounding = nearest (bits 10-11 = 00)
            assert_eq!(control_word() & 0x0f00, 0x0300);
        }

        #[test]
        fn arithmetic_matches_the_fpu() {
            let mut g = SplitMix(20261110);
            let ops = [Op::Add, Op::Sub, Op::Mul, Op::Div];
            // exponent bands: everyday values, near overflow, near and below underflow
            for &(centre, spread) in &[(16383, 3), (16383, 70), (16383, 2000), (32700, 80), (60, 80)] {
                for _ in 0..150_000 {
                    let a = random_f80(&mut g, centre, spread);
                    let b = random_f80(&mut g, centre, spread);
                    for &op in &ops {
                        let want = x87(op, a, b);
                        let got = match op {
                            Op::Add => a + b,
                            Op::Sub => a - b,
                            Op::Mul => a * b,
                            Op::Div => a / b,
                        };
                        assert!(same(got, want), "{:?} vs x87 {:?} for {:?} and {:?}",
                                got.to_bits(), want.to_bits(), a.to_bits(), b.to_bits());
                    }
                }
            }
        }

        #[test]
        fn sums_of_doubles_match_the_fpu() {
            // the operations R performs: long double accumulators plus doubles, and differences
            let mut g = SplitMix(7);
            for _ in 0..300_000 {
                let s = random_f80(&mut g, 16383, 60);
                let fx = F80::from_f64(random_f80(&mut g, 16383, 60).to_f64());
                assert!(same(s + fx, x87(Op::Add, s, fx)));
                assert!(same(fx - s, x87(Op::Sub, fx, s)));
            }
        }

        #[test]
        fn conversions_match_the_fpu() {
            let mut g = SplitMix(42);
            for i in 0..400_000 {
                let x = f64::from_bits(g.next());
                assert!(same(F80::from_f64(x), x87_from_f64(x)), "from_f64({x:e})");
                // values near the double range limits, rounding boundaries, and anywhere
                let a = match i % 4 {
                    0 => random_f80(&mut g, 16383, 1100),
                    1 => random_f80(&mut g, 16383 - 1022, 70),
                    2 => random_f80(&mut g, 16383 + 1023, 3),
                    _ => {
                        // significand exactly halfway between two doubles, or one off it
                        let m = (g.next() | INT_BIT) & !0x7ff | 0x400;
                        let m = m.wrapping_add(g.next() % 3).wrapping_sub(1) | INT_BIT;
                        F80 { se: (16383 + (g.next() % 2100) as i32 - 1100) as u16, m }
                    }
                };
                let (got, want) = (a.to_f64(), x87_to_f64(a));
                assert!(got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                        "to_f64({:?}): {:e} vs {:e}", a.to_bits(), got, want);
                assert!(same(a.round_ties_even(), x87_rint(a)), "rint({:?})", a.to_bits());
            }
        }
    }
}
