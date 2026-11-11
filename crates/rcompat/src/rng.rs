//! R's default random-number generator: Mersenne-Twister seeded by `set.seed`, uniforms with
//! `unif_rand`'s fixup, and `R_unif_index` for `sample.kind = "Rejection"`, ported from R 4.5.3's
//! `src/main/RNG.c`. Normal deviates (`norm_rand`, inversion) are in [`crate::normal`].

use crate::libm;
use std::fmt;

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const TEMPERING_MASK_B: u32 = 0x9d2c_5680;
const TEMPERING_MASK_C: u32 = 0xefc6_0000;
/// 1/(2^32 - 1), as written in RNG.c.
const I2_32M1: f64 = 2.328306437080797e-10;

/// `.Random.seed[1]` for R's defaults: Mersenne-Twister (3) + 100 * Inversion (4)
/// + 10000 * Rejection (1).
pub const RANDOM_SEED_KIND: i32 = 10403;
/// Length of `.Random.seed` for Mersenne-Twister: the kind code, `mti` and 624 state words.
pub const RANDOM_SEED_LEN: usize = 626;

/// R's random-number stream with its default kinds: `RNGkind("Mersenne-Twister", "Inversion",
/// "Rejection")`.
///
/// One generator is passed explicitly (`&mut RRng`) through every call, like R's single global
/// stream (plan A1). Every function consumes exactly as many uniforms as its R counterpart, so the
/// stream continues identically afterwards.
#[derive(Clone)]
pub struct RRng {
    /// The Mersenne-Twister state, `.Random.seed[3:626]`.
    mt: [u32; N],
    /// Position in the state, `.Random.seed[2]`; 624 means "generate the next block first".
    mti: u32,
}

/// A `.Random.seed` vector that R would reject or reinitialise from the clock.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedError {
    Length(usize),
    /// Only Mersenne-Twister with Inversion and Rejection (10403) is implemented.
    Kind(i32),
    /// An all-zero state: R would reseed from the time of day.
    AllZero,
    /// `mti` above 625 makes R read outside the state.
    Position(u32),
}

impl fmt::Display for SeedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SeedError::Length(n) => write!(f, ".Random.seed has length {n}, not {RANDOM_SEED_LEN}"),
            SeedError::Kind(k) => write!(f, "RNG kind {k} is not {RANDOM_SEED_KIND} (Mersenne-Twister, Inversion, Rejection)"),
            SeedError::AllZero => write!(f, "the Mersenne-Twister state is all zeros"),
            SeedError::Position(p) => write!(f, "Mersenne-Twister position {p} is out of range"),
        }
    }
}

impl std::error::Error for SeedError {}

impl fmt::Debug for RRng {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RRng {{ mti: {}, mt[0..3]: {:?} }}", self.mti, &self.mt[..3])
    }
}

/// `fixup()`: keep uniforms strictly inside (0, 1).
fn fixup(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.5 * I2_32M1;
    }
    if 1.0 - x <= 0.0 {
        return 1.0 - 0.5 * I2_32M1;
    }
    x
}

impl RRng {
    /// A generator after `set.seed(seed)`.
    pub fn new(seed: i32) -> RRng {
        let mut rng = RRng { mt: [0; N], mti: N as u32 + 1 };
        rng.set_seed(seed);
        rng
    }

    /// `set.seed(seed)`: `RNG_Init`'s initial scrambling (50 rounds of `seed = 69069 * seed + 1`),
    /// then 625 more rounds fill `mti` and the 624 state words, and `FixupSeeds` sets `mti = 624`.
    pub fn set_seed(&mut self, seed: i32) {
        let mut s = seed as u32;
        for _ in 0..50 {
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        // The first value fills .Random.seed[2] (mti) "for historical consistency", and is then
        // overwritten by FixupSeeds.
        s = s.wrapping_mul(69069).wrapping_add(1);
        for w in self.mt.iter_mut() {
            s = s.wrapping_mul(69069).wrapping_add(1);
            *w = s;
        }
        self.mti = N as u32;
        // An all-zero state would make R reseed from the clock; a linear congruential sequence
        // with an odd increment never yields 624 zeros in a row.
        debug_assert!(self.mt.iter().any(|&w| w != 0));
    }

    /// The state as R's `.Random.seed` (626 integers).
    pub fn random_seed(&self) -> Vec<i32> {
        let mut v = Vec::with_capacity(RANDOM_SEED_LEN);
        v.push(RANDOM_SEED_KIND);
        v.push(self.mti as i32);
        v.extend(self.mt.iter().map(|&w| w as i32));
        v
    }

    /// A generator continuing from an R `.Random.seed` (as `GetRNGstate` loads it, including
    /// `FixupSeeds`' replacement of `mti = 0` by 624).
    pub fn from_random_seed(seed: &[i32]) -> Result<RRng, SeedError> {
        if seed.len() != RANDOM_SEED_LEN {
            return Err(SeedError::Length(seed.len()));
        }
        if seed[0] != RANDOM_SEED_KIND {
            return Err(SeedError::Kind(seed[0]));
        }
        let mut mti = seed[1] as u32;
        if mti == 0 {
            mti = N as u32;
        }
        if mti > N as u32 + 1 {
            return Err(SeedError::Position(mti));
        }
        let mut mt = [0u32; N];
        for (w, &s) in mt.iter_mut().zip(&seed[2..]) {
            *w = s as u32;
        }
        if mt.iter().all(|&w| w == 0) {
            return Err(SeedError::AllZero);
        }
        Ok(RRng { mt, mti })
    }

    /// `MT_sgenrand`, used only when the state was never seeded (`mti == 625`).
    fn mt_sgenrand(&mut self, mut seed: u32) {
        for w in self.mt.iter_mut() {
            *w = seed & 0xffff_0000;
            seed = seed.wrapping_mul(69069).wrapping_add(1);
            *w |= (seed & 0xffff_0000) >> 16;
            seed = seed.wrapping_mul(69069).wrapping_add(1);
        }
        self.mti = N as u32;
    }

    /// `MT_genrand`: the next tempered word, scaled to [0, 1) by 2^-32.
    fn mt_genrand(&mut self) -> f64 {
        const MAG01: [u32; 2] = [0, MATRIX_A];
        let mut mti = self.mti as usize;
        if mti >= N {
            if mti == N + 1 {
                self.mt_sgenrand(4357);
            }
            let mt = &mut self.mt;
            let mut kk = 0;
            while kk < N - M {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
                kk += 1;
            }
            while kk < N - 1 {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
                kk += 1;
            }
            let y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
            mti = 0;
        }
        let mut y = self.mt[mti];
        mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & TEMPERING_MASK_B;
        y ^= (y << 15) & TEMPERING_MASK_C;
        y ^= y >> 18;
        self.mti = mti as u32;
        y as f64 * 2.3283064365386963e-10
    }

    /// `unif_rand()`: one uniform in (0, 1).
    pub fn unif_rand(&mut self) -> f64 {
        fixup(self.mt_genrand())
    }

    /// `rbits(bits)`: a random integer below 2^bits, built from 16-bit chunks. The loop runs for
    /// n = 0, 16, 32, ... while n <= bits, so it uses one uniform for 0-15 bits, two for 16-31
    /// (including exactly 16), three for 32-47 and four for 48-52.
    fn rbits(&mut self, bits: i32) -> f64 {
        let mut v: i64 = 0;
        let mut n = 0;
        while n <= bits {
            let v1 = (self.unif_rand() * 65536.0).floor() as i64;
            v = v.wrapping_mul(65536).wrapping_add(v1);
            n += 16;
        }
        (v & ((1i64 << bits) - 1)) as f64
    }

    /// `R_unif_index(dn)` with `sample.kind = "Rejection"`: a uniform integer in [0, dn), by
    /// rejection sampling of `rbits(ceil(log2(dn)))`. Returns 0 without drawing if `dn <= 0`.
    pub fn unif_index(&mut self, dn: f64) -> f64 {
        if dn <= 0.0 {
            return 0.0;
        }
        let bits = libm::log2(dn).ceil() as i32;
        loop {
            let dv = self.rbits(bits);
            if dn > dv {
                return dv;
            }
        }
    }

    /// `runif(1, a, b)`: NaN if a or b is not finite or b < a; `a` if a == b; else
    /// `a + (b - a) * u`.
    pub fn runif(&mut self, a: f64, b: f64) -> f64 {
        if !a.is_finite() || !b.is_finite() || b < a {
            return f64::NAN;
        }
        if a == b {
            return a;
        }
        let mut u;
        loop {
            u = self.unif_rand();
            if u > 0.0 && u < 1.0 {
                break;
            }
        }
        a + (b - a) * u
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_seed_layout() {
        let r = RRng::new(1);
        let s = r.random_seed();
        assert_eq!(s.len(), RANDOM_SEED_LEN);
        assert_eq!(s[0], RANDOM_SEED_KIND);
        assert_eq!(s[1], 624);
        let back = RRng::from_random_seed(&s).unwrap();
        assert_eq!(back.random_seed(), s);
    }

    #[test]
    fn fixup_bounds() {
        assert_eq!(fixup(0.0), 0.5 * I2_32M1);
        assert_eq!(fixup(1.0), 1.0 - 0.5 * I2_32M1);
        assert_eq!(fixup(0.25), 0.25);
    }

    #[test]
    fn mti_zero_is_replaced_and_bad_seeds_are_rejected() {
        let mut s = RRng::new(7).random_seed();
        s[1] = 0;
        assert_eq!(RRng::from_random_seed(&s).unwrap().mti, 624);
        s[0] = 10402;
        assert_eq!(RRng::from_random_seed(&s).unwrap_err(), SeedError::Kind(10402));
        assert_eq!(RRng::from_random_seed(&s[..10]).unwrap_err(), SeedError::Length(10));
    }

    #[test]
    fn unif_index_draw_counts() {
        // dn = 2^k needs k bits: one uniform up to 15 bits, two from exactly 16 bits on
        for (dn, uniforms) in [(1.0, 1), (2.0, 1), (32768.0, 1), (32769.0, 2), (65536.0, 2),
                               (65537.0, 2), (2f64.powi(31), 2), (2f64.powi(32) + 1.0, 3)] {
            let mut a = RRng::new(99);
            let mut b = a.clone();
            // the first candidate is accepted when below dn; compare positions with a stream that
            // draws the uniforms directly
            let v = a.rbits(libm::log2(dn).ceil() as i32);
            for _ in 0..uniforms {
                b.unif_rand();
            }
            assert_eq!(a.mti, b.mti, "dn = {dn}");
            assert!(v < 2f64.powi(libm::log2(dn).ceil() as i32));
        }
        assert_eq!(RRng::new(3).unif_index(1.0), 0.0);
    }
}
