//! `Kiss64Random`, the generator RcppAnnoy 0.0.23 gives every Annoy index (`kissrandom.h`,
//! Annoy 1.17.3; docs/PLAN.md C2).
//!
//! *Ported from Annoy 1.17.3's `kissrandom.h` (Copyright (c) 2013 Spotify AB, Apache License
//! 2.0, `third_party/licenses/Apache-2.0.txt`): translated into Rust; NOTICE lists the changes.*
//!
//! George Marsaglia's 64-bit KISS: a linear congruential generator, a xorshift and a
//! multiply-with-carry, added. Every operation is on `uint64_t` and wraps modulo 2^64.
//! Seurat never calls `setSeed`, so every index starts from [`Kiss64Random::DEFAULT_SEED`].

/// The generator's state, as `struct Kiss64Random` holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Kiss64Random {
    x: u64,
    y: u64,
    z: u64,
    c: u64,
}

impl Kiss64Random {
    /// `Kiss64Random::default_seed`, which every Annoy index uses unless `set_seed` is called.
    pub const DEFAULT_SEED: u64 = 1234567890987654321;

    /// `Kiss64Random(seed)`: the seed goes into `x`, the other words take fixed values.
    pub fn new(seed: u64) -> Kiss64Random {
        Kiss64Random { x: seed, y: 362436362436362436, z: 1066149217761810, c: 123456123456123456 }
    }

    /// `kiss()`: the next 64-bit value.
    #[inline]
    pub fn kiss(&mut self) -> u64 {
        // linear congruence: z = 6906969069LL * z + 1234567 (the signed constant converts to
        // uint64_t, so the product wraps)
        self.z = self.z.wrapping_mul(6906969069).wrapping_add(1234567);
        // xorshift
        self.y ^= self.y << 13;
        self.y ^= self.y >> 17;
        self.y ^= self.y << 43;
        // multiply-with-carry: t = (x << 58) + c; c = x >> 6; x += t; c += (x < t)
        let t = (self.x << 58).wrapping_add(self.c);
        self.c = self.x >> 6;
        self.x = self.x.wrapping_add(t);
        self.c = self.c.wrapping_add((self.x < t) as u64);
        self.x.wrapping_add(self.y).wrapping_add(self.z)
    }

    /// `flip()`: `kiss() & 1`.
    #[inline]
    pub fn flip(&mut self) -> usize {
        (self.kiss() & 1) as usize
    }

    /// `index(n)`: `kiss() % n`, for `n >= 1`.
    #[inline]
    pub fn index(&mut self, n: usize) -> usize {
        (self.kiss() % n as u64) as usize
    }

    /// `set_seed(seed)`: replaces `x` only, as the C++ does.
    pub fn set_seed(&mut self, seed: u64) {
        self.x = seed;
    }
}

impl Default for Kiss64Random {
    fn default() -> Self {
        Kiss64Random::new(Kiss64Random::DEFAULT_SEED)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_values_of_the_default_seed() {
        // the start of tests/fixtures/clustering/k0/kiss64_kiss.tsv (RcppAnnoy's own header)
        let mut r = Kiss64Random::default();
        assert_eq!(r.kiss(), 8932985056925012148);
        assert_eq!(r.kiss(), 5710300428094272059);
        assert_eq!(r.kiss(), 18342510866933518593);
    }

    #[test]
    fn set_seed_replaces_x_only() {
        let mut a = Kiss64Random::new(7);
        a.kiss();
        let mut b = a.clone();
        b.set_seed(99);
        assert_eq!((b.y, b.z, b.c), (a.y, a.z, a.c));
        assert_eq!(b.x, 99);
    }
}
