//! `JavaRandom`, the generator of Seurat's modularity optimizer (docs/PLAN.md C2 and C3c):
//! `java.util.Random`'s 48-bit linear congruential generator, as Seurat 5.5.1's
//! `ModularityOptimizer.cpp` implements it in C++, and `Arrays2::generateRandomPermutation`,
//! which every local-moving pass draws.
//!
//! *Ported from Seurat 5.5.1's `src/ModularityOptimizer.cpp` (Copyright 2021 Seurat authors, MIT),
//! itself a translation of Ludo Waltman and Nees Jan van Eck's ModularityOptimizer 1.3.0 (Java);
//! translated into Rust statement by statement. NOTICE lists the changes.*
//!
//! - `RunModularityClusteringCpp` takes the seed as an `int` and passes it to the constructor's
//!   `uint64_t`, so a negative seed sign-extends ([`JavaRandom::from_r_seed`]).
//! - `nextInt(n)` takes Java's multiply-and-shift branch when n is a power of two. Otherwise it
//!   rejects as Java does: `while (bits - val + (n - 1) < 0)` relies on 32-bit signed
//!   wrap-around, which is undefined behaviour in C++, but the compiled `Seurat.so` keeps the
//!   check (an `add` and a `js` on 32-bit registers), so the port wraps.

/// Java's multiplier.
const MULTIPLIER: u64 = 0x5DEECE66D;
/// Java's increment.
const ADDEND: u64 = 0xB;
/// The 48-bit state mask.
const MASK: u64 = (1u64 << 48) - 1;

/// `ModularityOptimizer::JavaRandom`: a `java.util.Random` stream. It also counts its `nextInt`
/// calls, and those that took the power-of-two branch (not part of the C++; they change nothing).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JavaRandom {
    seed: u64,
    draws: u64,
    pow2_draws: u64,
}

impl JavaRandom {
    /// `JavaRandom(uint64_t seed)`.
    pub fn new(seed: u64) -> JavaRandom {
        let mut r = JavaRandom { seed: 0, draws: 0, pow2_draws: 0 };
        r.set_seed(seed);
        r
    }

    /// The generator `RunModularityClusteringCpp` makes from FindClusters' `random.seed`, an R
    /// integer passed as C++ `int`: converted to `uint64_t`, so negative seeds sign-extend.
    pub fn from_r_seed(seed: i32) -> JavaRandom {
        JavaRandom::new(seed as i64 as u64)
    }

    /// `setSeed`: `(seed ^ 0x5DEECE66D) & (2^48 - 1)`.
    pub fn set_seed(&mut self, seed: u64) {
        self.seed = (seed ^ MULTIPLIER) & MASK;
    }

    /// `next(bits)`: advances the state and returns its top `bits` bits (as `int`).
    #[inline]
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        // static_cast<int>(seed >> (48 - bits)): below 2^31 for bits = 31
        (self.seed >> (48 - bits)) as i32
    }

    /// `nextInt(n)`: a uniform integer in `[0, n)`.
    ///
    /// # Panics
    /// If `n <= 0` (the C++ throws `std::out_of_range`; the optimizer never asks for it).
    #[inline]
    pub fn next_int(&mut self, n: i32) -> i32 {
        assert!(n > 0, "n must be positive");
        self.draws += 1;
        if (n & n.wrapping_neg()) == n {
            // a power of two
            self.pow2_draws += 1;
            return ((n as u64).wrapping_mul(self.next(31) as u64) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % n;
            // while (bits - val + (n - 1) < 0), in 32-bit two's complement
            if bits.wrapping_sub(val).wrapping_add(n - 1) >= 0 {
                return val;
            }
        }
    }

    /// `nextInt` calls so far, and those whose bound was a power of two.
    pub fn draws(&self) -> (u64, u64) {
        (self.draws, self.pow2_draws)
    }

    /// `Arrays2::generateRandomPermutation(nElements, random)`: the identity, then for each
    /// position i in order a swap with position `nextInt(nElements)`. Draws n numbers.
    pub fn permutation(&mut self, n: usize) -> Vec<u32> {
        let mut permutation: Vec<u32> = (0..n as u32).collect();
        let ni = i32::try_from(n).expect("more nodes than an int holds");
        for i in 0..n {
            let j = self.next_int(ni) as usize;
            permutation.swap(i, j);
        }
        permutation
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_values() {
        // new java.util.Random(0).nextInt(10), ten times, as Java gives them (M8)
        let mut r = JavaRandom::new(0);
        let v: Vec<i32> = (0..10).map(|_| r.next_int(10)).collect();
        assert_eq!(v, vec![0, 8, 9, 7, 5, 3, 1, 1, 9, 4]);
    }

    #[test]
    fn negative_seeds_sign_extend() {
        assert_eq!(JavaRandom::from_r_seed(-1), JavaRandom::new(u64::MAX));
        assert_eq!(JavaRandom::from_r_seed(12345), JavaRandom::new(12345));
    }

    #[test]
    fn counts_draws() {
        let mut r = JavaRandom::new(7);
        r.permutation(8);
        r.permutation(5);
        assert_eq!(r.draws(), (13, 8));
    }

    #[test]
    fn permutation_is_a_permutation() {
        let mut r = JavaRandom::new(12345);
        for n in [1usize, 2, 3, 16, 300] {
            let mut p = r.permutation(n);
            p.sort_unstable();
            assert_eq!(p, (0..n as u32).collect::<Vec<_>>());
        }
    }

    #[test]
    #[should_panic]
    fn zero_bound_panics() {
        JavaRandom::new(1).next_int(0);
    }
}
