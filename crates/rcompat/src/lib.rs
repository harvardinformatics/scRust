//! R-compatible random numbers and R arithmetic semantics (docs/PLAN.md A5; milestone M1).
//!
//! Exact agreement with R depends on copying R's behaviour in a few places. Each piece here is a
//! port of R 4.5.3's own code (the sources fetched by `third_party/fetch_sources.sh`), tested bit
//! for bit against vectors that R generates (`validation/R/rcompat_vectors.R`, fixtures in
//! `tests/fixtures/rcompat/`):
//!
//! - [`RRng`]: R's Mersenne-Twister stream, including `set.seed`'s initial scrambling and
//!   `unif_rand`'s fixup, `R_unif_index` with `sample.kind = "Rejection"`, and `.Random.seed`
//!   import and export. One generator is passed explicitly (`&mut RRng`) through every call, like
//!   R's single global stream;
//! - [`sample`]: `sample()` and `sample.int()` with and without replacement, with and without
//!   weights (`revsort`, Walker's alias method, `ProbSampleNoReplace`), and the hash-based
//!   algorithm R switches to above 1e7 elements;
//! - [`normal`]: `qnorm` (AS 241 and R's extreme-tail formulas) and `rnorm` by inversion;
//! - [`summary`]: `sum()` and `mean()` in 80-bit long double, two-pass `var()` and `sd()`, all
//!   through the exact x87 emulation in [`ldouble`];
//! - [`arith`]: `round()` half to even and `round(x, digits)`, `^` (`x^2` is `x*x`, the rest is
//!   glibc's `pow` with R's special cases), and `seq()` by step and by length;
//! - [`order`]: stable `order()`, ties to the lower index in both directions;
//! - [`format`]: numbers as `as.character()` and `paste()` write them (15 significant digits,
//!   fixed or scientific: `5e-04` but `0.001`);
//! - [`strtod`]: R's own decimal parser (`R_strtod`, not correctly rounded) and `as.numeric()`
//!   of strings;
//! - [`fft`]: `fft()`, Singleton's mixed-radix transform, as R computes it;
//! - [`normal::dnorm`] and [`approx`]: `dnorm()`, and `approx()`/`approxfun()` interpolation;
//! - [`blas`]: the conda env's OpenBLAS through FFI (`dgemm`, `dgemv`, `ddot`, `dnrm2`) and R's
//!   `tcrossprod`, so products round exactly as R's do on the same kernel (M4);
//! - [`quantile`]: `quantile(type = 7)`, `IQR`, `median`, `mad` (M4);
//! - [`density`]: `density()` with `bw = "nrd"` and its binning (M4);
//! - [`gamma`]: R's `lgammafn`, `digamma`, `trigamma`, and glibc's `lgamma` (M4);
//! - [`port`]: `nlminb()`, R's PORT library (`DRMNHB` and its call tree) (M4);
//! - [`matrixstats`]: matrixStats' `rowMeans2`, `colMeans2`, `rowVars`, `weightedMedian`, and
//!   R's `R_qsort_I` (M4).
//!
//! Indices are 0-based throughout (R's minus one). Functions that call the C math library ([`libm`])
//! use the system glibc, as R does.

pub mod approx;
pub mod arith;
pub mod blas;
pub mod density;
pub mod fft;
pub mod format;
pub mod gamma;
pub mod ldouble;
pub mod libm;
pub mod matrixstats;
pub mod nbinom;
pub mod normal;
pub mod optimize;
pub mod order;
pub mod port;
mod powl10;
pub mod quantile;
pub mod rng;
pub mod sample;
pub mod strtod;
pub mod summary;

pub use arith::{r_pow, round, round_digits, seq_by, seq_length_out};
pub use format::as_character;
pub use ldouble::F80;
pub use normal::{dnorm, qnorm, rnorm};
pub use order::order;
pub use rng::RRng;
pub use sample::{sample_int, sample_int_weighted};
pub use strtod::{as_numeric, r_strtod};
pub use summary::{is_na, mean, sd, sum, var, NA_REAL};
