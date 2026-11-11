//! The C library's math functions, called directly: R calls these same functions, and glibc's
//! versions are not all correctly rounded, so only they give R's exact results. (Rust's `f64`
//! methods may lower to other implementations or be folded at compile time.) The binary links
//! the system glibc at run time, as R does (plan A5 and A5b).

mod sys {
    extern "C" {
        pub fn pow(x: f64, y: f64) -> f64;
        pub fn log(x: f64) -> f64;
        pub fn log1p(x: f64) -> f64;
        pub fn log2(x: f64) -> f64;
        pub fn log10(x: f64) -> f64;
        pub fn logb(x: f64) -> f64;
        pub fn exp(x: f64) -> f64;
        pub fn expm1(x: f64) -> f64;
    }
}

macro_rules! wrap {
    ($($name:ident($($arg:ident),*)),* $(,)?) => {$(
        #[inline]
        pub fn $name($($arg: f64),*) -> f64 {
            // SAFETY: pure C math functions of plain doubles.
            unsafe { sys::$name($($arg),*) }
        }
    )*};
}

wrap!(pow(x, y), log(x), log1p(x), log2(x), log10(x), logb(x), exp(x), expm1(x));
