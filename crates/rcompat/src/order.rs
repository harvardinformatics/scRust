//! `order()` for one numeric key, as R 4.5.3 computes it with its default method: the radix sort
//! of `src/main/radixsort.c` (`order(x)` picks method "radix" for numeric, factor and logical
//! keys), which is stable.
//!
//! The sort key is `dtwiddle(order * x)` with no rounding (`order()` passes `retGrp = FALSE`):
//! `+0` and `-0` are equal, and with `na.last = TRUE` every NA and NaN sorts last, tied with each
//! other. `decreasing = TRUE` negates the key rather than reversing the result, so ties keep their
//! original order (lower index first) in both directions.

/// The permutation `order(x, decreasing = decreasing)` (with `na.last = TRUE`), 0-based.
pub fn order(x: &[f64], decreasing: bool) -> Vec<usize> {
    let sign = if decreasing { -1.0 } else { 1.0 };
    let keys: Vec<u64> = x.iter().map(|&v| twiddle(sign * v)).collect();
    let mut idx: Vec<usize> = (0..x.len()).collect();
    idx.sort_by_key(|&i| keys[i]); // stable
    idx
}

/// `dtwiddle` with `dround = 0` and `nalast = 1`: an unsigned key whose order is the numeric
/// order, with NA and NaN largest.
fn twiddle(v: f64) -> u64 {
    if v.is_nan() {
        return u64::MAX;
    }
    let bits = if v == 0.0 { 0 } else { v.to_bits() };
    if bits & (1 << 63) != 0 {
        !bits
    } else {
        bits ^ (1 << 63)
    }
}

/// R's `NA_integer_`.
pub const NA_INTEGER: i32 = i32::MIN;

/// The permutation `order(x, decreasing = decreasing)` for integers (with `na.last = TRUE`),
/// 0-based: stable, `NA_integer_` last.
pub fn order_int(x: &[i32], decreasing: bool) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..x.len()).collect();
    // key: NA last; otherwise the value, negated when decreasing (as i64, so -i32::MAX is fine)
    let key = |v: i32| -> (bool, i64) {
        if v == NA_INTEGER {
            (true, 0)
        } else if decreasing {
            (false, -(v as i64))
        } else {
            (false, v as i64)
        }
    };
    idx.sort_by_key(|&i| key(x[i]));
    idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_in_both_directions() {
        let x = [2.0, 1.0, 2.0, f64::NAN, -0.0, 0.0, 1.0];
        assert_eq!(order(&x, false), vec![4, 5, 1, 6, 0, 2, 3]);
        assert_eq!(order(&x, true), vec![0, 2, 1, 6, 4, 5, 3]);
        assert_eq!(order_int(&[3, NA_INTEGER, 1, 3], true), vec![0, 3, 2, 1]);
    }
}
