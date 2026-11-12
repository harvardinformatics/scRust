//! Seurat 5.5.1's `RowVar` (`src/data_manipulation.cpp`), which decides `PrepDR`'s
//! zero-variance filter and `RunPCA`'s `total.variance`:
//!
//! ```cpp
//! Eigen::ArrayXd r = x.row(i).array();
//! double rowMean = r.mean();
//! out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
//! ```
//!
//! Seurat was built against RcppEigen 0.3.4.0.2 (Eigen 3.4.0, its conda recipe) with
//! `-march=nocona`, so Eigen vectorises with SSE2 packets of two doubles. `r` is a fresh,
//! aligned `ArrayXd`, and both sums take Eigen's `redux_impl<..., LinearVectorizedTraversal,
//! NoUnrolling>` (`Eigen/src/Core/Redux.h`): two packet accumulators over blocks of four, then a
//! last whole packet, the packet's two lanes added (`predux`), then the odd last element. The
//! squares are formed element by element, as the packets form them.

/// Eigen 3.4's vectorised sum of `size` values `x(0..size)` with SSE2 packets (alignment 0).
fn eigen_sum(size: usize, x: impl Fn(usize) -> f64) -> f64 {
    const PS: usize = 2; // doubles per SSE2 packet
    let aligned_size2 = (size / (2 * PS)) * (2 * PS);
    let aligned_size = (size / PS) * PS;
    if aligned_size == 0 {
        // too small to vectorise: one element (`size` is at least 1)
        let mut res = x(0);
        for i in 1..size {
            res += x(i);
        }
        return res;
    }
    let mut p0 = [x(0), x(1)];
    if aligned_size > PS {
        let mut p1 = [x(2), x(3)];
        let mut i = 2 * PS;
        while i < aligned_size2 {
            p0 = [p0[0] + x(i), p0[1] + x(i + 1)];
            p1 = [p1[0] + x(i + 2), p1[1] + x(i + 3)];
            i += 2 * PS;
        }
        p0 = [p0[0] + p1[0], p0[1] + p1[1]];
        if aligned_size > aligned_size2 {
            p0 = [p0[0] + x(aligned_size2), p0[1] + x(aligned_size2 + 1)];
        }
    }
    let mut res = p0[0] + p0[1];
    for i in aligned_size..size {
        res += x(i);
    }
    res
}

/// `RowVar` of one row (at least one value): its mean, the sum of squared deviations, divided by
/// `ncol - 1` (converted from Eigen's integer `Index`).
pub fn row_var(row: &[f64]) -> f64 {
    let n = row.len();
    assert!(n > 0, "row_var: empty row");
    let mean = eigen_sum(n, |i| row[i]) / n as f64;
    let ss = eigen_sum(n, |i| {
        let d = row[i] - mean;
        d * d
    });
    ss / (n as i64 - 1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigen_order() {
        // the order is visible with values whose sums round: 1e16 absorbs a 1 unless 1s pair up
        let x = [1e16, 1.0, 1.0, 1.0, -1e16, 1.0, 1.0];
        // lanes: (1e16 + 1 + -1e16 + ...) blocks of four: p0 = [1e16, 1] + ...; see eigen_sum
        let s = eigen_sum(x.len(), |i| x[i]);
        // p0 = [1e16, 1], p1 = [1, 1]; no second block (aligned_size2 = 4);
        // p0 = [1e16 + 1, 1 + 1] = [1e16, 2]; extra packet [-1e16, 1]: [0, 3]; res = 3; + 1 = 4
        assert_eq!(s, 4.0);
        let seq: f64 = x.iter().sum();
        assert_ne!(s, seq); // sequential summation gives a different value here
        assert_eq!(eigen_sum(1, |_| 5.0), 5.0);
        assert_eq!(eigen_sum(3, |i| [1.0, 2.0, 4.0][i]), 7.0);
    }

    #[test]
    fn variances() {
        assert_eq!(row_var(&[1.0, 2.0, 3.0, 4.0]), 5.0 / 3.0);
        assert_eq!(row_var(&[2.5; 9]), 0.0);
        // a constant row whose mean does not round has zero variance; one whose mean rounds
        // (0.1 x 316) does not, which is why the summation order matters for PrepDR's filter
        // (checked against Seurat's own RowVar in tests/r_vectors.rs)
        assert_eq!(row_var(&[0.0009765625; 316]), 0.0);
        assert!(row_var(&[0.1; 316]) > 0.0);
    }
}
