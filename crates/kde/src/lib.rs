//! Kernel density estimation and moments, as DoubletFinder uses them to choose pK
//! (docs/PLAN.md B3; milestone M3).
//!
//! This crate will hold:
//!
//! - `KernSmooth::bkde` (KernSmooth 2.23-27): linear binning, a Normal kernel truncated at 4
//!   bandwidths, and convolution through R's `fft`;
//! - R's `approxfun` (linear, `rule = 1`, `ties = mean`);
//! - DoubletFinder's skewness, kurtosis and bimodality coefficient, with R's `sum` and `mean`
//!   semantics.

#[cfg(test)]
mod tests {
    /// Placeholder until M3: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "kde");
    }
}
