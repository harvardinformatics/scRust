//! DoubletFinder's parameter grids and the values that pass through strings (docs/PLAN.md B1).
//!
//! `paramSweep` builds its pK and pN grids with `seq()`, names its results with `paste()`, and
//! `summarizeSweep` parses pN and pK back with `as.numeric()`; `find.pK` labels its rows with
//! the factor levels (`as.character()` of the parsed values) and the workflow takes the chosen
//! pK as `as.numeric(as.character(...))`. The parsed doubles can differ in the last bit from the
//! grid's: the sweep itself uses the grid values, the final run the parsed one.

use rcompat::arith::{round, seq_by};
use rcompat::{as_character, as_numeric};

/// paramSweep's pK grid before filtering: `c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))`.
pub fn pk_grid_full() -> Vec<f64> {
    let mut g = vec![0.0005, 0.001, 0.005];
    g.extend(seq_by(0.01, 0.3, 0.01).expect("valid seq"));
    g
}

/// paramSweep's pN grid: `seq(0.05, 0.3, by = 0.05)`.
pub fn pn_grid() -> Vec<f64> {
    seq_by(0.05, 0.3, 0.05).expect("valid seq")
}

/// The pK values paramSweep keeps for `n` input cells: those with `round(pK * min.cells) >= 1`,
/// where `min.cells = round(n/(1 - 0.05) - n)`.
pub fn pk_grid(n: usize) -> Vec<f64> {
    let nf = n as f64;
    let min_cells = round(nf / (1.0 - 0.05) - nf);
    pk_grid_full().into_iter().filter(|&pk| round(pk * min_cells) >= 1.0).collect()
}

/// `round(n_real/(1 - pN) - n_real)`: the number of artificial doublets for pN.
pub fn n_doublets(n_real: usize, pn: f64) -> usize {
    let nf = n_real as f64;
    round(nf / (1.0 - pn) - nf) as usize
}

/// A value as `paste()` writes it, and as `as.numeric()` reads that string back.
#[derive(Debug, Clone, PartialEq)]
pub struct Labelled {
    /// The value used in the computation.
    pub value: f64,
    /// `as.character(value)`, as `paste()` writes it.
    pub label: String,
    /// `as.numeric(label)`, R's parser.
    pub parsed: f64,
}

impl Labelled {
    pub fn new(value: f64) -> Labelled {
        let label = as_character(value);
        let parsed = as_numeric(&label);
        Labelled { value, label, parsed }
    }

    /// The factor level of the parsed value: `as.character(as.numeric(label))`.
    pub fn level(&self) -> String {
        as_character(self.parsed)
    }
}

/// `paste("pN", pN, "pK", pK, sep = "_")`: the name of one sweep result.
pub fn sweep_name(pn: &Labelled, pk: &Labelled) -> String {
    format!("pN_{}_pK_{}", pn.label, pk.label)
}

/// `paste(prefix, pN, pK, nExp, sep = "_")`: doubletFinder's column names, `pANN_...` and
/// `DF.classifications_...`.
pub fn column_name(prefix: &str, pn: f64, pk: f64, nexp: f64) -> String {
    format!("{prefix}_{}_{}_{}", as_character(pn), as_character(pk), as_character(nexp))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grids_and_labels() {
        let g = pk_grid_full();
        assert_eq!(g.len(), 33);
        assert_eq!(Labelled::new(g[0]).label, "5e-04");
        assert_eq!(Labelled::new(g[31]).label, "0.29");
        // the grid's 0.29 is 0.01 + 28 * 0.01; R reads "0.29" back as the closest double
        assert_eq!(g[31], 0.29000000000000004);
        assert_eq!(Labelled::new(g[31]).parsed, 0.29);
        // tier T (300 cells): min.cells = 16, so pK from 0.04 on
        assert_eq!(pk_grid(300).first(), Some(&0.04));
        assert_eq!(pk_grid(300).len(), 27);
        // the P objects: min.cells 164, 623 and 4,179 for 3,108, 11,831 and 79,409 cells
        assert_eq!(pk_grid(3108).len(), 31);
        assert_eq!(pk_grid(11831).len(), 32);
        assert_eq!(pk_grid(79409).len(), 33);
        assert_eq!(n_doublets(300, 0.25), 100);
        assert_eq!(n_doublets(1122, 0.2), 280);
        assert_eq!(column_name("pANN", 0.25, 0.29, 24.0), "pANN_0.25_0.29_24");
        assert_eq!(column_name("pANN", 0.25, 0.0005, 1e5), "pANN_0.25_5e-04_1e+05");
        assert_eq!(pn_grid().iter().map(|&p| Labelled::new(p).label).collect::<Vec<_>>(),
                   ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3"]);
    }
}
