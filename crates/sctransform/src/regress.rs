//! SCTransform's second regression, `vars.to.regress` (milestone M15b; docs/PLAN.md C3).
//!
//! Seurat 5.5.1's `SCTransform.default` clips the variable features' residuals to `clip.range`,
//! then calls `ScaleData(scale.data, vars.to.regress, latent.data, model.use = "linear",
//! use.umi = FALSE, do.scale, do.center, scale.max = Inf, block.size = 750)`. With
//! `vars.to.regress` set, `ScaleData.default` takes the covariates for the cells
//! (`latent.data[colnames(object), , drop = FALSE]`, the meta.data columns `SCTransform.Seurat`
//! passed, in `vars.to.regress`' order) and, with the sequential `future` plan, calls
//! `RegressOutMatrix()` once on every row:
//!
//! ```text
//! fmla <- as.formula(paste("GENE ~", paste(vars.to.regress, collapse = "+")))
//! qr <- lm(fmla, data = cbind(latent.data, GENE = data.expr[1, ]), qr = TRUE)$qr
//! for each row x: data.resid[x, ] <- qr.resid(qr = qr, y = data.expr[x, ])
//! ```
//!
//! then centres each row (`FastRowScale`: `mat - rowMeans2(mat, na.rm = TRUE)`, in blocks of 750
//! rows, which changes nothing since rows are independent) and sets `NA` (and NaN) values to 0.
//! The variable features, the residual statistics and `correct()`'s corrected counts are vst's
//! and do not change; only `scale.data` does. In the v5 path (`SCTransform.StdAssay`),
//! `FetchResiduals()` afterwards returns the existing `scale.data` rows for the variable features
//! (none needs computing), so the regressed values stay.
//!
//! [`Regression`] is the QR, made once ([`rcompat::linpack::lm_qr`] of the intercept and the
//! covariates); [`Regression::regress_rows`] applies `qr.resid` to every row, the threads sharing
//! the rows (each row's arithmetic is the same whatever the thread count).

use rcompat::linpack::{lm_qr, model_matrix, Qr, QrResid};

/// The regression of `scale.data`'s rows on per-cell covariates.
#[derive(Debug, Clone, PartialEq)]
pub struct Regression {
    /// `vars.to.regress`, in the formula's order.
    pub names: Vec<String>,
    /// `lm(GENE ~ <names>, qr = TRUE)$qr`: the intercept, then the covariates.
    pub qr: Qr,
}

impl Regression {
    /// The QR of the model matrix for `covariates` (one value per cell each, in the cells' order,
    /// named `names`). Stops as R does on a non-finite value ("NA/NaN/Inf in 'x'").
    pub fn new(names: Vec<String>, covariates: &[Vec<f64>], ncells: usize) -> Result<Regression, String> {
        if names.len() != covariates.len() || names.is_empty() {
            return Err("vars.to.regress: one covariate per name, at least one".into());
        }
        let cols: Vec<&[f64]> = covariates.iter().map(Vec::as_slice).collect();
        if cols.iter().any(|c| c.len() != ncells) {
            return Err("vars.to.regress: a covariate does not have one value per cell".into());
        }
        let x = model_matrix(&cols, ncells);
        let qr = lm_qr(x, ncells, 1 + names.len()).map_err(|e| format!("vars.to.regress ({}): {e}", names.join(", ")))?;
        Ok(Regression { names, qr })
    }

    /// The columns of the model matrix `dqrdc2` judged linearly independent (`qr$rank`) and the
    /// number it has (`p`: the intercept and the covariates).
    pub fn rank(&self) -> (usize, usize) {
        (self.qr.rank, self.qr.p)
    }

    /// `RegressOutMatrix()` on `rows` (row-major, `ncells` per row), in place: each row becomes
    /// `qr.resid(qr, y = row)`.
    pub fn regress_rows(&self, rows: &mut [f64], ncells: usize, threads: usize) {
        if ncells == 0 || rows.is_empty() {
            return;
        }
        assert_eq!(ncells, self.qr.n, "regress_rows: the rows' cells are not the covariates'");
        assert_eq!(rows.len() % ncells, 0, "regress_rows: not whole rows");
        let nrow = rows.len() / ncells;
        let r = QrResid::new(&self.qr);
        let threads = threads.max(1).min(nrow);
        let per = nrow.div_ceil(threads);
        std::thread::scope(|scope| {
            for chunk in rows.chunks_mut(per * ncells) {
                let r = &r;
                scope.spawn(move || {
                    let (mut qty, mut rsd) = (vec![0.0f64; ncells], vec![0.0f64; ncells]);
                    for row in chunk.chunks_mut(ncells) {
                        r.resid(row, &mut qty, &mut rsd);
                        row.copy_from_slice(&rsd);
                    }
                });
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_regressed_alike_on_any_number_of_threads() {
        let n = 37;
        let v: Vec<f64> = (0..n).map(|i| ((i * 13) % 17) as f64 / 3.0).collect();
        let w: Vec<f64> = (0..n).map(|i| (i as f64).sqrt()).collect();
        let reg = Regression::new(vec!["v".into(), "w".into()], &[v, w], n).unwrap();
        assert_eq!(reg.rank(), (3, 3));
        let rows: Vec<f64> = (0..11 * n).map(|k| ((k * 7919) % 1009) as f64 / 100.0 - 5.0).collect();
        let mut one = rows.clone();
        reg.regress_rows(&mut one, n, 1);
        for t in [2, 3, 4, 11, 64] {
            let mut many = rows.clone();
            reg.regress_rows(&mut many, n, t);
            assert!(one.iter().zip(&many).all(|(a, b)| a.to_bits() == b.to_bits()), "{t} threads");
        }
        assert!(Regression::new(vec!["v".into()], &[vec![1.0, f64::NAN]], 2).unwrap_err().contains("NA/NaN/Inf"));
    }
}
