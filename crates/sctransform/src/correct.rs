//! The corrected counts (docs/PLAN.md C3f, milestone M13): sctransform 0.4.3's `correct()`, as
//! `vst()` calls it when Seurat 5.5.1's `SCTransform()` passes `return_corrected_umi = TRUE`
//! (`do.correct.umi`, the default), with no `scale_factor`:
//!
//! ```r
//! cell_attr[, latent_var] <- apply(cell_attr[, latent_var, drop=FALSE], 2,
//!                                  function(x) rep(median(x), length(x)))
//! regressor_data <- model.matrix(get_model_formula(x$model_str), cell_attr)   # (Intercept), log_umi
//! for (i in 1:max_bin) {                                       # bins of 500 genes, in vst's order
//!   pearson_residual <- data[genes_bin, ]                      # vst's residuals, not yet clipped
//!   coefs <- x$model_pars_fit[genes_bin, -1, drop=FALSE]
//!   theta <- x$model_pars_fit[genes_bin, 1]
//!   mu <- exp(tcrossprod(coefs, regressor_data))
//!   variance <- mu + mu^2 / theta
//!   corrected_data[genes_bin, ] <- mu + pearson_residual * sqrt(variance)
//! }
//! corrected_data <- round(corrected_data, 0)                   # nearbyint(): half to even
//! corrected_data[corrected_data < 0] <- 0
//! rv$umi_corrected <- make.sparse(mat = rv$umi_corrected)      # zeros dropped
//! ```
//!
//! `vst()` calls it on its own residual matrix `rv$y`: the Pearson residuals of every kept gene
//! and cell, `(y - mu) / sqrt(max(mu + mu^2 / theta, min_var))`, before it clips them at
//! `+-sqrt(N)`. The port makes them in the same pass as S2d ([`crate::residuals_corrected`]), one
//! bin's products at a time, and keeps only the corrected counts, sparse; it never holds a dense
//! genes x cells matrix. `mu`'s product is the OpenBLAS `dgemm` R calls, with R's shapes, so the
//! counts are R's bit for bit under R's kernel (a different kernel can move a value across .5).
//!
//! Seurat then makes the SCT assay's `counts` of them, its `data` of their `log1p` (glibc's), and
//! adding the assay computes `nCount_SCT` (`colSums`) and `nFeature_SCT` (the entries > 0).

use std::io;
use std::path::Path;

/// The corrected counts of one `SCTransform()` call: genes (the kept genes, in their order) x
/// cells, compressed by column (cell) as Seurat's `dgCMatrix` stores them, rows increasing
/// within each column; every stored value is a positive integer.
#[derive(Debug, Clone, PartialEq)]
pub struct Corrected {
    pub ngenes: usize,
    pub ncells: usize,
    /// R's `median()` of the cells' `log_umi`: the depth every cell is corrected to.
    pub median_log_umi: f64,
    /// vst's clip range of the residuals, and how many of the residuals `correct()` read lie
    /// beyond it (it reads them before vst clips them).
    pub res_clip: (f64, f64),
    pub beyond_clip: u64,
    pub p: Vec<usize>,
    pub i: Vec<u32>,
    pub x: Vec<u32>,
}

/// `log1p` of a corrected count, as R computes Seurat's `data` (glibc's `log1p`): small counts
/// from a table, the same values.
pub struct Log1p {
    table: Vec<f64>,
}

impl Log1p {
    pub fn new() -> Log1p {
        Log1p { table: (0..1u32 << 16).map(|v| rcompat::libm::log1p(v as f64)).collect() }
    }

    #[inline]
    pub fn of(&self, v: u32) -> f64 {
        match self.table.get(v as usize) {
            Some(&y) => y,
            None => rcompat::libm::log1p(v as f64),
        }
    }
}

impl Default for Log1p {
    fn default() -> Self {
        Log1p::new()
    }
}

impl Corrected {
    /// From each gene's corrected counts (cell, count), in the kept genes' order: the
    /// transposition into compressed columns. Each gene's list is freed once placed.
    pub(crate) fn from_genes(genes: Vec<Vec<(u32, u32)>>, ncells: usize, median_log_umi: f64, res_clip: (f64, f64), beyond_clip: u64) -> Corrected {
        let ngenes = genes.len();
        let mut p = vec![0usize; ncells + 1];
        for g in &genes {
            for &(c, _) in g {
                p[c as usize + 1] += 1;
            }
        }
        for c in 0..ncells {
            p[c + 1] += p[c];
        }
        let nnz = p[ncells];
        let mut i = vec![0u32; nnz];
        let mut x = vec![0u32; nnz];
        let mut next = p.clone();
        for (gi, g) in genes.into_iter().enumerate() {
            for (c, v) in g {
                let k = next[c as usize];
                i[k] = gi as u32;
                x[k] = v;
                next[c as usize] += 1;
            }
        }
        Corrected { ngenes, ncells, median_log_umi, res_clip, beyond_clip, p, i, x }
    }

    pub fn nnz(&self) -> usize {
        self.x.len()
    }

    /// `nCount_SCT`: each cell's column sum (`Matrix::colSums`; integers, so exact in any order).
    pub fn n_count(&self) -> Vec<f64> {
        (0..self.ncells).map(|c| self.x[self.p[c]..self.p[c + 1]].iter().map(|&v| v as u64).sum::<u64>() as f64).collect()
    }

    /// `nFeature_SCT`: each cell's entries above 0 (every stored entry is).
    pub fn n_feature(&self) -> Vec<u32> {
        (0..self.ncells).map(|c| (self.p[c + 1] - self.p[c]) as u32).collect()
    }

    /// The smallest stored count (R's reference runs record it; 1 when any entry is stored).
    pub fn min(&self) -> Option<u32> {
        self.x.iter().copied().min()
    }

    /// The SHA-256 of the stored values as doubles, little-endian, in stored order: R's
    /// `raw_sha256(counts@x)` (`validation/R/lib/clustering.R`), with `log1p` applied for the
    /// `data` layer's digest (`raw_sha256(data@x)`).
    pub fn x_sha256(&self, log1p: bool) -> String {
        let l = Log1p::new();
        let mut h = scio::sha256::Sha256::new();
        let mut buf: Vec<u8> = Vec::with_capacity(8 << 16);
        for chunk in self.x.chunks(1 << 16) {
            buf.clear();
            for &v in chunk {
                let y = if log1p { l.of(v) } else { v as f64 };
                buf.extend_from_slice(&y.to_le_bytes());
            }
            h.update(&buf);
        }
        h.hex()
    }

    /// Writes the SCT assay's `counts` and `data` layers as Matrix Market files in `dir`:
    /// `corrected_counts.mtx.gz` (integer) and `log1p_corrected_counts.mtx.gz` (real), entries in
    /// stored (column-major) order, values as R's reference runs write them (`%.17g`), each with its
    /// comment lines after the banner; and their rows and columns, `sct_features.tsv.gz` (`genes`,
    /// the kept genes' names, one per line) and `sct_cells.tsv.gz` (`cells`). Returns the file
    /// names written. The bytes do not depend on `threads`.
    pub fn write(&self, dir: &Path, genes: &[String], cells: &[String], counts_comments: &[String], data_comments: &[String], threads: usize) -> io::Result<Vec<String>> {
        if genes.len() != self.ngenes || cells.len() != self.ncells {
            return Err(io::Error::other(format!("corrected counts of {} genes x {} cells, given {} gene and {} cell names", self.ngenes, self.ncells, genes.len(), cells.len())));
        }
        std::fs::create_dir_all(dir)?;
        let x = &self.x;
        scio::mtx::write_gz_with(&dir.join("corrected_counts.mtx.gz"), self.ngenes, self.ncells, &self.p, &self.i, &|k| x[k] as f64, counts_comments, threads)?;
        // log1p of the counts: each distinct small count's text made once (%.17g, as R's fmt17)
        let l = Log1p::new();
        let top = self.x.iter().copied().max().unwrap_or(0).min(1 << 16) as usize;
        let text: Vec<String> = (0..=top).map(|v| scio::fmt17(l.of(v as u32))).collect();
        // R's mtx_write_plain: "integer" when every value is integer-valued (only with no entry)
        let field = if self.x.is_empty() { "integer" } else { "real" };
        scio::mtx::write_gz_formatted(
            &dir.join("log1p_corrected_counts.mtx.gz"),
            self.ngenes,
            self.ncells,
            &self.p,
            &self.i,
            field,
            &|k, out| match text.get(x[k] as usize) {
                Some(t) => out.extend_from_slice(t.as_bytes()),
                None => out.extend_from_slice(scio::fmt17(l.of(x[k])).as_bytes()),
            },
            data_comments,
            threads,
        )?;
        let lines = |v: &[String]| {
            let mut t = String::with_capacity(v.iter().map(|s| s.len() + 1).sum());
            for s in v {
                t.push_str(s);
                t.push('\n');
            }
            scio::deflate::gzip(t.as_bytes())
        };
        std::fs::write(dir.join("sct_features.tsv.gz"), lines(genes))?;
        std::fs::write(dir.join("sct_cells.tsv.gz"), lines(cells))?;
        Ok(["corrected_counts.mtx.gz", "log1p_corrected_counts.mtx.gz", "sct_features.tsv.gz", "sct_cells.tsv.gz"].iter().map(|s| s.to_string()).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transposition_and_column_statistics() {
        // 3 genes x 4 cells: gene 0 in cells 1 and 3, gene 1 in cell 1, gene 2 in cells 0 and 1
        let genes = vec![vec![(1, 2), (3, 7)], vec![(1, 1)], vec![(0, 4), (1, 9)]];
        let c = Corrected::from_genes(genes, 4, 2.5, (-2.0, 2.0), 3);
        assert_eq!(c.p, vec![0, 1, 4, 4, 5]);
        assert_eq!(c.i, vec![2, 0, 1, 2, 0]);
        assert_eq!(c.x, vec![4, 2, 1, 9, 7]);
        assert_eq!(c.n_count(), vec![4.0, 12.0, 0.0, 7.0]);
        assert_eq!(c.n_feature(), vec![1, 3, 0, 1]);
        assert_eq!((c.nnz(), c.min()), (5, Some(1)));
        // R's raw_sha256 of the doubles
        let bytes: Vec<u8> = [4.0f64, 2.0, 1.0, 9.0, 7.0].iter().flat_map(|v| v.to_le_bytes()).collect();
        assert_eq!(c.x_sha256(false), scio::sha256::hex(&bytes));
        let lb: Vec<u8> = [4.0f64, 2.0, 1.0, 9.0, 7.0].iter().flat_map(|&v| rcompat::libm::log1p(v).to_le_bytes()).collect();
        assert_eq!(c.x_sha256(true), scio::sha256::hex(&lb));
        let l = Log1p::new();
        assert_eq!(l.of(1).to_bits(), std::f64::consts::LN_2.to_bits());
        assert_eq!(l.of(70_000).to_bits(), rcompat::libm::log1p(70_000.0).to_bits());
    }
}
