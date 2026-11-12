//! The counts of the genes `vst` keeps (`umi <- umi[genes, ]`), row-major: the fits read dense
//! rows over the fitting cells (`as.matrix(umi[genes_bin, cells_step1])`) and the residuals
//! dense rows over all cells (`as.matrix(umi[genes_bin, ])`).

use scio::CsrMatrix;

/// The kept genes' rows, in the order of `genes` (row `k` is input row `genes[k]`), each with its
/// stored entries in column order, explicit zeros included.
#[derive(Debug, Clone)]
pub struct GeneRows {
    pub ncol: usize,
    pub p: Vec<usize>,
    pub j: Vec<u32>,
    pub x: Vec<f64>,
}

impl GeneRows {
    /// Rows `genes` of a row-major matrix.
    pub fn from_csr(csr: &CsrMatrix, genes: &[usize]) -> GeneRows {
        let nnz: usize = genes.iter().map(|&g| csr.p[g + 1] - csr.p[g]).sum();
        let mut out = GeneRows { ncol: csr.ncol, p: Vec::with_capacity(genes.len() + 1), j: Vec::with_capacity(nnz), x: Vec::with_capacity(nnz) };
        out.p.push(0);
        for &g in genes {
            let (cols, vals) = csr.row(g);
            out.j.extend_from_slice(cols);
            out.x.extend_from_slice(vals);
            out.p.push(out.j.len());
        }
        out
    }

    /// Rows `genes` of a column-major matrix.
    pub fn from_csc(csc: &scio::CscMatrix, genes: &[usize]) -> GeneRows {
        GeneRows::from_csr(&csc.to_csr(), genes)
    }

    pub fn nrow(&self) -> usize {
        self.p.len() - 1
    }

    pub fn row(&self, k: usize) -> (&[u32], &[f64]) {
        let (a, b) = (self.p[k], self.p[k + 1]);
        (&self.j[a..b], &self.x[a..b])
    }

    /// Row `k` over all cells, densely, into `out` (length `ncol`).
    pub fn dense_row(&self, k: usize, out: &mut [f64]) {
        out.fill(0.0);
        let (cols, vals) = self.row(k);
        for (&j, &v) in cols.iter().zip(vals) {
            out[j as usize] = v;
        }
    }

    /// Row `k` over the cells `cells` in that order. `pos[c]` is the position of cell `c` in
    /// `cells` (`usize::MAX` if absent), so the row's stored entries are scattered once.
    pub fn dense_row_cells(&self, k: usize, pos: &[usize], out: &mut [f64]) {
        out.fill(0.0);
        let (cols, vals) = self.row(k);
        for (&j, &v) in cols.iter().zip(vals) {
            let p = pos[j as usize];
            if p != usize::MAX {
                out[p] = v;
            }
        }
    }

    /// `median(umi@x)`: the median of every stored value of the kept genes, explicit zeros
    /// included (`get_nz_median2(umi)`), as R's `median` computes it.
    pub fn stored_median(&self) -> f64 {
        rcompat::quantile::median(&self.x)
    }
}

/// `pos[c]` = position of cell `c` in `cells`, `usize::MAX` for cells not in it.
pub fn positions(cells: &[usize], ncol: usize) -> Vec<usize> {
    let mut pos = vec![usize::MAX; ncol];
    for (k, &c) in cells.iter().enumerate() {
        pos[c] = k;
    }
    pos
}
