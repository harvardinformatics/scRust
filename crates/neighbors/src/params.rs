//! FindNeighbors' arguments (docs/PLAN.md A7b and C3, "Parameters"): every one the neighbour
//! search uses, with Seurat 5.5.1's default (`FindNeighbors.Seurat`), and the parsers for their
//! text forms (`--dims 1:30`, `--prune 1/15`).

use crate::annoy::Metric;

/// `FindNeighbors(object, reduction = "pca", dims, k.param, prune.SNN, nn.method = "annoy",
/// n.trees, annoy.metric, l2.norm = FALSE)`. Other code paths (`l2.norm = TRUE`,
/// `nn.method = "rann"` or `"hnsw"`, `annoy.metric = "hamming"`, a distance matrix) are not
/// ported, and have no field here.
#[derive(Debug, Clone, PartialEq)]
pub struct NeighborParams {
    /// `dims`: the embedding columns used, 1-based and in the given order (Seurat's
    /// `data.use[, dims]`). Seurat's default is 1:10; the workflow passes 1:30.
    pub dims: Vec<usize>,
    /// `k.param` (20). Seurat lowers it to the number of cells minus 1 when there are fewer cells.
    pub k_param: usize,
    /// `n.trees` (50).
    pub n_trees: i32,
    /// `annoy.metric` (euclidean).
    pub metric: Metric,
    /// `prune.SNN` (1/15): SNN entries `<` it are removed.
    pub prune_snn: f64,
    /// AnnoySearch's `search.k` (-1: k times the number of trees). FindNeighbors does not pass
    /// it, so Seurat always uses -1; it is a parameter here for the sweeps of plan A7b.
    pub search_k: i32,
}

impl Default for NeighborParams {
    fn default() -> Self {
        NeighborParams { dims: (1..=10).collect(), k_param: 20, n_trees: 50, metric: Metric::Euclidean, prune_snn: 1.0 / 15.0, search_k: -1 }
    }
}

impl NeighborParams {
    /// The workflow's call, `FindNeighbors(obj, dims = 1:30)`, with every other default.
    pub fn workflow() -> NeighborParams {
        NeighborParams { dims: (1..=30).collect(), ..NeighborParams::default() }
    }

    /// Checks the parameters against an embedding of `ncol` columns, as Seurat would.
    pub fn check(&self, ncol: usize) -> Result<(), String> {
        if self.dims.is_empty() {
            return Err("dims: no dimensions".into());
        }
        if let Some(&d) = self.dims.iter().find(|&&d| d == 0) {
            return Err(format!("dims: {d} is not a dimension (dims are 1-based)"));
        }
        let max = *self.dims.iter().max().unwrap();
        if max > ncol {
            // Seurat's message
            return Err(format!("More dimensions specified in dims ({max}) than have been computed ({ncol})"));
        }
        if self.k_param == 0 {
            return Err("k.param must be at least 1".into());
        }
        if self.n_trees < 1 {
            return Err(format!("n.trees = {}: at least one tree is needed", self.n_trees));
        }
        if self.prune_snn.is_nan() {
            return Err("prune.SNN is NaN".into());
        }
        if self.search_k < -1 || self.search_k == 0 {
            return Err(format!("search.k = {}: expected -1 or a positive number", self.search_k));
        }
        Ok(())
    }
}

/// `dims` as R would write them: `a:b` (decreasing when a > b, as R's `:`), single numbers and
/// comma-separated lists of both, e.g. `1:30` or `1:10,15`. Dimensions are 1-based.
pub fn parse_dims(s: &str) -> Result<Vec<usize>, String> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let part = part.trim();
        let num = |t: &str| t.trim().parse::<usize>().map_err(|_| format!("dims: '{t}' is not a positive integer"));
        if let Some((a, b)) = part.split_once(':') {
            let (a, b) = (num(a)?, num(b)?);
            if a <= b {
                out.extend(a..=b);
            } else {
                out.extend((b..=a).rev());
            }
        } else {
            out.push(num(part)?);
        }
    }
    if out.contains(&0) {
        return Err("dims: 0 is not a dimension (dims are 1-based)".into());
    }
    Ok(out)
}

/// A number as R reads it (`as.numeric`, R's own decimal parser), or a fraction `a/b`, divided
/// in double as R evaluates `1/15`.
pub fn parse_number(s: &str) -> Result<f64, String> {
    let read = |t: &str| {
        let v = rcompat::as_numeric(t.trim());
        if v.is_nan() {
            Err(format!("'{t}' is not a number"))
        } else {
            Ok(v)
        }
    };
    match s.split_once('/') {
        Some((a, b)) => Ok(read(a)? / read(b)?),
        None => read(s),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dims_like_r() {
        assert_eq!(parse_dims("1:3").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_dims("3:1").unwrap(), vec![3, 2, 1]);
        assert_eq!(parse_dims("1:2,5, 7").unwrap(), vec![1, 2, 5, 7]);
        assert!(parse_dims("0:3").is_err());
        assert!(parse_dims("a").is_err());
        assert_eq!(NeighborParams::workflow().dims.len(), 30);
    }

    #[test]
    fn fractions_are_divided() {
        assert_eq!(parse_number("1/15").unwrap().to_bits(), (1.0f64 / 15.0).to_bits());
        assert_eq!(parse_number("0.2").unwrap(), 0.2);
        assert_eq!(parse_number("0").unwrap(), 0.0);
        assert!(parse_number("x").is_err());
    }

    #[test]
    fn checks_like_seurat() {
        let p = NeighborParams::workflow();
        assert!(p.check(50).is_ok());
        assert!(p.check(20).unwrap_err().starts_with("More dimensions"));
    }
}
