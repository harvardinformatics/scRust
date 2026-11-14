//! FindClusters' arguments (docs/PLAN.md A7b and C3, "Parameters"): every one the modularity
//! optimizer path uses, with Seurat 5.5.1's default (`FindClusters.Seurat`), and their parsers.

use crate::optimizer::OptimizerParams;

/// `FindClusters(object, graph.name, modularity.fxn, resolution, algorithm, n.start, n.iter,
/// random.seed, group.singletons)`. Leiden (`algorithm = 4`), `initial.membership`,
/// `node.sizes` and `edge.file.name` take other code paths and are not ported.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterParams {
    /// `graph.name`: the SNN graph's name, the prefix of each column (Seurat's default is the
    /// default assay's, `SCT_snn` after SCTransform).
    pub graph_name: String,
    /// `modularity.fxn` (1).
    pub modularity_fxn: i32,
    /// `resolution` (0.8): one FindClusters call may take several.
    pub resolutions: Vec<f64>,
    /// `algorithm` (1, the original Louvain algorithm; 2 with multilevel refinement; 3 smart
    /// local moving).
    pub algorithm: i32,
    /// `n.start` (10).
    pub n_start: i32,
    /// `n.iter` (10).
    pub n_iter: i32,
    /// `random.seed` (0; the workflow passes 12345).
    pub random_seed: i32,
    /// `group.singletons` (TRUE).
    pub group_singletons: bool,
}

impl Default for ClusterParams {
    fn default() -> Self {
        ClusterParams {
            graph_name: "SCT_snn".into(),
            modularity_fxn: 1,
            resolutions: vec![0.8],
            algorithm: 1,
            n_start: 10,
            n_iter: 10,
            random_seed: 0,
            group_singletons: true,
        }
    }
}

impl ClusterParams {
    /// The workflow's call, `FindClusters(obj, random.seed = 12345)`, with every other default.
    pub fn workflow() -> ClusterParams {
        ClusterParams { random_seed: 12345, ..ClusterParams::default() }
    }

    /// The optimizer's arguments for one resolution.
    pub fn optimizer(&self, resolution: f64) -> OptimizerParams {
        OptimizerParams {
            modularity_function: self.modularity_fxn,
            resolution,
            algorithm: self.algorithm,
            n_random_starts: self.n_start,
            n_iterations: self.n_iter,
            random_seed: self.random_seed,
        }
    }

    /// `paste0("res.", r)`: the resolution as R's `as.character` writes it (15 significant
    /// digits).
    pub fn res_name(resolution: f64) -> String {
        format!("res.{}", rcompat::as_character(resolution))
    }

    /// The metadata column FindClusters.Seurat names: `paste(graph.name, "res.<r>", sep = "_")`.
    pub fn column_name(&self, resolution: f64) -> String {
        format!("{}_{}", self.graph_name, ClusterParams::res_name(resolution))
    }

    /// Checks the arguments as FindClusters and `RunModularityClusteringCpp` would, and refuses
    /// the paths that are not ported.
    pub fn check(&self) -> Result<(), String> {
        if self.algorithm == 4 {
            return Err("algorithm 4 (Leiden) runs through leidenbase, which is not ported (Leiden is not part of Step 2)".into());
        }
        if !(1..=3).contains(&self.algorithm) {
            return Err(format!("algorithm {}: expected 1 (Louvain), 2 (Louvain with multilevel refinement) or 3 (SLM)", self.algorithm));
        }
        if self.resolutions.is_empty() {
            return Err("no resolution".into());
        }
        if let Some(r) = self.resolutions.iter().find(|r| !r.is_finite()) {
            return Err(format!("resolution {r} is not a finite number"));
        }
        let names: Vec<String> = self.resolutions.iter().map(|&r| ClusterParams::res_name(r)).collect();
        for (i, a) in names.iter().enumerate() {
            if names[..i].contains(a) {
                return Err(format!("two resolutions give the column name {a}: FindClusters would overwrite the first; give each resolution once"));
            }
        }
        if self.graph_name.is_empty() || self.graph_name.contains(['\t', '\n', '/']) {
            return Err(format!("graph name {:?}: expected a plain name such as SCT_snn", self.graph_name));
        }
        for &r in &self.resolutions {
            self.optimizer(r).check().map_err(|e| e.to_string())?;
        }
        Ok(())
    }
}

/// Resolutions as R reads them: comma-separated numbers, each parsed by R's own decimal parser
/// (`as.numeric`), so `0.8` is the double R's `resolution = 0.8` gives; a fraction `a/b` is
/// divided in double.
pub fn parse_resolutions(s: &str) -> Result<Vec<f64>, String> {
    s.split(',').map(|t| neighbors::parse_number(t.trim()).map_err(|e| format!("resolution: {e}"))).collect()
}

/// An R integer argument (`n.start`, `n.iter`, `random.seed`, `algorithm`, `modularity.fxn`):
/// a whole number within the range of a C++ `int`, as Rcpp passes it.
pub fn parse_int(name: &str, s: &str) -> Result<i32, String> {
    s.trim().parse::<i32>().map_err(|_| format!("{name}: '{s}' is not an integer between -2147483648 and 2147483647"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_as_r_formats_them() {
        let p = ClusterParams::workflow();
        assert_eq!(p.column_name(0.8), "SCT_snn_res.0.8");
        assert_eq!(ClusterParams::res_name(2.0), "res.2");
        assert_eq!(ClusterParams::res_name(1.2), "res.1.2");
        assert_eq!(ClusterParams::res_name(0.1 + 0.2), "res.0.3");
        assert_eq!(ClusterParams::res_name(1e-5), "res.1e-05");
        assert_eq!(parse_resolutions("0.2, 0.4,0.8,1.2,2.0").unwrap(), vec![0.2, 0.4, 0.8, 1.2, 2.0]);
        assert!(parse_resolutions("0.8,x").is_err());
    }

    #[test]
    fn checks() {
        assert!(ClusterParams::workflow().check().is_ok());
        let p = |f: fn(&mut ClusterParams)| {
            let mut p = ClusterParams::workflow();
            f(&mut p);
            p.check()
        };
        assert!(p(|p| p.algorithm = 4).unwrap_err().contains("Leiden"));
        assert!(p(|p| p.algorithm = 0).is_err());
        assert!(p(|p| p.resolutions = vec![0.8, 0.8]).unwrap_err().contains("overwrite"));
        assert!(p(|p| p.resolutions = vec![0.3, 0.1 + 0.2]).is_err());
        assert!(p(|p| p.resolutions = vec![f64::NAN]).is_err());
        assert_eq!(p(|p| p.n_start = 0).unwrap_err(), "Have to have at least one start");
        assert_eq!(p(|p| {
            p.modularity_fxn = 2;
            p.resolutions = vec![1.5];
        })
        .unwrap_err(), "error: resolution<1 for alternative modularity");
        assert_eq!(parse_int("seed", "-5").unwrap(), -5);
        assert!(parse_int("seed", "1.5").is_err());
    }
}
