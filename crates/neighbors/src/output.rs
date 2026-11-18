//! `scrust neighbors`' output layout (docs/PLAN.md C3), written by the command and by the
//! clustering chain (`scrust cluster`'s `neighbors/`, M11), so that both give the same bytes:
//!
//! - `knn_index.tsv`: barcode, nn_1 .. nn_k: each cell's neighbours, 1-based, nearest first;
//! - `knn_distance.tsv`: barcode, nn_1 .. nn_k: their distances as Annoy returns them (floats,
//!   written with 17 significant digits, so exactly);
//! - `nn.mtx.gz`: `SCT_nn` (Matrix Market, cells x cells, entries in stored order);
//! - `snn.mtx.gz`: `SCT_snn` (Matrix Market, values with 17 significant digits);
//! - with `annoy`, `annoy.ann`: the index in Annoy's file format.

use crate::{Kiss64Random, Metric, NeighborParams, Neighbors};
use scio::{fmt17, Json};
use std::fmt::Write as _;
use std::io;
use std::path::Path;

fn write(path: &Path, text: &str) -> io::Result<()> {
    std::fs::write(path, text).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
}

fn knn_table(barcodes: &[String], k: usize, value: impl Fn(usize, usize) -> String) -> String {
    let mut t = String::from("barcode");
    for r in 1..=k {
        let _ = write!(t, "\tnn_{r}");
    }
    t.push('\n');
    for (i, b) in barcodes.iter().enumerate() {
        t.push_str(b);
        for r in 0..k {
            t.push('\t');
            t.push_str(&value(i, r));
        }
        t.push('\n');
    }
    t
}

/// Writes the layout into `dir` (which must exist). Returns the names written.
pub fn write_outputs(dir: &Path, barcodes: &[String], r: &Neighbors, annoy: bool, threads: usize) -> io::Result<Vec<String>> {
    assert_eq!(barcodes.len(), r.n);
    let k = r.k;
    write(&dir.join("knn_index.tsv"), &knn_table(barcodes, k, |i, j| (r.nn_idx[i * k + j] + 1).to_string()))?;
    write(&dir.join("knn_distance.tsv"), &knn_table(barcodes, k, |i, j| fmt17(r.nn_dist[i * k + j] as f64)))?;
    let version = env!("CARGO_PKG_VERSION");
    let describe = |what: &str| vec![format!(" scrust {version} neighbors: {what}"), format!(" rows and columns are the cells, in the order of knn_index.tsv ({} cells)", r.n)];
    scio::mtx::write_gz(&dir.join("nn.mtx.gz"), &r.nn, &describe("Seurat's SCT_nn graph (FindNeighbors): 1 at (cell, neighbour)"), threads)?;
    scio::mtx::write_gz(&dir.join("snn.mtx.gz"), &r.snn, &describe("Seurat's SCT_snn graph (FindNeighbors, ComputeSNN): s / (k + (k - s)) for s shared neighbours"), threads)?;
    let mut files: Vec<String> = ["knn_index.tsv", "knn_distance.tsv", "nn.mtx.gz", "snn.mtx.gz"].iter().map(|s| s.to_string()).collect();
    if annoy {
        r.index.save(&dir.join("annoy.ann")).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", dir.join("annoy.ann").display())))?;
        files.push("annoy.ann".into());
    }
    Ok(files)
}

/// The parameters as `summary.json` records them.
pub fn params_json(params: &NeighborParams) -> Json {
    let int = |v: usize| Json::Int(v as i64);
    Json::obj(vec![
        ("dims", Json::Vector(params.dims.iter().map(|&d| int(d)).collect())),
        ("k_param", int(params.k_param)),
        ("n_trees", Json::Int(params.n_trees as i64)),
        ("metric", Json::Str(params.metric.seurat_name().into())),
        (
            "annoy_index",
            Json::Str(format!(
                "Annoy{}",
                match params.metric {
                    Metric::Euclidean => "Euclidean",
                    Metric::Angular => "Angular",
                    Metric::Manhattan => "Manhattan",
                }
            )),
        ),
        ("prune_snn", Json::Num(params.prune_snn)),
        ("search_k", Json::Int(params.search_k as i64)),
        ("annoy_seed", Json::Str(Kiss64Random::DEFAULT_SEED.to_string())),
    ])
}

/// The `result` section of `summary.json` (deterministic). `n_components` is the embeddings'
/// columns.
pub fn result_json(r: &Neighbors, n_components: usize) -> Json {
    let int = |v: usize| Json::Int(v as i64);
    let k = r.k;
    let self_first = (0..r.n).filter(|&i| k > 0 && r.nn_idx[i * k] == i as i32).count();
    let self_missing = (0..r.n).filter(|&i| !r.nn_idx[i * k..(i + 1) * k].contains(&(i as i32))).count();
    let mut distinct: Vec<u64> = r.snn.x.iter().map(|v| v.to_bits()).collect();
    distinct.sort_unstable();
    distinct.dedup();
    let st = r.index.stats();
    Json::obj(vec![
        ("n_cells", int(r.n)),
        ("n_components", int(n_components)),
        ("k", int(k)),
        ("annoy_f", int(r.index.f())),
        ("annoy_node_bytes", int(r.index.node_bytes())),
        ("annoy_n_nodes", int(r.index.n_nodes())),
        ("annoy_n_trees", int(r.index.n_trees())),
        ("annoy_file_bytes", int(r.index.n_nodes() * r.index.node_bytes())),
        (
            "annoy_build",
            Json::obj(vec![
                ("split_nodes", int(st.split_nodes)),
                ("leaves", int(st.leaves)),
                ("accepted_at_attempt", Json::Vector(st.accepted_at_attempt.iter().map(|&v| int(v)).collect())),
                ("kept_imbalanced", int(st.kept_imbalanced)),
                ("random_fallbacks", int(st.random_fallbacks)),
                ("random_rounds", int(st.random_rounds)),
                ("zero_margins", int(st.zero_margins)),
                ("draws", Json::Int(st.draws as i64)),
                ("max_depth", int(st.max_depth)),
            ]),
        ),
        ("knn_self_first", int(self_first)),
        ("knn_self_missing", int(self_missing)),
        ("nn_nnz", int(r.nn.nnz())),
        ("snn_nnz", int(r.snn.nnz())),
        ("snn_distinct_values", int(distinct.len())),
    ])
}

/// The manifest's entries for this layout (docs/PLAN.md C3f, M13), in the directory `dir`
/// relative to the output directory (`""` or `"neighbors"`): `n` cells, the graphs' stored
/// entries when known (else read from the files), `annoy.ann` when written. The graphs' rows and
/// columns are the cells of `cells` (a path relative to the output directory and its barcode
/// column), else of this layout's `knn_index.tsv`.
pub fn manifest_entries(dir: &str, n: usize, nnz: Option<(usize, usize)>, annoy: bool, cells: Option<scio::manifest::Ids>) -> Vec<scio::manifest::Entry> {
    use scio::manifest::{Entry, Ids};
    let own = Ids::named("cells", "knn_index.tsv", "barcode");
    let graph_cells = cells.unwrap_or_else(|| Ids::named("cells", &format!("{}knn_index.tsv", if dir.is_empty() { String::new() } else { format!("{}/", dir.trim_end_matches('/')) }), "barcode"));
    let cells = own;
    let graph = |file: &str, role: &str, what: &str, k: Option<usize>| {
        let mut e = Entry::new(file, role, role, "matrix_market", what);
        e.rows = Some(graph_cells.clone());
        e.columns = Some(graph_cells.clone());
        if let Some(k) = k {
            e = e.shape(&[n, n]).nnz(k);
        }
        e
    };
    let mut v = vec![
        Entry::new("knn_index.tsv", "knn_index", "neighbour_indices", "tsv", "each cell's k nearest neighbours by Annoy (FindNeighbors), 1-based cell numbers as in R: barcode, nn_1, ..., nn_k").header(true).rows(cells.clone()),
        Entry::new("knn_distance.tsv", "knn_distance", "neighbour_distances", "tsv", "their distances, Annoy's floats: barcode, nn_1, ..., nn_k").header(true).rows(cells.clone()),
        graph("nn.mtx.gz", "graph_nn", "Seurat's SCT_nn graph: 1 at (cell, neighbour); cells x cells, entries in stored order", nnz.map(|x| x.0)),
        graph("snn.mtx.gz", "graph_snn", "Seurat's SCT_snn graph (ComputeSNN): s / (k + (k - s)) for s shared neighbours, pruned below prune.SNN; cells x cells, entries in stored order", nnz.map(|x| x.1)),
    ];
    if annoy {
        v.push(Entry::new("annoy.ann", "annoy_index", "none", "annoy", "Annoy's index (RcppAnnoy 0.0.23's file format), as save() writes it"));
    }
    // `under` moves the identifier files too: the graphs' are given relative to the output
    // directory already, so they are put back after it
    v.into_iter()
        .map(|e| {
            let graph = e.role.starts_with("graph_");
            let mut e = e.under(dir);
            if graph {
                e.rows = Some(graph_cells.clone());
                e.columns = Some(graph_cells.clone());
            }
            e
        })
        .collect()
}
