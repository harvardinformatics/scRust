//! FindClusters' result as plain files (docs/PLAN.md A3 and C3, "Outputs"), and the `result`
//! section of `summary.json`.
//!
//! - `seurat_clusters.tsv` (`barcode`, `seurat_clusters`) and `seurat_clusters_levels.txt`: the
//!   layout `scrust doubletfinder --clusters` reads;
//! - `clusters.tsv`: one column per resolution, named as Seurat names it (`SCT_snn_res.0.8`), and
//!   each column's levels in `<column>_levels.txt`;
//! - `cluster_ids.txt` and `nclusters.txt`, as the workflow's `write_cluster_metadata` writes
//!   them (the identities' levels, and their number);
//! - `metadata.tsv`: the cells' metadata columns the clustering makes (each `res.` column and
//!   `seurat_clusters`), in cell order;
//! - `optimizer_ids.tsv`: each resolution's optimizer result before GroupSingletons (0-based,
//!   numbered by size), the ids M8's instrumented runs record as `K7_optimizer_ids`;
//! - with the starts kept, `stages/<column>_partitions.npy`: every start's own partition
//!   (starts × cells, int32), as the instrumented optimizer prints them.
//!
//! Every file is the same whatever the thread count.

use crate::find_clusters::FindClusters;
use scio::Json;
use std::fmt::Write as _;
use std::io;
use std::path::Path;

fn write(path: &Path, text: &str) -> io::Result<()> {
    std::fs::write(path, text).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))
}

/// Lines, each ending in `\n` (R's `writeLines`).
fn lines(v: &[String]) -> String {
    let mut s = String::new();
    for l in v {
        s.push_str(l);
        s.push('\n');
    }
    s
}

/// A table with a `barcode` column and the given columns of text.
fn table(barcodes: &[String], header: &[&str], cols: &[Vec<&str>]) -> String {
    let mut s = String::from("barcode");
    for h in header {
        s.push('\t');
        s.push_str(h);
    }
    s.push('\n');
    for (i, b) in barcodes.iter().enumerate() {
        s.push_str(b);
        for c in cols {
            s.push('\t');
            s.push_str(c[i]);
        }
        s.push('\n');
    }
    s
}

/// Writes the label files into `dir` (which must exist). Returns the names written.
pub fn write_labels(dir: &Path, barcodes: &[String], r: &FindClusters) -> io::Result<Vec<String>> {
    assert_eq!(barcodes.len(), r.n_cells);
    let mut written = Vec::new();
    let mut put = |name: &str, text: &str| -> io::Result<()> {
        write(&dir.join(name), text)?;
        written.push(name.to_string());
        Ok(())
    };
    let sc = r.seurat_clusters.labels();
    put("seurat_clusters.tsv", &table(barcodes, &["seurat_clusters"], std::slice::from_ref(&sc)))?;
    put("seurat_clusters_levels.txt", &lines(&r.seurat_clusters.levels))?;
    let names: Vec<&str> = r.resolutions.iter().map(|x| x.column.name.as_str()).collect();
    let cols: Vec<Vec<&str>> = r.resolutions.iter().map(|x| x.column.labels()).collect();
    put("clusters.tsv", &table(barcodes, &names, &cols))?;
    for x in &r.resolutions {
        put(&format!("{}_levels.txt", x.column.name), &lines(&x.column.levels))?;
    }
    // write_cluster_metadata: levels(Idents), without NA or empty ones, and their number
    let ids: Vec<String> = r.cluster_ids().iter().filter(|s| !s.is_empty()).cloned().collect();
    put("cluster_ids.txt", &lines(&ids))?;
    put("nclusters.txt", &format!("{}\n", ids.len()))?;
    let mut meta_names = names.clone();
    meta_names.push("seurat_clusters");
    let mut meta_cols = cols.clone();
    meta_cols.push(sc);
    put("metadata.tsv", &table(barcodes, &meta_names, &meta_cols))?;
    let opt: Vec<Vec<String>> = r.resolutions.iter().map(|x| x.optimizer.ids.iter().map(|v| v.to_string()).collect()).collect();
    let opt_ref: Vec<Vec<&str>> = opt.iter().map(|c| c.iter().map(String::as_str).collect()).collect();
    put("optimizer_ids.tsv", &table(barcodes, &names, &opt_ref))?;
    Ok(written)
}

/// Writes every start's partition, `<dir>/<column>_partitions.npy` (starts × cells, int32), for
/// the resolutions whose partitions were kept. Returns the names written.
pub fn write_partitions(dir: &Path, r: &FindClusters) -> io::Result<Vec<String>> {
    let mut written = Vec::new();
    for x in &r.resolutions {
        let starts = &x.optimizer.starts;
        if starts.iter().any(|s| s.partition.is_none()) {
            continue;
        }
        let n = r.n_cells;
        let mut data = Vec::with_capacity(starts.len() * n);
        // Fortran order (the npy writer's layout): column by column
        for c in 0..n {
            for s in starts {
                data.push(s.partition.as_ref().unwrap()[c] as i32);
            }
        }
        let name = format!("{}_partitions.npy", x.column.name);
        scio::npy::write_i32(&dir.join(&name), &data, &[starts.len(), n])?;
        written.push(name);
    }
    Ok(written)
}

fn int(v: usize) -> Json {
    Json::Int(v as i64)
}

/// The `result` section of `summary.json`: deterministic (no timings), so equal between runs
/// that agree.
pub fn result_json(r: &FindClusters) -> Json {
    let per = r
        .resolutions
        .iter()
        .map(|x| {
            let o = &x.optimizer;
            let g = &x.grouped;
            let mut printed = x.printed();
            printed.pop(); // the elapsed time, which is in the run section
            Json::obj(vec![
                ("resolution", Json::Num(x.resolution)),
                ("name", Json::Str(x.res_name.clone())),
                ("column", Json::Str(x.column.name.clone())),
                ("resolution_used", Json::Num(o.resolution2)),
                ("start_iterations", Json::Vector(o.starts.iter().map(|s| int(s.iterations)).collect())),
                ("start_modularity", Json::Vector(o.starts.iter().map(|s| Json::Num(s.modularity)).collect())),
                ("best_start", int(o.best_start + 1)),
                ("max_modularity", Json::Num(o.max_modularity)),
                ("max_modularity_printed", Json::Str(o.max_modularity_printed())),
                ("communities", int(o.n_communities)),
                ("next_int_draws", Json::Int(o.draws as i64)),
                ("next_int_power_of_two_draws", Json::Int(o.pow2_draws as i64)),
                ("singletons", Json::Vector(g.singletons.iter().map(|v| Json::Str(v.to_string())).collect())),
                ("singleton_cells", Json::Vector(g.singleton_cells.iter().map(|&c| int(c + 1)).collect())),
                ("singleton_ties", Json::Vector(g.n_tied.iter().map(|&t| int(t)).collect())),
                ("singleton_connectivity", Json::nums(&g.max_connectivity)),
                ("singleton_joined", Json::Vector(g.joined.iter().map(|v| Json::Str(v.to_string())).collect())),
                ("n_singletons", int(g.singletons.len())),
                ("final_clusters", int(g.n_final())),
                ("levels", Json::strs(&x.column.levels)),
                ("printed", Json::strs(&printed)),
                ("message", g.message().map(Json::Str).unwrap_or(Json::Null)),
                ("rng_after_groupsingletons", g.rng_after.as_ref().map(|v| Json::Vector(v.iter().map(|&s| Json::Int(s as i64)).collect())).unwrap_or(Json::Null)),
            ])
        })
        .collect();
    Json::obj(vec![
        ("n_cells", int(r.n_cells)),
        ("graph", Json::Str(r.params.graph_name.clone())),
        ("nodes", int(r.n_nodes)),
        ("edges", int(r.n_edges)),
        ("total_edge_weight", Json::Num(r.total_edge_weight)),
        ("columns", Json::Vector(r.resolutions.iter().map(|x| Json::Str(x.column.name.clone())).collect())),
        ("resolutions", Json::List(per)),
        ("seurat_clusters_levels", Json::strs(&r.seurat_clusters.levels)),
        ("n_clusters", int(r.seurat_clusters.levels.len())),
        ("rng_after", r.rng_after().map(|v| Json::Vector(v.iter().map(|&s| Json::Int(s as i64)).collect())).unwrap_or(Json::Null)),
    ])
}

/// The parameters as `summary.json` records them (`run.parameters`).
pub fn params_json(p: &crate::ClusterParams) -> Json {
    Json::obj(vec![
        ("graph_name", Json::Str(p.graph_name.clone())),
        ("modularity_fxn", Json::Int(p.modularity_fxn as i64)),
        ("resolution", Json::nums(&p.resolutions)),
        ("algorithm", Json::Int(p.algorithm as i64)),
        ("n_start", Json::Int(p.n_start as i64)),
        ("n_iter", Json::Int(p.n_iter as i64)),
        ("random_seed", Json::Int(p.random_seed as i64)),
        ("group_singletons", Json::Bool(p.group_singletons)),
    ])
}

/// One line per resolution for a log: the optimizer's numbers and GroupSingletons'.
pub fn describe(r: &FindClusters) -> String {
    let mut s = String::new();
    for x in &r.resolutions {
        let o = &x.optimizer;
        let _ = writeln!(
            s,
            "{}: {} nodes, {} edges; maximum modularity {} (start {} of {}); {} communities; {} singletons; {} clusters",
            x.column.name,
            o.n_nodes,
            o.n_edges,
            o.max_modularity_printed(),
            o.best_start + 1,
            o.starts.len(),
            o.n_communities,
            x.grouped.singletons.len(),
            x.column.levels.len()
        );
    }
    s
}
