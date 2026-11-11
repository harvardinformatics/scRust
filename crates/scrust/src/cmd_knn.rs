//! `scrust knn`: exact nearest neighbours as DoubletFinder ranks them (`order()` of
//! `fields::rdist` distances), the shared stage on its own (docs/PLAN.md A1 and B3).

use crate::cli::{self, Args, Phases};
use knn::{neighbours, Digest, Order, Points};
use scio::npy;
use scio::{Json, Npy};
use std::path::PathBuf;

pub const USAGE: &str = "\
Usage: scrust knn (--embeddings <file.npy> | --from-stages <dir> --context <ctx>)
                  --queries <n> --neighbours <L> --output count|digest|positions --out <file.npy>
                  [--dims <d>] [--threads N]

For each of the first n points (the queries; DoubletFinder's real cells), the first L positions
of R's order() of its Euclidean distances to all points, computed as fields::rdist does, ties to
the lower index. Points are the rows of a column-major (Fortran order) matrix of doubles, e.g.
an instrumented R run's S3_embeddings block for one context.

  --output count       int32 per query: points with index >= n among positions 2..L
                       (pANN's numerator); no sort is needed
  --output digest      S4 digest (h1, h2) per query, as the instrumented R runs dump it
  --output positions   the L x n matrix of 1-based positions (int32)
";

pub fn main(argv: &[String]) -> Result<(), String> {
    let args = Args::parse(argv)?;
    args.only(&["embeddings", "from-stages", "context", "queries", "neighbours", "output", "out", "dims", "threads"])?;
    let mut phases = Phases::new();
    let threads = cli::threads(&args)?;
    let (data, n, dim) = if let Some(f) = args.get("embeddings") {
        let a = Npy::read(&PathBuf::from(f)).map_err(|e| e.to_string())?;
        (a.block_f64(0, a.nrow()), a.nrow(), a.ncol())
    } else {
        let dir = args.path("from-stages")?;
        let ctx = args.required("context")?;
        let st = doubletfinder::stages::StageDir::open(&dir).map_err(|e| e.to_string())?;
        let b = st.block("S3_embeddings.npy", ctx).ok_or_else(|| format!("no S3_embeddings block for context {ctx}"))?.clone();
        let a = st.read("S3_embeddings.npy").map_err(|e| e.to_string())?;
        (a.block_f64(b.row_start, b.nrow), b.nrow, b.ncol)
    };
    let dims: usize = args.number("dims", dim)?;
    if dims == 0 || dims > dim {
        return Err(format!("--dims {dims}: the points have {dim} columns"));
    }
    let queries: usize = args.number("queries", 0)?;
    let l: usize = args.number("neighbours", 0)?;
    if queries == 0 || queries > n || l == 0 {
        return Err(format!("need 0 < --queries <= {n} and --neighbours > 0"));
    }
    let out = args.path("out")?;
    phases.mark("load");
    let points = Points::new(&data[..dims * n], n, dims);
    match args.required("output")? {
        "count" => {
            let c = neighbours(&points, queries, l, Order::FirstThenSet, threads, |_, pos| {
                pos[1..].iter().filter(|&&p| p as usize >= queries).count() as i32
            });
            phases.mark("knn");
            npy::write_i32(&out, &c, &[queries]).map_err(|e| e.to_string())?;
        }
        "digest" => {
            let dg = Digest::new(l.min(n));
            let d = neighbours(&points, queries, l, Order::Sorted, threads, |_, pos| dg.of(pos));
            phases.mark("knn");
            let mut v: Vec<f64> = d.iter().map(|x| x.0 as f64).collect();
            v.extend(d.iter().map(|x| x.1 as f64));
            npy::write_f64(&out, &v, &[queries, 2]).map_err(|e| e.to_string())?;
        }
        "positions" => {
            let p = neighbours(&points, queries, l, Order::Sorted, threads, |_, pos| pos.iter().map(|&x| x as i32 + 1).collect::<Vec<i32>>());
            phases.mark("knn");
            let lv = l.min(n);
            let flat: Vec<i32> = p.into_iter().flatten().collect();
            npy::write_i32(&out, &flat, &[lv, queries]).map_err(|e| e.to_string())?;
        }
        o => return Err(format!("--output: unknown kind '{o}'")),
    }
    phases.mark("write");
    let info = Json::obj(vec![
        ("points", Json::Int(n as i64)),
        ("dims", Json::Int(dims as i64)),
        ("queries", Json::Int(queries as i64)),
        ("neighbours", Json::Int(l as i64)),
        ("threads", Json::Int(threads as i64)),
        ("phases", phases.to_json()),
    ]);
    eprint!("{}", scio::json::to_file_text(&info));
    Ok(())
}
