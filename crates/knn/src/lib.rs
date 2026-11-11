//! Exact nearest neighbours, memory-bounded and parallel over cells
//! (docs/PLAN.md B3; milestone M3).
//!
//! It will reproduce DoubletFinder's neighbour ranking without building the N x N distance
//! matrix. For each cell, the distances to all cells, computed with the arithmetic of
//! `fields::rdist` (fields 17.3), go into an O(N) buffer per thread; the nearest are then
//! selected by (distance, index), which is exactly R's stable `order()`. Work is split over
//! cells, and results do not depend on thread scheduling.

#[cfg(test)]
mod tests {
    /// Placeholder until M3: shows the crate builds and its tests run.
    #[test]
    fn skeleton() {
        assert_eq!(env!("CARGO_PKG_NAME"), "knn");
    }
}
