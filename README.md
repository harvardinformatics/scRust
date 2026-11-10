# scRust

Rust reimplementations of single-cell RNA-seq tools, aimed at using less memory and running faster.

Each reimplementation is validated against the original tool's outputs, exactly or near-exactly, before it is offered as a replacement in the scRNA-seq preprocessing workflow.

**Planned tools, in order:**
1. DoubletFinder
2. a Seurat-compatible clustering chain
3. scDblFinder
4. emptyDrops
5. decontX
6. SoupX
7. CellBender

**Status:** planning; no code yet.

- Plan and validation framework: [docs/PLAN.md](docs/PLAN.md)
- License: GPL-3.0 ([LICENSE](LICENSE)). Acknowledgements and upstream licenses: [NOTICE](NOTICE)
