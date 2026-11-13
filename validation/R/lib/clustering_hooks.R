# Stage dumps for instrumented_clustering.R (docs/PLAN.md C4 and C5; milestone M8). Requires
# lib/io.R, lib/sct_outputs.R, lib/clustering.R and lib/doubletfinder_hooks.R, whose machinery
# (.dfi, dfi_add, dfi_rng, dfi_scalar, dfi_finish: index.tsv, rng.npy/rng.tsv, stages.json) and
# SCTransform and RunPCA hooks (stages S2a-S2d and S3, M2) are used unchanged, in one context,
# "chain". The mechanism is the same as M2's: trace() adds passive hooks to the installed functions
# that read values and write files; they draw nothing from R's random stream and change no
# argument or local variable, so the instrumented run's final files equal the direct run's byte
# for byte (checked by compare_clustering.R).
#
# Stages, in <run>/stages/ (index conventions as M2's: genes are 0-based rows of the input matrix,
# cells 0-based columns, in input order):
#   S2a-S2d  SCTransform at seed 12345 (lib/doubletfinder_hooks.R), plus vst/: SCTransform.default's
#            result in `scrust sctransform`'s layout (lib/sct_outputs.R; scale_data.npy except at
#            level light)
#   S3       RunPCA (50 PCs): embeddings, loadings, sdev, irlba's d; RNG at irlba's entry and exit
#   K2       RNG at RunUMAP's entry and exit
#   K3_items         the items Annoy holds, getItemsVector(i) for every cell (cells x dims, the
#                    floats as doubles); stages.json says whether they equal R's own float rounding
#                    of the embeddings (writeBin(size = 4))
#   annoy.ann        K4: the index of a second build from the same rows (new index, addItem in
#                    cell order, build(n.trees)), saved with save(). save() reloads the file and
#                    re-derives the roots, so the first build, which FindNeighbors searches, is left
#                    alone; the reloaded index is searched for every cell as well, and stages.json
#                    says whether it returns the same neighbours and distances
#   K5_nn_index, K5_nn_dist  every cell's k neighbours (1-based, as R holds them) and distances
#                    (Annoy's floats as doubles)
#   K7_optimizer_ids RunModularityClustering's result before GroupSingletons (0-based, numbered by
#                    size); K7_final_ids GroupSingletons' result; RNG at both functions' entry and
#                    exit (GroupSingletons calls set.seed(1) when there are singletons)
#   K8_silhouette    approxSilhouette: width, the nearest other cluster (0-based level index);
#                    K8_centroids (clusters x dims) and K8_cluster_var, in level order
#   K9_*             neighborPurity: dist50 (each cell's distance to its 50th nearest other cell,
#                    findDistance() on neighborPurity's own index, the same call), the radius (their
#                    median, stages.json), neighbor_count (cells within it), at level full or with
#                    --purity-lists yes the lists (neighbors_ptr, neighbors_idx: 1-based, in
#                    findNeighbors' order), weights, sums (cells x clusters), totals, purity and
#                    maximum (0-based level index)
#   K10              sctransform::correct(), which makes the corrected counts inside vst(): the
#                    median log_umi it regresses to (stages.json), its genes, and how many of vst's
#                    residuals it reads lie beyond vst's clip range (vst clips only afterwards)

.cli <- new.env(parent = emptyenv())

# No RNG checkpoint at SCTransform's entry: in a fresh session R's stream does not exist yet, and
# SCTransform.Seurat's first statement reseeds it (the vst_entry checkpoint follows the reseeding).
cli_sct_seurat_entry <- function(env) {
  dfi_count("SCTransform.Seurat")
}

cli_sct_default_exit <- function(env, rv) {
  dfi_sct_default_exit(env, rv)
  id_of <- setNames(.dfi$features, .dfi$features)
  sct_write_outputs(rv, file.path(.dfi$dir, "vst"), id_of, scale_data = .dfi$level != "light")
  dfi_scalar(vst_umi_corrected_nnz = length(rv$umi_corrected@x))
}

cli_correct_entry <- function(env) {
  dfi_count("correct")
  x <- env$x
  y <- x$y
  clip <- x$arguments$res_clip_range
  beyond <- 0
  for (s in seq(1, nrow(y), by = 500)) {
    b <- y[s:min(nrow(y), s + 499), , drop = FALSE]
    beyond <- beyond + sum(b < clip[1] | b > clip[2])
  }
  dfi_add("K10_genes", gene_index(rownames(y), .dfi$features), "K10")
  dfi_scalar(correct_n_genes = nrow(y), correct_n_cells = ncol(y), correct_res_clip_range = clip,
             correct_residuals_beyond_clip = beyond, correct_bin_size = x$arguments$bin_size,
             correct_y_is_data = identical(env$data, "y"))
  dfi_rng("K10", "correct_entry")
}

cli_correct_exit <- function(env, rv) {
  lv <- env$x$arguments$latent_var
  # apply() may leave a one-column matrix in the column: compare values, not attributes
  med <- unique(as.vector(env$cell_attr[, lv]))
  if (length(med) != 1) stop("cli: correct() did not set ", lv, " to one value")
  dfi_scalar(correct_latent_var = lv, correct_median_log_umi = med,
             correct_median_is_R_median = identical(med, median(as.vector(env$x$cell_attr[, lv]))),
             correct_regressor_row = env$regressor_data[1, ],
             correct_scale_factor = env$scale_factor, correct_do_round = env$do_round,
             correct_do_pos = env$do_pos)
}

cli_umap_entry <- function(env) {
  dfi_count("RunUMAP.Seurat")
  dfi_rng("K2", "umap_entry")
}
cli_umap_exit <- function(env, rv) dfi_rng("K2", "umap_exit")

cli_new_index <- function(metric, f) {
  switch(metric,
         euclidean = new(RcppAnnoy::AnnoyEuclidean, f),
         cosine = new(RcppAnnoy::AnnoyAngular, f),
         manhattan = new(RcppAnnoy::AnnoyManhattan, f),
         stop("cli: metric ", metric))
}

cli_annoy_build_exit <- function(env, rv) {
  dfi_count("AnnoyBuildIndex")
  dfi_rng("K3", "annoy_build_exit")
  data <- env$data
  n <- nrow(data)
  f <- ncol(data)
  items <- matrix(0, n, f)
  for (i in seq_len(n)) items[i, ] <- rv$getItemsVector(i - 1L)
  as_float <- readBin(writeBin(as.vector(data), raw(), size = 4), "double", n = n * f, size = 4)
  dim(as_float) <- dim(data)
  dfi_add("K3_items", items, "K3")
  b <- cli_new_index(env$metric, f)
  for (ii in seq_len(n)) b$addItem(ii - 1, data[ii, ])
  b$build(env$n.trees)
  path <- file.path(.dfi$dir, "annoy.ann")
  saved <- b$save(path)
  node_bytes <- 16 + 4 * f
  dfi_scalar(annoy_metric = env$metric, annoy_f = f, annoy_n_trees_arg = env$n.trees,
             annoy_n_items = rv$getNItems(), annoy_n_trees = rv$getNTrees(),
             annoy_items_are_float_of_embeddings = identical(unname(items), unname(as_float)),
             annoy_saved = saved, annoy_file_bytes = file.size(path), annoy_node_bytes = node_bytes,
             annoy_n_nodes = file.size(path) / node_bytes, annoy_reloaded_n_items = b$getNItems(),
             annoy_reloaded_n_trees = b$getNTrees())
  .cli$reloaded <- b
}

cli_annoy_search_exit <- function(env, rv) {
  dfi_count("AnnoySearch")
  idx <- rv$nn.idx
  if (any(idx != round(idx))) stop("cli: non-integer neighbour ids")
  storage.mode(idx) <- "integer"
  n <- nrow(idx)
  dfi_add("K5_nn_index", idx, "K5")
  dfi_add("K5_nn_dist", rv$nn.dists, "K5")
  self_any <- vapply(seq_len(n), function(i) any(idx[i, ] == i), logical(1))
  dfi_scalar(knn_k = env$k, knn_search_k = env$search.k, knn_self_first = sum(idx[, 1] == seq_len(n)),
             knn_self_missing = sum(!self_any))
  b <- .cli$reloaded
  same <- NA
  if (!is.null(b)) {
    q <- env$query
    same <- TRUE
    for (x in seq_len(n)) {
      r <- b$getNNsByVectorList(q[x, ], env$k, env$search.k, TRUE)
      if (!identical(as.integer(r$item + 1), idx[x, ]) || !identical(r$distance, rv$nn.dists[x, ])) {
        same <- FALSE
        break
      }
    }
  }
  dfi_scalar(annoy_reloaded_search_identical = same)
  dfi_rng("K5", "annoy_search_exit")
}

cli_rmc_entry <- function(env) {
  dfi_count("RunModularityClustering")
  dfi_rng("K7", "optimizer_entry")
}

cli_rmc_exit <- function(env, rv) {
  dfi_add("K7_optimizer_ids", as.integer(rv), "K7")
  dfi_scalar(optimizer_modularity_function = env$modularity, optimizer_resolution = env$resolution,
             optimizer_algorithm = env$algorithm, optimizer_n_start = env$n.start,
             optimizer_n_iter = env$n.iter, optimizer_random_seed = env$random.seed,
             optimizer_communities = length(unique(rv)))
  dfi_rng("K7", "optimizer_exit")
}

cli_gs_entry <- function(env) {
  dfi_count("GroupSingletons")
  dfi_rng("K7", "groupsingletons_entry")
}

cli_gs_exit <- function(env, rv) {
  ids <- as.character(rv)
  if (any(!grepl("^[0-9]+$", ids))) stop("cli: non-numeric labels after GroupSingletons")
  dfi_add("K7_final_ids", as.integer(ids), "K7")
  dfi_scalar(singletons = I(as.character(env$singletons)), n_singletons = length(env$singletons),
             final_ids_are_character = is.character(rv))
  dfi_rng("K7", "groupsingletons_exit")
}

cli_sil_exit <- function(env, rv) {
  dfi_count("approxSilhouette")
  uc <- as.character(env$uclust)
  dfi_add("K8_silhouette_width", rv$width, "K8")
  dfi_add("K8_silhouette_other", match(as.character(rv$other), uc) - 1L, "K8")
  dfi_add("K8_centroids", do.call(rbind, env$averaged), "K8")
  dfi_add("K8_cluster_var", env$clust.var, "K8")
  dfi_scalar(silhouette_levels = I(uc), silhouette_dims = ncol(env$x))
}

cli_purity_exit <- function(env, rv) {
  dfi_count("neighborPurity")
  d50 <- BiocNeighbors::findDistance(env$idx, k = env$k, num.threads = env$num.threads, BNPARAM = env$BNPARAM)
  radius <- env$dist
  nout <- env$nout
  cnt <- lengths(nout)
  uc <- as.character(env$uclust)
  dfi_add("K9_dist50", d50, "K9")
  dfi_add("K9_neighbor_count", as.integer(cnt), "K9")
  if (.cli$purity_lists) {
    dfi_add("K9_neighbors_ptr", as.integer(c(0, cumsum(cnt))), "K9")
    dfi_add("K9_neighbors_idx", as.integer(unlist(nout, use.names = FALSE)), "K9")
  }
  dfi_add("K9_weights", env$w, "K9")
  dfi_add("K9_sums", t(env$aggregated), "K9")
  dfi_add("K9_totals", env$totals, "K9")
  dfi_add("K9_purity", rv$purity, "K9")
  dfi_add("K9_maximum", match(as.character(rv$maximum), uc) - 1L, "K9")
  dfi_scalar(purity_k = env$k, purity_radius = radius, purity_radius_is_median_of_dist50 = identical(median(d50), radius),
             purity_index_class = class(env$idx)[1], purity_levels = I(uc),
             purity_neighbor_lists = .cli$purity_lists, purity_mean_neighbors = mean(cnt))
}

cli_install <- function(stage_dir, level, cells, features, purity_lists) {
  dfi_init(stage_dir, level, cells, features)
  .dfi$ctx <- "chain"
  .dfi$ctxs <- "chain"
  .cli$purity_lists <- purity_lists
  .cli$reloaded <- NULL
  q <- function(fun) bquote(.(as.name(fun))(environment()))
  qx <- function(fun) bquote(.(as.name(fun))(environment(), returnValue()))
  ns <- asNamespace
  suppressMessages({
    trace("SCTransform.Seurat", where = ns("Seurat"), print = FALSE,
          tracer = q("cli_sct_seurat_entry"), exit = qx("dfi_sct_seurat_exit"))
    trace("SCTransform.default", where = ns("Seurat"), print = FALSE, exit = qx("cli_sct_default_exit"))
    vst_at <- c(2L, dfi_step(sctransform::vst, "n_genes < length(genes_step1)"))
    trace("vst", where = ns("sctransform"), print = FALSE, at = vst_at,
          tracer = q("dfi_vst_step"), exit = qx("dfi_vst_exit"))
    trace("correct", where = ns("sctransform"), print = FALSE,
          tracer = q("cli_correct_entry"), exit = qx("cli_correct_exit"))
    trace("get_model_pars", where = ns("sctransform"), print = FALSE,
          tracer = q("dfi_gmp_entry"), exit = qx("dfi_gmp_exit"))
    trace("reg_model_pars", where = ns("sctransform"), print = FALSE, exit = qx("dfi_reg_exit"))
    trace("future_lapply", where = ns("future.apply"), print = FALSE, tracer = q("dfi_future_lapply_entry"))
    trace("estimate_dispersions_roughly", where = ns("glmGamPoi"), print = FALSE, exit = qx("dfi_glm_disp_init"))
    trace("estimate_betas_roughly_group_wise", where = ns("glmGamPoi"), print = FALSE, exit = qx("dfi_glm_beta_init"))
    trace("estimate_betas_group_wise", where = ns("glmGamPoi"), print = FALSE, exit = qx("dfi_glm_beta"))
    trace("glm_gp", where = ns("glmGamPoi"), print = FALSE, exit = qx("dfi_glm_exit"))
    trace("RunPCA.default", where = ns("Seurat"), print = FALSE,
          tracer = q("dfi_pca_entry"), exit = qx("dfi_pca_exit"))
    trace("irlba", where = ns("irlba"), print = FALSE, tracer = q("dfi_irlba_entry"), exit = qx("dfi_irlba_exit"))
    trace("RunUMAP.Seurat", where = ns("Seurat"), print = FALSE,
          tracer = q("cli_umap_entry"), exit = qx("cli_umap_exit"))
    trace("AnnoyBuildIndex", where = ns("Seurat"), print = FALSE, exit = qx("cli_annoy_build_exit"))
    trace("AnnoySearch", where = ns("Seurat"), print = FALSE, exit = qx("cli_annoy_search_exit"))
    trace("RunModularityClustering", where = ns("Seurat"), print = FALSE,
          tracer = q("cli_rmc_entry"), exit = qx("cli_rmc_exit"))
    trace("GroupSingletons", where = ns("Seurat"), print = FALSE,
          tracer = q("cli_gs_entry"), exit = qx("cli_gs_exit"))
    trace("approxSilhouette", where = ns("bluster"), print = FALSE, exit = qx("cli_sil_exit"))
    trace("neighborPurity", where = ns("bluster"), print = FALSE, exit = qx("cli_purity_exit"))
  })
  invisible(TRUE)
}

cli_finish <- function(seu, umap) {
  dfi_rng("K9", "chain_exit")
  .cli$reloaded <- NULL
  dfi_finish(expected = list(
    SCTransform.Seurat = 1L, SCTransform.default = 1L, vst = 1L, vst_before_gene_sampling = 1L,
    get_model_pars = 1L, reg_model_pars = 1L, correct = 1L, RunPCA.default = 1L, irlba = 1L,
    RunUMAP.Seurat = if (umap) 1L else 0L, AnnoyBuildIndex = 1L, AnnoySearch = 1L,
    RunModularityClustering = 1L, GroupSingletons = 1L, approxSilhouette = 1L, neighborPurity = 1L,
    glm_gp = .dfi$counts$fit_bins))
}
