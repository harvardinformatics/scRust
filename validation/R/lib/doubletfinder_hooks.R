# Stage dumps for instrumented_doubletfinder.R (docs/PLAN.md B5, stages S1-S7). Requires lib/io.R.
#
# Mechanism: base R's trace() inserts passive hooks into the installed functions, in their
# namespaces: at entry, at exit (on.exit, reading the frame's locals and returnValue()) and, for
# sctransform::vst, before one statement. trace() also updates the copies imported by other
# packages (Seurat imports vst, irlba and future_lapply; DoubletFinder imports bkde and Seurat's
# generics; S3 methods are updated in the registration table), so the call path is the
# installed code, unchanged. The hooks only read values and write files: they draw no random
# numbers, change no argument or local variable, and run no BLAS. The proof is in the outputs:
# the instrumented run's cells.tsv, sweep_stats.tsv, bcmvn.tsv and summary.json "result" (which
# includes a digest of the final .Random.seed) must equal the direct run's byte for byte.
#
# trace() in a namespace does not touch the copy exported to the attached package
# (package:DoubletFinder), so the instrumented script calls DoubletFinder::paramSweep etc.,
# which are the traced namespace functions.
#
# Output, in <run>/stages/:
#   index.tsv          one row per (array, context): file, ctx, row_start, nrow, ncol, dtype,
#                      stage. Arrays of one quantity from all contexts are stacked along rows
#                      ("ragged"); row_start/nrow give each context's block. Large matrices have
#                      one file per context (<name>.<ctx>.npy) and row_start 0.
#   rng.npy, rng.tsv   .Random.seed checkpoints (int32, one row of 626 per checkpoint) and labels
#   stages.json        scalars per context (sizes, min_variance, bandwidths, irlba iterations,
#                      checks) and hook call counts
# Contexts: sweep_pN<pN> for the six paramSweep iterations, summary (summarizeSweep), final
# (doubletFinder). Index conventions: genes are 0-based rows of the input matrix
# (features.tsv); cells are 0-based columns of the matrix SCTransform sees: the real cells (all
# input cells, or the sweep subsample) followed by doublets X1..Xn; doublet pairs are 0-based
# positions in the real cells; S4 positions are R's 1-based order() positions.
#
# Dump levels: "full" (tier T): everything, including the clipped residual matrix y (genes x
# cells) and the S4 order positions; "standard" (small and medium P): no y and no S4
# positions (their per-cell digests remain); scale.data (the PCA input) is kept; "light"
# (larger P): no residual-sized matrices at all.

.dfi <- new.env(parent = emptyenv())

dfi_init <- function(stage_dir, level, input_barcodes, features) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  .dfi$dir <- stage_dir
  .dfi$level <- match.arg(level, c("full", "standard", "light"))
  .dfi$input_barcodes <- input_barcodes
  .dfi$features <- features
  .dfi$ctx <- "setup"
  .dfi$ctxs <- character()
  .dfi$arrays <- list()
  .dfi$index <- list()
  .dfi$rng <- list()
  .dfi$scalars <- list()
  .dfi$counts <- list()
  .dfi$in_gmp <- FALSE
  .dfi$in_summary <- FALSE
  .dfi$bin <- 0L
  .dfi$glm_rows <- list()
  .dfi$glm_cur <- list()
  .dfi$sweep_real <- NULL
  invisible(TRUE)
}

# Bookkeeping ---------------------------------------------------------------------------------

dfi_count <- function(name) {
  .dfi$counts[[name]] <- (if (is.null(.dfi$counts[[name]])) 0L else .dfi$counts[[name]]) + 1L
}

dfi_rng <- function(stage, label) {
  s <- rng_state()
  if (is.null(s)) s <- rep(NA_integer_, 626)
  .dfi$rng[[length(.dfi$rng) + 1]] <- list(ctx = .dfi$ctx, stage = stage, label = label, seed = s)
}

dfi_scalar <- function(...) {
  l <- list(...)
  cur <- .dfi$scalars[[.dfi$ctx]]
  if (is.null(cur)) cur <- list()
  cur[names(l)] <- l
  .dfi$scalars[[.dfi$ctx]] <- cur
}

# Appends one context's block of a stacked ("ragged") array.
dfi_add <- function(name, value, stage) {
  if (is.logical(value)) storage.mode(value) <- "integer"
  dimnames(value) <- NULL
  names(value) <- NULL
  .dfi$arrays[[name]] <- c(.dfi$arrays[[name]], list(list(ctx = .dfi$ctx, value = value, stage = stage)))
}

# Writes one large matrix per context right away.
dfi_write_big <- function(name, value, stage) {
  if (is.logical(value)) storage.mode(value) <- "integer"
  if (is.double(value) && all(value == round(value)) && max(abs(value)) < .Machine$integer.max &&
      grepl("positions", name)) storage.mode(value) <- "integer"
  dimnames(value) <- NULL
  f <- sprintf("%s.%s.npy", name, .dfi$ctx)
  npy_write(file.path(.dfi$dir, f), value)
  .dfi$index[[length(.dfi$index) + 1]] <- list(
    file = f, ctx = .dfi$ctx, row_start = 0L, nrow = NROW(value), ncol = NCOL(value),
    dtype = if (is.integer(value)) "<i4" else "<f8", stage = stage)
}

gene_index <- function(names, universe) {
  i <- match(names, universe)
  if (anyNA(i)) stop("dfi: gene names not found")
  i - 1L
}

# S4 digest of order positions (one column per real cell; R's 1-based positions), exact in
# double arithmetic: h1 = sum_r r * pos[r]; h2 = (sum_r (pos[r] * w[r]) mod P) mod P with
# P = 2^31 - 1, w[r] = 16807^r mod P. All intermediate values are integers below 2^53.
dfi_weights <- function(L) {
  P <- 2147483647
  w <- numeric(L)
  v <- 1
  for (r in seq_len(L)) {
    v <- (v * 16807) %% P
    w[r] <- v
  }
  w
}

dfi_digest_matrix <- function(M) {
  P <- 2147483647
  L <- nrow(M)
  r <- as.double(seq_len(L))
  w <- dfi_weights(L)
  D <- as.double(M)
  dim(D) <- dim(M)
  h1 <- colSums(D * r)
  h2 <- colSums((D * w) %% P) %% P
  cbind(h1, h2)
}

# Hooks (called from trace()) ------------------------------------------------------------------

dfi_paramSweep_entry <- function(env) {
  dfi_count("paramSweep")
  .dfi$ctx <- "sweep"
  dfi_rng("S1", "paramSweep_entry")
}

dfi_pps_entry <- function(env) {
  dfi_count("parallel_paramSweep")
  .dfi$ctx <- paste0("sweep_pN", env$pN[env$n])
  .dfi$ctxs <- c(.dfi$ctxs, .dfi$ctx)
  dfi_rng("S1", "before_doublets")
  idx <- match(env$real.cells, .dfi$input_barcodes)
  if (anyNA(idx)) stop("dfi: sweep cells not in the input")
  if (is.null(.dfi$sweep_real)) {
    .dfi$sweep_real <- idx
    if (length(idx) < length(.dfi$input_barcodes)) {
      npy_write(file.path(.dfi$dir, "S1_subsample.npy"), idx - 1L)
      .dfi$index[[length(.dfi$index) + 1]] <- list(file = "S1_subsample.npy", ctx = "sweep",
        row_start = 0L, nrow = length(idx), ncol = 1L, dtype = "<i4", stage = "S1")
    }
  } else if (!identical(.dfi$sweep_real, idx)) {
    stop("dfi: sweep cells differ between pN values")
  }
  dfi_scalar(pN = env$pN[env$n], n_real = env$n.real.cells, subsampled = length(idx) < length(.dfi$input_barcodes))
}

dfi_pps_exit <- function(env) {
  rc <- env$real.cells
  dfi_add("S1_doublets", cbind(match(env$real.cells1, rc) - 1L, match(env$real.cells2, rc) - 1L), "S1")
  M <- env$dist.mat
  storage.mode(M) <- "integer"
  dfi_add("S4_digest", dfi_digest_matrix(M), "S4")
  if (.dfi$level == "full") dfi_write_big("S4_positions", M, "S4")
  P <- vapply(env$sweep.res.list, function(d) d$pANN, numeric(env$n.real.cells))
  dfi_add("S5_pann_sweep", P, "S5")
  dfi_scalar(n_doublets = env$n_doublets, nCells = env$nCells, ind = env$ind,
             pK = env$pK, pk_temp = round(env$nCells * env$pK),
             first_position_is_self = sum(M[1, ] == seq_len(ncol(M))))
  dfi_rng("S1", "sweep_iteration_exit")
}

dfi_sct_seurat_entry <- function(env) {
  dfi_count("SCTransform.Seurat")
  dfi_rng("S1", "after_doublets")
}

dfi_sct_seurat_exit <- function(env, rv) {
  vf <- VariableFeatures(rv)
  dfi_add("S2d_variable_features", gene_index(vf, .dfi$features), "S2d")
  dfi_scalar(variable_features_equal_top_features = identical(vf, .dfi$top_features))
  if (.dfi$level != "light") {
    sd <- GetAssayData(rv, assay = "SCT", layer = "scale.data")
    dfi_scalar(scale_data_equals_sctransform_default_y =
                 identical(sd[rownames(.dfi$y_default), , drop = FALSE], .dfi$y_default, num.eq = FALSE))
    rm(sd)
  }
  .dfi$y_default <- NULL
  dfi_rng("S2", "sctransform_exit")
}

dfi_sct_default_exit <- function(env, rv) {
  dfi_count("SCTransform.default")
  top <- rv$variable_features
  .dfi$top_features <- top
  dfi_add("S2d_top_features", gene_index(top, .dfi$features), "S2d")
  dfi_scalar(scale_data_rows_are_top_features = identical(rownames(rv$y), top),
             clip_range = sqrt(ncol(rv$y) / 30))
  if (.dfi$level != "light") .dfi$y_default <- rv$y
}

dfi_vst_step <- function(env) {
  if (!exists("genes_step1", envir = env, inherits = FALSE)) {
    dfi_count("vst")
    .dfi$vst_features <- rownames(env$umi)
    .dfi$vst_cells <- colnames(env$umi)
    dfi_scalar(vst_rows_are_input_features = identical(.dfi$vst_features, .dfi$features),
               vst_n_cells = ncol(env$umi))
    dfi_rng("S2a", "vst_entry")
  } else {
    dfi_count("vst_before_gene_sampling")
    dfi_rng("S2a", "after_cell_sampling")
    dfi_add("S2a_genes_step1_pool", gene_index(env$genes_step1, env$genes), "S2a")
  }
}

dfi_vst_exit <- function(env, rv) {
  genes <- env$genes
  dfi_add("S2a_genes", gene_index(genes, .dfi$vst_features), "S2a")
  ca <- rv$cell_attr
  if (!identical(rownames(ca), .dfi$vst_cells)) stop("dfi: cell_attr rows differ from the vst input")
  dfi_add("S2a_cells", cbind(ca$umi, ca$log_umi), "S2a")
  ga <- rv$gene_attr
  if (!identical(rownames(ga), genes)) stop("dfi: gene_attr rows differ from genes")
  dfi_add("S2a_gene_stats", cbind(env$genes_log_gmean, env$genes_amean, env$genes_var,
                                  ga$detection_rate, ga$gmean), "S2a")
  dfi_add("S2a_cells_step1", match(env$cells_step1, .dfi$vst_cells) - 1L, "S2a")
  dfi_add("S2a_genes_step1", gene_index(env$genes_step1, genes), "S2a")
  dfi_add("S2a_log_gmean_step1", env$genes_log_gmean_step1, "S2a")
  sampled <- exists("sampling_prob", envir = env, inherits = FALSE)
  if (sampled) {
    dfi_add("S2a_sampling_prob", env$sampling_prob, "S2a")
    d <- env$log_gmean_dens
    dfi_add("S2a_density", cbind(d$x, d$y), "S2a")
    dfi_scalar(density_bw = d$bw, density_n = d$n)
  }
  mp <- rv$model_pars
  if (!identical(rownames(mp), env$genes_step1)) stop("dfi: model_pars rows differ from genes_step1")
  dfi_add("S2b_model_pars", mp[, c("theta", "(Intercept)", "log_umi"), drop = FALSE], "S2b")
  mpf <- rv$model_pars_fit
  if (!identical(rownames(mpf), genes)) stop("dfi: model_pars_fit rows differ from genes")
  dfi_add("S2c_model_pars_fit", mpf[, c("theta", "(Intercept)", "log_umi"), drop = FALSE], "S2c")
  dfi_add("S2c_outliers", as.integer(rv$model_pars_outliers), "S2c")
  dfi_add("S2d_residual_stats", cbind(ga$residual_mean, ga$residual_variance), "S2d")
  if (.dfi$level == "full") dfi_write_big("S2d_y", rv$y, "S2d")
  dfi_scalar(min_variance = rv$arguments$set_min_var, n_genes = length(genes),
             n_genes_step1 = length(env$genes_step1), n_cells_step1 = length(env$cells_step1),
             n_genes_step1_pool = if (sampled) length(env$sampling_prob) else length(env$genes_step1),
             gene_sampling = sampled, n_cells = nrow(ca),
             res_clip_range = rv$arguments$res_clip_range, n_bins_fit = .dfi$bin)
  dfi_rng("S2", "vst_exit")
}

dfi_gmp_entry <- function(env) {
  dfi_count("get_model_pars")
  .dfi$in_gmp <- TRUE
  .dfi$bin <- 0L
  .dfi$glm_rows <- list()
  .dfi$glm_cur <- list()
  dfi_rng("S2b", "fits_entry")
}

dfi_gmp_exit <- function(env, rv) {
  .dfi$in_gmp <- FALSE
  dfi_add("S2b_model_pars_raw", rv[, c("theta", "(Intercept)", "log_umi"), drop = FALSE], "S2b")
  dfi_add("S2b_glm", do.call(rbind, .dfi$glm_rows), "S2b")
  dfi_scalar(n_theta_below_1e7 = sum(rv[, "theta"] < 1e-7), n_theta_inf = sum(!is.finite(rv[, "theta"])))
  dfi_rng("S2b", "fits_exit")
}

dfi_future_lapply_entry <- function(env) {
  if (!.dfi$in_gmp) return(invisible())
  .dfi$bin <- .dfi$bin + 1L
  dfi_count("fit_bins")
  dfi_rng("S2b", sprintf("bin%d_entry", .dfi$bin))
}

dfi_glm_disp_init <- function(env, rv) if (.dfi$in_gmp) .dfi$glm_cur$disp_init <- rv
dfi_glm_beta_init <- function(env, rv) if (.dfi$in_gmp) .dfi$glm_cur$beta_init <- rv[, 1]
dfi_glm_beta <- function(env, rv) {
  if (!.dfi$in_gmp) return(invisible())
  if (is.null(.dfi$glm_cur$beta1)) {
    .dfi$glm_cur$beta1 <- rv$Beta[, 1]
    .dfi$glm_cur$iter1 <- rv$iterations
  } else {
    .dfi$glm_cur$iter2 <- rv$iterations
    .dfi$glm_cur$dev2 <- rv$deviances
  }
}
dfi_glm_exit <- function(env, rv) {
  if (!.dfi$in_gmp) return(invisible())
  dfi_count("glm_gp")
  g <- .dfi$glm_cur
  trend <- rv$overdispersion_shrinkage_list$dispersion_trend
  if (is.null(trend)) trend <- rep(NA_real_, nrow(rv$Beta))
  .dfi$glm_rows <- c(.dfi$glm_rows, list(cbind(
    disp_init = g$disp_init, beta_init = g$beta_init, beta1 = g$beta1, iter1 = g$iter1,
    disp_mle = rv$overdispersions, disp_trend = trend, beta2 = rv$Beta[, 1],
    iter2 = g$iter2, deviance2 = g$dev2)))
  .dfi$glm_cur <- list()
}

dfi_reg_exit <- function(env, rv) {
  dfi_count("reg_model_pars")
  dfi_add("S2c_poisson_genes", gene_index(env$all_poisson_genes, env$genes), "S2c")
  dfi_add("S2c_smoothing_genes", gene_index(rownames(env$model_pars), env$genes), "S2c")
  dfi_scalar(reg_bw = env$bw, mean_cell_sum = env$mean_cell_sum,
             n_outliers = sum(env$outliers), n_all_poisson_genes = length(env$all_poisson_genes),
             n_smoothing_genes = nrow(env$model_pars))
}

dfi_pca_entry <- function(env) {
  dfi_count("RunPCA.default")
  x <- env$object
  dfi_add("S2d_pca_features", gene_index(rownames(x), .dfi$features), "S2d")
  dfi_scalar(pca_n_features = nrow(x), pca_cells_are_vst_cells = identical(colnames(x), .dfi$vst_cells))
  if (.dfi$level != "light") {
    dfi_write_big("S2d_scale_data", x, "S2d")
  }
  dfi_rng("S3", "pca_entry")
}

dfi_pca_exit <- function(env, rv) {
  dfi_add("S3_embeddings", rv@cell.embeddings, "S3")
  dfi_add("S3_loadings", rv@feature.loadings, "S3")
  dfi_add("S3_sdev", matrix(rv@stdev, nrow = 1), "S3")
  dfi_scalar(total_variance = rv@misc$total.variance)
  dfi_rng("S3", "pca_exit")
}

dfi_irlba_entry <- function(env) {
  dfi_count("irlba")
  dfi_rng("S3", "irlba_entry")
}

dfi_irlba_exit <- function(env, rv) {
  dfi_add("S3_irlba_d", matrix(rv$d, nrow = 1), "S3")
  dfi_scalar(irlba_iter = rv$iter, irlba_mprod = rv$mprod)
  dfi_rng("S3", "irlba_exit")
}

dfi_df_entry <- function(env) {
  dfi_count("doubletFinder")
  .dfi$ctx <- "final"
  .dfi$ctxs <- c(.dfi$ctxs, "final")
  dfi_rng("S1", "before_doublets")
}

dfi_df_exit <- function(env, rv) {
  rc <- env$real.cells
  dfi_add("S1_doublets", cbind(match(env$real.cells1, rc) - 1L, match(env$real.cells2, rc) - 1L), "S1")
  n <- env$n_real.cells
  k <- env$k
  L <- k + 1
  D <- env$dist.mat
  P <- 2147483647
  r <- as.double(seq_len(L))
  w <- dfi_weights(L)
  h <- matrix(0, n, 2)
  full <- .dfi$level == "full"
  if (full) pos_all <- matrix(0L, L, n)
  self <- 0L
  for (i in seq_len(n)) {
    pos <- order(D[, i])[seq_len(L)]
    pd <- as.double(pos)
    h[i, 1] <- sum(pd * r)
    h[i, 2] <- sum((pd * w) %% P) %% P
    if (pos[1] == i) self <- self + 1L
    if (full) pos_all[, i] <- pos
  }
  dfi_add("S4_digest", h, "S4")
  if (full) dfi_write_big("S4_positions", pos_all, "S4")
  dfi_add("S5_pann_final", env$pANN$pANN, "S5")
  dfi_add("S7_doublet", as.integer(env$classifications == "Doublet"), "S7")
  dfi_scalar(n_doublets = env$n_doublets, nCells = env$nCells, k = k, pK = env$pK, nExp = env$nExp,
             first_position_is_self = self)
  dfi_rng("S7", "doubletFinder_exit")
}

dfi_summary_entry <- function(env) {
  dfi_count("summarizeSweep")
  .dfi$ctx <- "summary"
  .dfi$in_summary <- TRUE
  .dfi$sweep_names <- names(env$sweep.list)
  .dfi$sweep_item <- 0L
}
dfi_summary_exit <- function(env, rv) .dfi$in_summary <- FALSE

dfi_bkde_exit <- function(env, rv) {
  if (!.dfi$in_summary) return(invisible())
  dfi_count("bkde")
  # one block per (pN, pK), in sweep order, labelled with the sweep list name
  .dfi$sweep_item <- .dfi$sweep_item + 1L
  .dfi$ctx <- paste0("summary:", .dfi$sweep_names[.dfi$sweep_item])
  dfi_add("S6_bkde", cbind(rv$x, rv$y), "S6")
}

dfi_bc_entry <- function(env) {
  if (!.dfi$in_summary) return(invisible())
  dfi_count("bimodality_coefficient")
  dfi_add("S6_gkde", env$x, "S6")
}

# Installation ---------------------------------------------------------------------------------

# Index of the statement in a function body (as.list(body(f)), as trace(at = ) counts) whose
# text contains `pattern`; it must be unique.
dfi_step <- function(f, pattern) {
  b <- as.list(body(f))
  hit <- which(vapply(b, function(e) any(grepl(pattern, deparse(e), fixed = TRUE)), logical(1)))
  hit <- hit[hit > 1]
  if (length(hit) != 1) stop("dfi: statement '", pattern, "' found ", length(hit), " times")
  hit
}

dfi_install <- function() {
  q <- function(fun) bquote(.(as.name(fun))(environment()))
  qx <- function(fun) bquote(.(as.name(fun))(environment(), returnValue()))
  ns <- asNamespace
  suppressMessages({
    trace("paramSweep", where = ns("DoubletFinder"), print = FALSE, tracer = q("dfi_paramSweep_entry"))
    trace("parallel_paramSweep", where = ns("DoubletFinder"), print = FALSE,
          tracer = q("dfi_pps_entry"), exit = q("dfi_pps_exit"))
    trace("doubletFinder", where = ns("DoubletFinder"), print = FALSE,
          tracer = q("dfi_df_entry"), exit = q("dfi_df_exit"))
    trace("summarizeSweep", where = ns("DoubletFinder"), print = FALSE,
          tracer = q("dfi_summary_entry"), exit = qx("dfi_summary_exit"))
    trace("bimodality_coefficient", where = ns("DoubletFinder"), print = FALSE, tracer = q("dfi_bc_entry"))
    trace("bkde", where = ns("KernSmooth"), print = FALSE, exit = qx("dfi_bkde_exit"))
    trace("SCTransform.Seurat", where = ns("Seurat"), print = FALSE,
          tracer = q("dfi_sct_seurat_entry"), exit = qx("dfi_sct_seurat_exit"))
    trace("SCTransform.default", where = ns("Seurat"), print = FALSE, exit = qx("dfi_sct_default_exit"))
    vst_at <- c(2L, dfi_step(sctransform::vst, "n_genes < length(genes_step1)"))
    trace("vst", where = ns("sctransform"), print = FALSE, at = vst_at,
          tracer = q("dfi_vst_step"), exit = qx("dfi_vst_exit"))
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
  })
  invisible(TRUE)
}

# Writes the stacked arrays, RNG checkpoints, scalars and the index.
dfi_finish <- function(expected) {
  for (name in names(.dfi$arrays)) {
    blocks <- .dfi$arrays[[name]]
    vals <- lapply(blocks, `[[`, "value")
    is_mat <- vapply(vals, function(v) !is.null(dim(v)), logical(1))
    if (any(is_mat) && !all(is_mat)) stop("dfi: mixed vector/matrix blocks for ", name)
    nrows <- vapply(vals, NROW, integer(1))
    ncols <- vapply(vals, NCOL, integer(1))
    if (length(unique(ncols)) != 1) stop("dfi: blocks of ", name, " differ in columns")
    all_int <- all(vapply(vals, is.integer, logical(1)))
    value <- if (all(is_mat)) do.call(rbind, vals) else unlist(vals, use.names = FALSE)
    if (all_int) storage.mode(value) <- "integer" else storage.mode(value) <- "double"
    f <- paste0(name, ".npy")
    npy_write(file.path(.dfi$dir, f), value)
    starts <- cumsum(c(0L, nrows))[seq_along(nrows)]
    for (b in seq_along(blocks)) {
      .dfi$index[[length(.dfi$index) + 1]] <- list(
        file = f, ctx = blocks[[b]]$ctx, row_start = starts[b], nrow = nrows[b], ncol = ncols[b],
        dtype = if (all_int) "<i4" else "<f8", stage = blocks[[b]]$stage)
    }
  }
  idx <- do.call(rbind, lapply(.dfi$index, as.data.frame, stringsAsFactors = FALSE))
  idx <- idx[order(idx$file, idx$row_start), ]
  tsv_write(idx, file.path(.dfi$dir, "index.tsv"))
  seeds <- do.call(rbind, lapply(.dfi$rng, `[[`, "seed"))
  npy_write(file.path(.dfi$dir, "rng.npy"), seeds, dtype = "<i4")
  tsv_write(data.frame(row = seq_len(nrow(seeds)) - 1L,
                       ctx = vapply(.dfi$rng, `[[`, "", "ctx"),
                       stage = vapply(.dfi$rng, `[[`, "", "stage"),
                       label = vapply(.dfi$rng, `[[`, "", "label")),
            file.path(.dfi$dir, "rng.tsv"))
  counts <- .dfi$counts
  problems <- character()
  for (nm in names(expected)) {
    got <- if (is.null(counts[[nm]])) 0L else counts[[nm]]
    if (!identical(as.integer(got), as.integer(expected[[nm]]))) {
      problems <- c(problems, sprintf("%s called %d times, expected %d", nm, got, expected[[nm]]))
    }
  }
  json_write(list(level = .dfi$level, contexts = .dfi$ctxs, scalars = .dfi$scalars,
                  hook_counts = counts, expected_counts = expected, problems = I(problems)),
             file.path(.dfi$dir, "stages.json"))
  if (length(problems)) stop("dfi: ", paste(problems, collapse = "; "))
  invisible(TRUE)
}
