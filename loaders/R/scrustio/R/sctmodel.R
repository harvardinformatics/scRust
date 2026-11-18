# The SCTModel, rebuilt from the model files as SCTransform makes it: Seurat 5.5.1's
# PrepVSTResults() (R/objects.R) on what sctransform 0.4.3's vst() returns (R/vst.R), with the
# arguments vst() records when Seurat 5.5.1's SCTransform() calls it (R/preprocessing.R:
# SCTransform.default and SCTransform.Assay). The same steps, in the same order, on the same
# values, give the same data frames, column types and row names.

# The model's tables from the output: vst's gene attributes (role sct_feature_attributes), the
# regularised parameters (sct_model_pars_fit), the step-1 fits (sct_model_pars) and the cell
# attributes (sct_cell_attributes), with Seurat's feature names.
.sct_tables <- function(m, feat, to = "seurat") {
  ga <- .read_tsv(.full_path(m, .entry(m, "sct_feature_attributes")),
                  numeric = c("detection_rate", "gmean", "amean", "variance", "residual_mean", "residual_variance"))
  pars <- c("theta", "(Intercept)", "log_umi")
  fit <- .read_tsv(.full_path(m, .entry(m, "sct_model_pars_fit")), numeric = pars)
  step1 <- .read_tsv(.full_path(m, .entry(m, "sct_model_pars")), numeric = pars)
  ca <- .read_tsv(.full_path(m, .entry(m, "sct_cell_attributes")), numeric = c("umi", "log_umi"))
  ga_names <- .map_features(feat, ga$name, ga$id, "model genes (gene_attr)", to)
  if (!identical(fit$name, ga$name)) stop("scrustio: model_pars_fit's genes differ from gene_attr's", call. = FALSE)
  step1_names <- .map_features(feat, step1$name, step1$id, "step-1 genes (model_pars)", to)
  list(gene_attr = ga, gene_names = ga_names, fit = fit, step1 = step1, step1_names = step1_names,
       cell_attr = ca, pars = pars)
}

# Which cells vst fitted on in its first step (SCTModel cell.attributes$cells_step1). Seurat
# 5.5.1's SCTransform.default calls set.seed(seed.use) and passes vst n_cells = min(ncells,
# cells); vst() takes every cell unless n_cells is below the number of cells, else
# sample(x = colnames(umi), size = n_cells), the first draw from R's generator after that
# set.seed() (sctransform 0.4.3, R/vst.R: 216-218 and 233; Seurat's PrepVSTResults then marks the
# drawn cells TRUE). The export does not list the cells: it records the seed and ncells (the
# manifest's parameters and seeds). So the draw is made again here, from those, on the cells in
# counts/ order (vst's colnames(umi)), with R's own sample() under R's default generator kinds
# set for the draw; the caller's generator kinds and .Random.seed are restored afterwards (and
# .Random.seed removed again if there was none).
#
# The self-check, given the raw counts: vst's next draw, from the same stream, picks the step-1
# genes (vst.R: 225-266) among the model's genes detected in at least min_cells of the drawn
# cells and overdispersed, with probabilities from the density of their log10 geometric means.
# The export lists those genes in that order (sct_model_pars), so the gene draw is made too, with
# sctransform's own row_gmean() and row_var(), and must give that list; a different set of cells
# would not. It needs no file the export does not have.
#
# The gene draw is weighted, so a last-bit change in its inputs can change it: in sctransform's
# internal functions, in R's density() (R 4.4.0 changed its values by about 0.1%) or in the C
# library's log() and exp(). So the check runs, and refuses an output whose list differs, only
# with the versions it mirrors and was validated with (`tested`: R 4.5.3, sctransform 0.4.3, on
# Linux). Otherwise, as without sctransform, the cells are drawn without it (the draw needs only
# R's sample()) and one warning says why.
.cells_step1 <- function(m, cells, counts = NULL, tabs = NULL, tested = .STEP1_TESTED) {
  s <- .step1_settings(m, length(cells))
  if (s$n_cells >= length(cells)) return(rep(TRUE, length(cells)))
  if (is.null(s$seed)) stop("scrustio: the manifest does not give SCTransform's seed", call. = FALSE)
  st <- NULL
  if (!is.null(counts) && !is.null(tabs)) {
    # (loading sctransform's namespace may not touch the caller's generator either)
    why <- .keep_rng(.step1_check_skip(tested))
    if (is.null(why)) {
      st <- .keep_rng(.step1_gene_stats(counts, tabs$gene_names))
    } else {
      warning("scrustio: the SCTModel's cells_step1 was drawn again without its check against the output's step-1 genes: ",
              why, call. = FALSE)
    }
  }
  if (!is.null(st) && is.null(s$min_cells)) stop("scrustio: the manifest does not give vst's min_cells", call. = FALSE)
  drawn <- .with_seed(s$seed, {
    chosen <- sample(x = cells, size = s$n_cells)
    list(cells = chosen, genes = if (is.null(st)) NULL else .step1_genes(st, counts, chosen, s))
  })
  if (!is.null(st) && !identical(unname(drawn$genes), unname(tabs$step1_names))) {
    stop("scrustio: SCTransform's step-1 cells, drawn again from the manifest's seed (", s$seed, ") and ncells (",
         s$n_cells, " of ", length(cells), " cells), do not give the step-1 genes the output lists in vst's sampling ",
         "order (sct_model_pars): the SCTModel's cells_step1 cannot be rebuilt", call. = FALSE)
  }
  cells %in% drawn$cells
}

# The step-1 settings from the manifest: n_cells = min(ncells, cells) from its parameters, which
# the summary's n_cells_step1 (where given) must equal; the seed from its parameters, which its
# seeds (where given) must equal; n_genes (0 or absent: every gene) and min_cells.
.step1_settings <- function(m, n_total) {
  par <- .sct_parameters(m)
  one <- function(x) if (is.null(x) || length(x) == 0) NULL else as.numeric(unlist(x))[1]
  ncells <- one(par$ncells)
  if (is.null(ncells)) stop("scrustio: the manifest does not give SCTransform's ncells", call. = FALSE)
  n_cells <- min(ncells, n_total)
  res <- m$summary$result
  sct <- if (!is.null(res$sctransform)) res$sctransform else res
  n_step1 <- one(sct$n_cells_step1)
  if (!is.null(n_step1) && n_step1 != n_cells) {
    stop("scrustio: the summary's n_cells_step1 (", n_step1, ") is not min(ncells, cells) = ", n_cells,
         " from the manifest's parameters", call. = FALSE)
  }
  seed <- one(par$seed_use)
  recorded <- one(m$provenance$seeds[["SCTransform seed.use"]])
  if (!is.null(seed) && !is.null(recorded) && seed != recorded) {
    stop("scrustio: the manifest's SCTransform seeds differ (parameters ", seed, ", seeds ", recorded, ")", call. = FALSE)
  }
  if (is.null(seed)) seed <- recorded
  n_genes <- one(par$n_genes)
  if (!is.null(n_genes) && n_genes == 0) n_genes <- NULL
  list(n_cells = n_cells, seed = seed, n_genes = n_genes, min_cells = one(par$min_cells))
}

# The versions the self-check mirrors and was validated with (M14: K11, K13, M17).
.STEP1_TESTED <- c(R = "4.5.3", sctransform = "0.4.3", os = "linux")

# Why the self-check cannot run in this session (NULL when it can): sctransform is not installed,
# or R, sctransform or the operating system is not the one in `tested`.
.step1_check_skip <- function(tested = .STEP1_TESTED) {
  if (!requireNamespace("sctransform", quietly = TRUE)) return("sctransform is not installed")
  have <- c(R = as.character(getRversion()), sctransform = as.character(utils::packageVersion("sctransform")),
            os = tolower(Sys.info()[["sysname"]]))
  got <- unname(have[names(tested)])
  differ <- names(tested)[is.na(got) | got != tested]
  if (length(differ) == 0) return(NULL)
  paste0("the check is made only under ", paste(names(tested), tested, collapse = ", "),
         " (what it mirrors and was validated with), and this session has ",
         paste(differ, have[differ], collapse = ", "))
}

# vst()'s per-gene statistics over all cells (vst.R: 199-207, 242-246): log10 of the geometric
# mean (row_gmean, eps 1) and whether the variance exceeds the mean (row_var, rowMeans), for the
# model's genes, which are the raw counts' genes vst kept, in their order. sctransform's own
# functions, row by row, so on the whole matrix they give vst's values on umi[genes, ]. Called
# only when .step1_check_skip() finds nothing.
.step1_gene_stats <- function(counts, genes) {
  idx <- match(genes, rownames(counts))
  if (anyNA(idx) || is.unsorted(idx, strictly = TRUE)) {
    stop("scrustio: the model's genes (gene_attr) are not the raw counts' genes in their order", call. = FALSE)
  }
  row_gmean <- utils::getFromNamespace("row_gmean", "sctransform")
  row_var <- utils::getFromNamespace("row_var", "sctransform")
  list(idx = idx, genes = genes, log_gmean = log10(row_gmean(counts, eps = 1))[idx],
       overdispersed = (row_var(counts) - Matrix::rowMeans(counts))[idx] > 0)
}

# vst()'s step-1 genes for the drawn cells (vst.R: 225-226, 242-257, 262-266), drawing from the
# generator as vst does; the caller has set the seed and made the cells' draw.
.step1_genes <- function(st, counts, chosen, s) {
  det <- Matrix::rowSums(counts[st$idx, match(chosen, colnames(counts)), drop = FALSE] > 0)
  keep <- det >= s$min_cells & st$overdispersed
  genes_step1 <- st$genes[keep]
  log_gmean_step1 <- st$log_gmean[keep]
  if (!is.null(s$n_genes) && s$n_genes < length(genes_step1)) {
    dens <- stats::density(x = log_gmean_step1, bw = "nrd", adjust = 1)
    prob <- 1 / (stats::approx(x = dens$x, y = dens$y, xout = log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = s$n_genes, prob = prob)
  }
  genes_step1
}

# Evaluates expr, then puts the caller's generator kinds and .Random.seed back (and removes
# .Random.seed again if there was none).
.keep_rng <- function(expr) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old_seed <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    if (had) assign(".Random.seed", old_seed, envir = env) else if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env)
  })
  expr
}

# Evaluates expr after set.seed(seed) under R's default generator kinds, whatever the caller's.
.with_seed <- function(seed, expr) {
  .keep_rng({
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
    set.seed(as.integer(seed))
    expr
  })
}

# SCTransform's parameters as the manifest's provenance records them (scrust cluster: under
# "sctransform"; scrust sctransform: at the top).
.sct_parameters <- function(m) {
  p <- m$provenance$parameters
  if (!is.null(p$sctransform)) p$sctransform else p
}

# The SCTModel's slots, as a list (no Seurat needed): feature.attributes, cell.attributes, clips,
# umi.assay, model, arguments, median_umi. `counts`, the raw counts (when the output holds them),
# lets .cells_step1() check its draw.
.sct_model_parts <- function(m, tabs, cells, counts = NULL) {
  res <- m$summary$result
  sct <- if (!is.null(res$sctransform)) res$sctransform else res
  par <- .sct_parameters(m)
  num <- function(x, what) {
    if (is.null(x)) stop("scrustio: the summary does not give ", what, call. = FALSE)
    as.numeric(unlist(x))
  }
  n_cells_total <- length(cells)
  # vst()'s result, as far as PrepVSTResults reads it
  cell_attr <- data.frame(umi = tabs$cell_attr$umi, log_umi = tabs$cell_attr$log_umi)
  rownames(cell_attr) <- cells
  gene_attr <- tabs$gene_attr[, c("detection_rate", "gmean", "amean", "variance", "residual_mean", "residual_variance")]
  rownames(gene_attr) <- tabs$gene_names
  model_pars_fit <- as.matrix(tabs$fit[, tabs$pars])
  dimnames(model_pars_fit) <- list(tabs$gene_names, tabs$pars)
  model_pars <- as.matrix(tabs$step1[, tabs$pars])
  dimnames(model_pars) <- list(tabs$step1_names, tabs$pars)
  cells_step1 <- cells[.cells_step1(m, cells, counts, tabs)]

  # PrepVSTResults (Seurat 5.5.1, R/objects.R)
  cell.names <- intersect(x = cells, y = rownames(x = cell_attr))
  cell.cols <- c("umi", "gene", "log_umi", "log_gene", "umi_per_gene", "log_umi_per_gene")
  cell.cols <- intersect(x = cell.cols, y = colnames(x = cell_attr))
  cell.attrs <- cell_attr[cell.names, cell.cols, drop = FALSE]
  colnames(x = cell.attrs) <- gsub(pattern = "gene", replacement = "feature", x = colnames(x = cell.attrs))
  cell.attrs[, "cells_step1"] <- FALSE
  cells_step1 <- intersect(x = cells_step1, y = rownames(x = cell.attrs))
  cell.attrs[cells_step1, "cells_step1"] <- TRUE
  feature.cols <- c("detection_rate", "gmean", "variance", "residual_mean", "residual_variance")
  feature.cols <- intersect(x = feature.cols, y = colnames(x = gene_attr))
  feature.attrs <- gene_attr[, feature.cols, drop = FALSE]
  feature.attrs <- cbind(feature.attrs, model_pars_fit[rownames(feature.attrs), , drop = FALSE])
  feature.attrs[, "genes_log_gmean_step1"] <- FALSE
  genes_step1 <- intersect(x = rownames(model_pars), y = rownames(feature.attrs))
  feature.attrs[genes_step1, "genes_log_gmean_step1"] <- TRUE
  feature.attrs[, paste0("step1_", colnames(model_pars))] <- NA
  feature.attrs[genes_step1, paste0("step1_", colnames(model_pars))] <- model_pars[genes_step1, ]

  res_clip <- num(sct$res_clip_range, "vst's res_clip_range")
  sct_clip <- num(sct$clip_range, "SCTransform's clip.range")
  arguments <- list(
    glmGamPoi_check = TRUE,
    latent_var = "log_umi",
    batch_var = NULL,
    latent_var_nonreg = NULL,
    n_genes = num(par$n_genes, "vst's n_genes"),
    n_cells = as.numeric(min(num(par$ncells, "SCTransform's ncells"), n_cells_total)),
    method = "glmGamPoi_offset",
    do_regularize = TRUE,
    theta_regularization = "od_factor",
    res_clip_range = res_clip,
    bin_size = num(par$bin_size, "vst's bin_size"),
    min_cells = num(par$min_cells, "vst's min_cells"),
    residual_type = "pearson",
    return_cell_attr = TRUE,
    return_gene_attr = TRUE,
    return_corrected_umi = TRUE,
    min_variance = "umi_median",
    bw_adjust = num(par$bw_adjust, "vst's bw_adjust"),
    gmean_eps = 1,
    theta_estimation_fun = "theta.ml",
    theta_given = NULL,
    exclude_poisson = TRUE,
    use_geometric_mean = TRUE,
    use_geometric_mean_offset = FALSE,
    fix_intercept = FALSE,
    fix_slope = FALSE,
    scale_factor = NA,
    vst.flavor = "v2",
    verbosity = 0,
    set_min_var = num(sct$min_variance, "vst's set_min_var"),
    sct.clip.range = sct_clip
  )
  clips <- list(vst = arguments$res_clip_range, sct = arguments$sct.clip.range)
  median_umi <- NA
  if ("scale_factor" %in% names(x = arguments)) median_umi <- arguments$scale_factor
  if (is.na(x = median_umi)) {
    if ("umi" %in% colnames(x = cell.attrs)) {
      median_umi <- median(cell.attrs$umi)
    } else if ("log_umi" %in% colnames(x = cell.attrs)) {
      median_umi <- median(10^cell.attrs$log_umi)
    }
  }
  list(feature.attributes = feature.attrs, cell.attributes = cell.attrs, clips = clips,
       umi.assay = "RNA", model = "y ~ log_umi", arguments = arguments, median_umi = median_umi)
}

# The SCTModel object (Seurat's class), as PrepVSTResults() makes it.
.sct_model <- function(m, tabs, cells, counts = NULL) {
  p <- .sct_model_parts(m, tabs, cells, counts)
  methods::new(methods::getClass("SCTModel", where = asNamespace("Seurat")),
               feature.attributes = p$feature.attributes, cell.attributes = p$cell.attributes,
               clips = p$clips, umi.assay = p$umi.assay, model = p$model, arguments = p$arguments,
               median_umi = p$median_umi)
}
