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

# Which cells vst fitted on in its first step (SCTModel cell.attributes$cells_step1). vst()
# takes every cell when there are at most n_cells (Seurat passes min(ncells, cells)), else
# sample(colnames(umi), n_cells): the first draw from R's generator after SCTransform's
# set.seed(seed.use) (sctransform 0.4.3, R/vst.R). The export records the seed and n_cells, so
# the draw is made again here with R's own sample(), under R's default generator kinds; the
# caller's random-number state is saved and restored.
.cells_step1 <- function(m, cells) {
  sct <- m$summary$result$sctransform
  if (is.null(sct)) sct <- m$summary$result
  n_step1 <- sct$n_cells_step1
  if (is.null(n_step1)) stop("scrustio: the summary does not give vst's number of step-1 cells", call. = FALSE)
  if (n_step1 >= length(cells)) return(rep(TRUE, length(cells)))
  seed <- .sct_parameters(m)$seed_use
  if (is.null(seed)) stop("scrustio: the manifest does not give SCTransform's seed", call. = FALSE)
  chosen <- .with_seed(as.integer(seed), sample(x = cells, size = n_step1))
  cells %in% chosen
}

.with_seed <- function(seed, expr) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old_seed <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had) assign(".Random.seed", old_seed, envir = env) else if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env)
  })
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(seed)
  expr
}

# SCTransform's parameters as the manifest's provenance records them (scrust cluster: under
# "sctransform"; scrust sctransform: at the top).
.sct_parameters <- function(m) {
  p <- m$provenance$parameters
  if (!is.null(p$sctransform)) p$sctransform else p
}

# The SCTModel's slots, as a list (no Seurat needed): feature.attributes, cell.attributes, clips,
# umi.assay, model, arguments, median_umi.
.sct_model_parts <- function(m, tabs, cells) {
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
  cells_step1 <- cells[.cells_step1(m, cells)]

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
.sct_model <- function(m, tabs, cells) {
  p <- .sct_model_parts(m, tabs, cells)
  methods::new(methods::getClass("SCTModel", where = asNamespace("Seurat")),
               feature.attributes = p$feature.attributes, cell.attributes = p$cell.attributes,
               clips = p$clips, umi.assay = p$umi.assay, model = p$model, arguments = p$arguments,
               median_umi = p$median_umi)
}
