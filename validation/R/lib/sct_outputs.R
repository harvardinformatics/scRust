# SCTransform's results in the layout of `scrust sctransform`'s outputs (milestone M4), so the two
# compare byte for byte: cell_attr.tsv, gene_attr.tsv, model_pars.tsv, model_pars_fit.tsv,
# variable_features.tsv (sprintf("%.17g") text through fmt17()), scale_data.npy (npy_write) and
# scale_data.npy.sha256. Needs lib/io.R.

# Captures what Seurat:::SCTransform.default returns, through a passive trace(); call before
# SCTransform() and read env$vst_out after it.
sct_capture <- function() {
  env <- new.env()
  suppressMessages(trace("SCTransform.default", where = asNamespace("Seurat"), print = FALSE,
                         exit = bquote(assign("vst_out", returnValue(), envir = .(env)))))
  env
}

# v: SCTransform.default's result; id_of: feature id by feature name. scale_data = FALSE leaves out
# scale_data.npy and its checksum (M8's light dump level, for objects whose residuals are too large).
sct_write_outputs <- function(v, out, id_of, scale_data = TRUE) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tsv <- function(name, header, cols) {
    lines <- do.call(paste, c(cols, sep = "\t"))
    writeLines(c(paste(header, collapse = "\t"), lines), file.path(out, name), useBytes = TRUE)
  }
  ca <- v$cell_attr
  tsv("cell_attr.tsv", c("barcode", "umi", "log_umi"), list(rownames(ca), fmt17(ca$umi), fmt17(ca$log_umi)))
  ga <- v$gene_attr
  tsv("gene_attr.tsv", c("id", "name", "detection_rate", "gmean", "amean", "variance", "residual_mean",
                         "residual_variance"),
      list(id_of[rownames(ga)], rownames(ga), fmt17(ga$detection_rate), fmt17(ga$gmean), fmt17(ga$amean),
           fmt17(ga$variance), fmt17(ga$residual_mean), fmt17(ga$residual_variance)))
  pars <- function(name, m) {
    tsv(name, c("id", "name", "theta", "(Intercept)", "log_umi"),
        list(id_of[rownames(m)], rownames(m), fmt17(m[, "theta"]), fmt17(m[, "(Intercept)"]),
             fmt17(m[, "log_umi"])))
  }
  pars("model_pars.tsv", v$model_pars)
  pars("model_pars_fit.tsv", v$model_pars_fit)
  vf <- v$variable_features
  tsv("variable_features.tsv", c("id", "name"), list(id_of[vf], vf))
  if (!identical(rownames(v$y), vf)) stop("scale.data rows are not the variable features")
  if (scale_data) {
    npy_write(file.path(out, "scale_data.npy"), unname(v$y))
    writeLines(unname(tools::sha256sum(file.path(out, "scale_data.npy"))), file.path(out, "scale_data.npy.sha256"))
  }
  invisible(out)
}
