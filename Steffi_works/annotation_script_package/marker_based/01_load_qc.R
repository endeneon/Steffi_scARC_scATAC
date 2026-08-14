# ============================================================================
# 01_load_qc.R — Load integrated Seurat object and apply QC filter
# ============================================================================
#
# Step 1 of the marker-based HCC annotation pipeline.
#
# Original analysis context:
#   Input  : merged_integrated_seurat_obj.qs2  (19 HCC tumor samples,
#            Harmony-integrated; assays: SCT + RNA; clustering: SCT_snn_res.0.5)
#   Filter : percent.mt <= 25  ->  223,521 cells -> 165,115 cells (73.9% kept)
#
# Usage:
#   Rscript 01_load_qc.R <input.qs2> <output.qs2> [mt_threshold]
#   # or source() and call load_and_qc()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

#' Load a qs2 Seurat object and apply the mitochondrial QC filter
#'
#' @param input_path   Path to .qs2 Seurat object
#' @param mt_threshold Maximum percent.mt (default 25, as in the original)
#' @param mt_pattern   Regex for mitochondrial genes (default human "^MT-")
#' @return Filtered Seurat object with percent.mt in meta.data
load_and_qc <- function(input_path,
                        mt_threshold = 25,
                        mt_pattern   = "^MT-") {

  message("Loading: ", input_path)
  obj <- qs_read(input_path, nthreads = .nthreads)
  message(sprintf("Loaded: %d cells, %d features", ncol(obj), nrow(obj)))

  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    message("Computing percent.mt with pattern '", mt_pattern, "'")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  }

  n_before <- ncol(obj)
  obj <- subset(obj, subset = percent.mt <= mt_threshold)
  n_after <- ncol(obj)
  message(sprintf("QC filter percent.mt <= %g: %d -> %d cells (%.1f%% retained)",
                  mt_threshold, n_before, n_after, 100 * n_after / n_before))

  obj
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 01_load_qc.R <input.qs2> <output.qs2> [mt_threshold]")
  }
  mt <- if (length(args) >= 3) as.numeric(args[3]) else 25

  obj <- load_and_qc(args[1], mt_threshold = mt)
  qs_save(obj, args[2], nthreads = .nthreads)
  message("Saved: ", args[2])
}
