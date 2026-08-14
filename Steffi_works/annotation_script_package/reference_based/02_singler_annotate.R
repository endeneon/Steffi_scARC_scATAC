# ============================================================================
# 02_singler_annotate.R — Annotate a query object with SingleR + GSE149614
# ============================================================================
#
# Step 2 of the reference-based annotation pipeline.
#
# Runs SingleR against the Lu et al. 2022 HCC atlas reference built by
# 01_fetch_reference.R. Produces two annotation tiers:
#   - singler_label_cell    : per-cell label
#   - singler_label_cluster : per-cluster majority label (matches the
#                             cluster-level logic of the marker-based route)
#
# Usage:
#   Rscript 02_singler_annotate.R <query.qs2> <reference.rds> <output.qs2> [cluster_col]
#   # or source() and call singler_annotate()
setwd("/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/annotation_script_package/reference_based")
suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(SingleR)
  library(SingleCellExperiment)
})

# Threads granted by the scheduler (exported as NTHREADS in the bsub script).
nthreads <- as.integer(Sys.getenv("NTHREADS", unset = "8"))
if (is.na(nthreads) || nthreads < 1L) nthreads <- 8L

#' Annotate query cells with SingleR
#'
#' @param obj         Query Seurat object (RNA assay will be log-normalized)
#' @param ref         SingleR reference (from build_reference()) or path to .rds
#' @param cluster_col Cluster column for cluster-level labels
#' @return Seurat object with singler_label_cell / singler_label_cluster
singler_annotate <- function(obj, ref, cluster_col = "SCT_snn_res.0.5") {

  if (is.character(ref)) ref <- readRDS(ref)

  DefaultAssay(obj) <- "RNA"
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- NormalizeData(obj, verbose = FALSE)

  query <- GetAssayData(obj, assay = "RNA", layer = "data")

  message(sprintf("Running SingleR (cell level, %d threads) ...", nthreads))
  pred_cell <- SingleR(test = query, ref = ref,
                       labels = if (is(ref, "SingleCellExperiment")) ref$label
                                else colnames(ref),
                       num.threads = nthreads)
  obj$singler_label_cell <- pred_cell$labels

  if (cluster_col %in% colnames(obj@meta.data)) {
    message("Running SingleR (cluster level) ...")
    pred_cl <- SingleR(test = query, ref = ref,
                       labels = if (is(ref, "SingleCellExperiment")) ref$label
                                else colnames(ref),
                       clusters = obj@meta.data[[cluster_col]],
                       num.threads = nthreads)
    cl_labels <- pred_cl$labels[as.character(obj@meta.data[[cluster_col]])]
    obj$singler_label_cluster <- cl_labels
  }

  message("Cell-level SingleR labels:")
  print(sort(table(obj$singler_label_cell), decreasing = TRUE))

  obj
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript 02_singler_annotate.R <query.qs2> <reference.rds> <output.qs2> [cluster_col]")
  }
  obj <- qs_read(args[1], nthreads = nthreads)
  cc  <- if (length(args) >= 4) args[4] else "SCT_snn_res.0.5"
  obj <- singler_annotate(obj, ref = args[2], cluster_col = cc)
  qs_save(obj, args[3], nthreads = nthreads)
  message("Saved: ", args[3])
}
