# ============================================================================
# 02_marker_de.R — Per-cluster marker DE (presto Wilcoxon on RNA assay)
# ============================================================================
#
# Step 2 of the marker-based HCC annotation pipeline.
#
# Original analysis context:
#   - RNA assay layers joined (JoinLayers) after the SCT-based integration
#   - Log-normalized RNA data (NormalizeData) used for DE, NOT the SCT assay
#   - Wilcoxon rank-sum via presto::wilcoxauc (fast at 165k-cell scale)
#   - Clustering column: SCT_snn_res.0.5 (16 clusters)
#   - Outputs: cluster_markers_res0.5.csv (all), cluster_markers_top10.csv
#
# Usage:
#   Rscript 02_marker_de.R <input.qs2> <out_prefix> [cluster_col]
#   # or source() and call run_marker_de()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(presto)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

#' Run per-cluster Wilcoxon marker DE on the RNA assay
#'
#' @param obj         Seurat object (post-QC)
#' @param cluster_col Metadata column with cluster IDs (default SCT_snn_res.0.5)
#' @param out_prefix  File prefix for CSV outputs (NULL = skip writing)
#' @return data.frame of markers (presto wilcoxauc output, sorted)
run_marker_de <- function(obj,
                          cluster_col = "cluster_final",
                          out_prefix  = "cluster_markers_res0.5") {

  stopifnot(cluster_col %in% colnames(obj@meta.data))

  DefaultAssay(obj) <- "RNA"
  # Seurat v5: join split RNA layers before DE
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

  if (!"data" %in% Layers(obj[["RNA"]]) ||
      max(obj[["RNA"]]$data[1:100, 1:10]) > 50) {
    # ensure log-normalized data exists
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  Idents(obj) <- obj@meta.data[[cluster_col]]
  message("Running presto::wilcoxauc across ",
          length(unique(Idents(obj))), " clusters ...")

  markers <- wilcoxauc(obj, seurat_assay = "RNA") %>%
    as.data.frame() %>%
    arrange(group, desc(logFC))

  if (!is.null(out_prefix)) {
    write.csv(markers, paste0(out_prefix, ".csv"), row.names = FALSE)
    top10 <- markers %>% group_by(group) %>% slice_max(logFC, n = 10)
    write.csv(top10, sub("\\.csv$", "", paste0(out_prefix, "_top10.csv")),
              row.names = FALSE)
    message("Wrote ", out_prefix, ".csv and ", out_prefix, "_top10.csv")
  }

  markers
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 02_marker_de.R <input.qs2> <out_prefix> [cluster_col]")
  }
  obj <- qs_read(args[1], nthreads = .nthreads)
  cc  <- if (length(args) >= 3) args[3] else "cluster_final"
  invisible(run_marker_de(obj, cluster_col = cc, out_prefix = args[2]))
}
