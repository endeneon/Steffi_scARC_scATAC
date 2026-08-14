# ============================================================================
# 03_module_scores.R — AddModuleScore on the 15 canonical HCC-TME marker sets
# ============================================================================
#
# Step 3 of the marker-based HCC annotation pipeline.
#
# Original analysis context:
#   - 15 lineage module sets (see data/marker_sets_15.R for provenance)
#   - Scored on the RNA assay (log-normalized), all cells
#   - Output columns: mod_<SetName>1 in meta.data
#   - Per-cluster means exported to canonical_module_scores_per_cluster.csv
#   - Canonical per-gene cluster means exported to
#     canonical_marker_avg_expression.csv (69 genes x 16 clusters)
#
# Usage:
#   Rscript 03_module_scores.R <input.qs2> <output.qs2> [cluster_col]
#   # or source() and call add_module_scores()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

# Marker sets live next to this script
.script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                              value = TRUE)[1]))
if (is.na(.script_dir)) .script_dir <- "."
source(file.path(.script_dir, "data", "marker_sets_15.R"))

#' Score the 15 canonical lineage modules on all cells
#'
#' @param obj         Seurat object (post-QC, RNA assay log-normalized)
#' @param cluster_col Cluster column for per-cluster summary (default
#'                    SCT_snn_res.0.5)
#' @return list(obj = object with mod_*1 columns,
#'              per_cluster = matrix of cluster-mean module scores,
#'              gene_avg = matrix of cluster-mean expression for the 69 genes)
add_module_scores <- function(obj, cluster_col = "cluster_final") {

  DefaultAssay(obj) <- "RNA"
  # Seurat v5: join split RNA layers and log-normalize before scoring
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- NormalizeData(obj, verbose = FALSE)

  # Warn about genes missing from this dataset (none were missing originally)
  all_genes <- unique(unlist(marker_sets_15))
  missing <- setdiff(all_genes, rownames(obj))
  if (length(missing) > 0) {
    warning("Genes absent from object and dropped from sets: ",
            paste(missing, collapse = ", "))
    marker_sets_15 <<- lapply(marker_sets_15, intersect, rownames(obj))
  }

  for (set_name in names(marker_sets_15)) {
    obj <- AddModuleScore(obj,
                          features = list(marker_sets_15[[set_name]]),
                          name     = paste0("mod_", set_name),
                          assay    = "RNA")
  }
  # AddModuleScore appends "1" -> columns are mod_<SetName>1

  mod_cols <- paste0("mod_", names(marker_sets_15), "1")
  clusters <- obj@meta.data[[cluster_col]]

  per_cluster <- sapply(mod_cols, function(mc)
    tapply(obj@meta.data[[mc]], clusters, mean))
  rownames(per_cluster) <- sort(unique(clusters))

  # Per-gene cluster means (log-normalized data) for the full 69-gene panel
  dat <- GetAssayData(obj, assay = "RNA", layer = "data")
  genes <- intersect(all_genes, rownames(dat))
  gene_avg <- sapply(sort(unique(clusters)), function(cl) {
    cells <- which(clusters == cl)
    Matrix::rowMeans(dat[genes, cells, drop = FALSE])
  })
  colnames(gene_avg) <- sort(unique(clusters))

  list(obj = obj, per_cluster = per_cluster, gene_avg = gene_avg)
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 03_module_scores.R <input.qs2> <output.qs2> [cluster_col]")
  }
  obj <- qs_read(args[1], nthreads = .nthreads)
  cc  <- if (length(args) >= 3) args[3] else "cluster_final"

  res <- add_module_scores(obj, cluster_col = cc)
  qs_save(res$obj, args[2], nthreads = .nthreads)
  write.csv(res$per_cluster, "canonical_module_scores_per_cluster.csv")
  write.csv(res$gene_avg,    "canonical_marker_avg_expression.csv")
  message("Saved object + module-score and gene-expression summaries")
}
