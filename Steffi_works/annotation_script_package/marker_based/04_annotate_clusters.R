# ============================================================================
# 04_annotate_clusters.R — Data-driven cluster -> cell-type labels
# ============================================================================
#
# Step 4 of the marker-based HCC annotation pipeline.
#
# Labels are derived automatically from the 15 canonical lineage module scores
# added in 03_module_scores.R (columns mod_<Set>1). Each cluster is assigned
# the lineage whose mean module score is highest (argmax), and that lineage is
# mapped to a broad cell-type label. This replaces the object-specific
# hard-coded cluster_label_map.csv so the step works for any clustering
# (e.g. the `cluster_final` column from 01b_choose_resolution.R).
#
# Outputs:
#   celltype_fine  : the winning lineage module (e.g. "T_CD8", "Hepatocyte")
#   celltype_broad : the broad lineage (e.g. "T cell", "Hepatocyte")
#   data/cluster_label_map_derived.csv : the derived cluster -> label table
#
# Usage:
#   Rscript 04_annotate_clusters.R <input.qs2> <output.qs2> [cluster_col]
#   # or source() and call annotate_clusters()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

.script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                              value = TRUE)[1]))
if (is.na(.script_dir)) .script_dir <- "."
source(file.path(.script_dir, "data", "marker_sets_15.R"))

# Map each of the 15 lineage modules to a broad cell-type label
module_to_broad <- c(
  Hepatocyte         = "Hepatocyte",
  Cholangiocyte      = "Cholangiocyte",
  T_pan              = "T cell",
  T_CD8              = "T cell",
  T_CD4              = "T cell",
  T_Treg             = "T cell",
  NK                 = "NK cell",
  B_cell             = "B/Plasma cell",
  Plasma             = "B/Plasma cell",
  Macrophage_Kupffer = "Macrophage",
  Monocyte           = "Monocyte",
  DC                 = "DC",
  Endothelial        = "Endothelial",
  Fibroblast_HSC     = "Fibroblast/HSC",
  Cycling            = "Cycling"
)

#' Assign cluster labels from the per-cluster argmax lineage module score
#'
#' @param obj         Seurat object with mod_<Set>1 columns (from step 03)
#' @param cluster_col Metadata column with cluster IDs
#' @param write_map   Optional path to write the derived label map CSV
#' @return Seurat object with celltype_fine and celltype_broad columns
annotate_clusters <- function(obj,
                              cluster_col = "cluster_final",
                              write_map = file.path(.script_dir, "data",
                                          "cluster_label_map_derived.csv")) {

  stopifnot(cluster_col %in% colnames(obj@meta.data))

  set_names <- names(marker_sets_15)
  mod_cols  <- paste0("mod_", set_names, "1")
  missing   <- setdiff(mod_cols, colnames(obj@meta.data))
  if (length(missing) > 0) {
    stop("Module-score columns missing (run 03_module_scores.R first): ",
         paste(missing, collapse = ", "))
  }

  clusters  <- as.character(obj@meta.data[[cluster_col]])
  cl_levels <- sort(unique(clusters))

  # Per-cluster mean module score, then argmax lineage per cluster
  cl_mean <- sapply(mod_cols, function(mc)
    tapply(obj@meta.data[[mc]], clusters, mean)[cl_levels])
  cl_mean <- matrix(cl_mean, nrow = length(cl_levels),
                    dimnames = list(cl_levels, set_names))

  fine_by_cluster  <- set_names[apply(cl_mean, 1, which.max)]
  broad_by_cluster <- unname(module_to_broad[fine_by_cluster])
  names(fine_by_cluster)  <- cl_levels
  names(broad_by_cluster) <- cl_levels

  idx <- match(clusters, cl_levels)
  obj$celltype_fine  <- factor(unname(fine_by_cluster[idx]))
  obj$celltype_broad <- factor(unname(broad_by_cluster[idx]))

  label_map <- data.frame(
    cluster     = cl_levels,
    fine_label  = fine_by_cluster,
    broad_label = broad_by_cluster,
    top_score   = round(apply(cl_mean, 1, max), 4),
    row.names   = NULL
  )
  message("Derived cluster labels:")
  print(label_map)
  message("\nBroad cell-type counts:")
  print(table(obj$celltype_broad))

  if (!is.null(write_map)) {
    write.csv(label_map, write_map, row.names = FALSE)
    message("Wrote derived label map: ", write_map)
  }

  obj
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 04_annotate_clusters.R <input.qs2> <output.qs2> [cluster_col]")
  }
  obj <- qs_read(args[1], nthreads = .nthreads)
  cc  <- if (length(args) >= 3) args[3] else "cluster_final"

  obj <- annotate_clusters(obj, cluster_col = cc)
  qs_save(obj, args[2], nthreads = .nthreads)

  summ <- obj@meta.data %>%
    count(!!sym(cc), celltype_fine, celltype_broad, name = "n_cells")
  write.csv(summ, "annotation_summary.csv", row.names = FALSE)
  message("Saved: ", args[2], " + annotation_summary.csv")
}
