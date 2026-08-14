# ============================================================================
# 03_compare_methods.R — Marker-based vs reference-based label agreement
# ============================================================================
#
# Step 3 of the reference-based annotation pipeline.
#
# Compares the marker-based labels (celltype_broad, from Package A) with the
# SingleR labels (singler_label_cell, from 02_singler_annotate.R) on the same
# object. The two taxonomies differ; a harmonization map is applied:
#
#   marker-based broad label   |  Lu et al. 2022 (GSE149614) celltype
#   ---------------------------+---------------------------------------
#   Hepatocyte                 |  Hepatocyte
#   T cell, T/NK cell          |  T/NK
#   Macrophage, Monocyte       |  Myeloid
#   B/Plasma cell              |  B
#   Endothelial                |  Endothelial
#   Fibroblast/HSC             |  Fibroblast
#   Cycling                    |  (no equivalent — lineage-agnostic state)
#
# Outputs: confusion matrix CSV + per-cluster agreement CSV.
#
# Usage:
#   Rscript 03_compare_methods.R <annotated.qs2> <out_prefix>
#   # or source() and call compare_methods()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
})

# Threads granted by the scheduler (exported as NTHREADS in the bsub script).
nthreads <- as.integer(Sys.getenv("NTHREADS", unset = "8"))
if (is.na(nthreads) || nthreads < 1L) nthreads <- 8L

# Harmonization map: marker-based broad label -> Lu et al. celltype
broad_to_lu <- c(
  "Hepatocyte"     = "Hepatocyte",
  "T cell"         = "T/NK",
  "T/NK cell"      = "T/NK",
  "Macrophage"     = "Myeloid",
  "Monocyte"       = "Myeloid",
  "B/Plasma cell"  = "B",
  "Endothelial"    = "Endothelial",
  "Fibroblast/HSC" = "Fibroblast",
  "Cycling"        = NA_character_   # no Lu-equivalent lineage
)

#' Compare marker-based and SingleR labels
#'
#' @param obj Seurat object with celltype_broad and singler_label_cell
#' @param cluster_col Cluster column for per-cluster agreement
#' @param out_prefix File prefix for CSV outputs (NULL = skip writing)
#' @return list(confusion = table, per_cluster = data.frame, overall_agreement)
compare_methods <- function(obj,
                            cluster_col = "SCT_snn_res.0.5",
                            out_prefix = "marker_vs_singler") {

  stopifnot(all(c("celltype_broad", "singler_label_cell") %in%
                colnames(obj@meta.data)))

  md <- obj@meta.data
  md$marker_harmonized <- unname(broad_to_lu[as.character(md$celltype_broad)])

  # Confusion matrix (harmonized marker label x SingleR label)
  confusion <- table(marker_based = md$marker_harmonized,
                     singler      = md$singler_label_cell,
                     useNA = "ifany")

  comparable <- !is.na(md$marker_harmonized)
  agreement <- mean(md$marker_harmonized[comparable] ==
                    md$singler_label_cell[comparable])
  message(sprintf("Overall cell-level agreement (Cycling excluded): %.1f%%",
                  100 * agreement))

  # Per-cluster agreement
  per_cluster <- md %>%
    filter(!is.na(marker_harmonized)) %>%
    mutate(match = marker_harmonized == singler_label_cell) %>%
    group_by(.data[[cluster_col]], celltype_fine) %>%
    summarise(n = n(),
              pct_agreement = 100 * mean(match, na.rm = TRUE),
              top_singler = names(which.max(table(singler_label_cell))),
              .groups = "drop")

  if (!is.null(out_prefix)) {
    write.csv(as.data.frame.matrix(confusion),
              paste0(out_prefix, "_confusion.csv"))
    write.csv(per_cluster, paste0(out_prefix, "_per_cluster.csv"),
              row.names = FALSE)
    message("Wrote ", out_prefix, "_confusion.csv and _per_cluster.csv")
  }

  list(confusion = confusion, per_cluster = per_cluster,
       overall_agreement = agreement)
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 03_compare_methods.R <annotated.qs2> <out_prefix>")
  }
  obj <- qs_read(args[1], nthreads = nthreads)
  invisible(compare_methods(obj, out_prefix = args[2]))
}
