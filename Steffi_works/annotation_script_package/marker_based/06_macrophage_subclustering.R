# ============================================================================
# 06_macrophage_subclustering.R — Macrophage re-integration + subtyping
# ============================================================================
#
# Step 6 of the marker-based HCC annotation pipeline.
#
# Original analysis context (see annotation_dictionary.md,
# "Macrophage subclustering"):
#   - Input: the 25,516 cells with celltype_broad == "Macrophage"
#     (whole-tissue clusters 4 + 9)
#   - Re-integrated with Harmony (per-sample batch = orig.ident)
#   - Clustered at resolution 0.3 (column mac_res.0.3)
#   - Reduced-dim slots on the macrophage-only object: harmony.mac, umap.mac
#   - Subtypes are assigned by MARKER MODULE SCORE, not by cluster ID: each
#     subcluster is labeled with the subtype whose marker panel scores highest
#     (mac_subtype_markers below). This is robust to the unstable Louvain
#     cluster numbering and to runs that yield a different number of
#     subclusters. Subtypes:
#       LYVE1_TRM                      LYVE1/F13A1/MRC1/MAF resident
#       TFRC_TAM                       TFRC/MSR1/ACP5 iron-recycling TAM
#       Doublet_hepatocyte             AFP/GPC3/HNF4A/CYP2C19
#       Doublet_hepatocyte_acute_phase SAA1/SAA2/HP/FGA/CRP
#       Doublet_Tcell                  THEMIS/CD3D/CD8A
#       Doublet_endothelial            LDB2/VEGFC/FLT1/PECAM1
#   - Doublets were FLAGGED, not removed (marker-based classification:
#     co-expression of non-macrophage lineage markers with macrophage counts;
#     no formal doublet detector was run)
#   - For downstream macrophage analysis, filter to
#     macrophage_subtype %in% c("LYVE1_TRM", "TFRC_TAM")
#
# Usage:
#   Rscript 06_macrophage_subclustering.R <input.qs2> <output.qs2> [sample_col]
#   # or source() and call subcluster_macrophages()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(harmony)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

# Subtype marker panels. Each subcluster is assigned to the subtype whose
# marker module score is highest, so labeling does not depend on the
# (run-unstable) Louvain cluster numbering or a fixed number of subclusters.
mac_subtype_markers <- list(
  LYVE1_TRM                      = c("LYVE1", "F13A1", "MRC1", "MAF"),
  TFRC_TAM                       = c("TFRC", "MSR1", "ACP5"),
  Doublet_hepatocyte             = c("AFP", "GPC3", "HNF4A", "CYP2C19"),
  Doublet_hepatocyte_acute_phase = c("SAA1", "SAA2", "HP", "FGA", "CRP"),
  Doublet_Tcell                  = c("THEMIS", "CD3D", "CD8A"),
  Doublet_endothelial            = c("LDB2", "VEGFC", "FLT1", "PECAM1")
)

#' Subset, re-integrate and subtype the macrophage compartment
#'
#' @param obj        Seurat object with celltype_broad column
#' @param sample_col Metadata column with sample IDs (default orig.ident)
#' @param resolution Clustering resolution (default 0.3, as original)
#' @return list(mac = macrophage-only Seurat object with macrophage_subtype,
#'              subtype_table = data.frame of subtype counts)
subcluster_macrophages <- function(obj,
                                   sample_col = "orig.ident",
                                   resolution = 0.3) {

  stopifnot("celltype_broad" %in% colnames(obj@meta.data))

  if (sum(obj$celltype_broad == "Macrophage", na.rm = TRUE) == 0) {
    warning("No cells labeled 'Macrophage'; skipping macrophage subclustering. ",
            "The clustering resolution may be too coarse to resolve myeloid cells.")
    return(list(mac = NULL, subtype_table = NULL, cluster_to_subtype = NULL))
  }

  mac <- subset(obj, subset = celltype_broad == "Macrophage")
  message(sprintf("Macrophage subset: %d cells", ncol(mac)))

  DefaultAssay(mac) <- "SCT"
  mac <- mac %>%
    RunPCA(assay = "SCT", verbose = FALSE) %>%
    RunHarmony(sample_col, plot_convergence = FALSE,
               reduction.save = "harmony.mac") %>%
    RunUMAP(reduction = "harmony.mac", dims = 1:30,
            reduction.name = "umap.mac") %>%
    FindNeighbors(reduction = "harmony.mac", dims = 1:30) %>%
    FindClusters(resolution = resolution,
                 cluster.name = paste0("mac_res.", resolution))

  cl_col <- paste0("mac_res.", resolution)
  clusters <- as.character(mac@meta.data[[cl_col]])

  # Marker-based labeling: score each subtype panel per cell, then assign
  # every subcluster to the subtype with the highest mean module score.
  DefaultAssay(mac) <- "RNA"
  mac[["RNA"]] <- JoinLayers(mac[["RNA"]])
  mac <- NormalizeData(mac, verbose = FALSE)

  for (sub in names(mac_subtype_markers)) {
    genes <- intersect(mac_subtype_markers[[sub]], rownames(mac))
    if (length(genes) == 0) {
      warning("No marker genes present for subtype '", sub, "'; skipped.")
      next
    }
    mac <- AddModuleScore(mac, features = list(genes),
                          name = paste0("subscore_", sub), assay = "RNA")
  }

  score_cols <- paste0("subscore_", names(mac_subtype_markers), "1")
  score_cols <- score_cols[score_cols %in% colnames(mac@meta.data)]
  cl_mean <- sapply(score_cols, function(sc)
    tapply(mac@meta.data[[sc]], clusters, mean))
  cl_mean <- matrix(cl_mean, nrow = length(unique(clusters)),
                    dimnames = list(sort(unique(clusters)), score_cols))
  subtype_names <- sub("1$", "", sub("^subscore_", "", score_cols))
  cluster_to_subtype <- setNames(subtype_names[apply(cl_mean, 1, which.max)],
                                 rownames(cl_mean))

  message("Subcluster -> subtype (argmax marker module score):")
  print(cluster_to_subtype)

  mac$macrophage_subtype <- factor(unname(cluster_to_subtype[clusters]),
                                   levels = names(mac_subtype_markers))

  subtype_table <- as.data.frame(table(mac$macrophage_subtype))
  colnames(subtype_table) <- c("subtype", "n_cells")
  subtype_table$pct <- 100 * subtype_table$n_cells / ncol(mac)
  print(subtype_table)

  list(mac = mac, subtype_table = subtype_table,
       cluster_to_subtype = cluster_to_subtype)
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 06_macrophage_subclustering.R <input.qs2> <output.qs2> [sample_col]")
  }
  obj <- qs_read(args[1], nthreads = .nthreads)
  sc  <- if (length(args) >= 3) args[3] else "orig.ident"

  res <- subcluster_macrophages(obj, sample_col = sc)

  # Write subtype back onto the full object (NA for non-macrophages)
  obj$macrophage_subtype <- NA_character_
  if (!is.null(res$mac)) {
    obj$macrophage_subtype[colnames(res$mac)] <-
      as.character(res$mac$macrophage_subtype)
  }
  obj$macrophage_subtype <- factor(obj$macrophage_subtype,
                                   levels = names(mac_subtype_markers))

  qs_save(obj, args[2], nthreads = .nthreads)
  if (!is.null(res$subtype_table)) {
    write.csv(res$subtype_table, "macrophage_subtype_summary.csv",
              row.names = FALSE)
  }
  message("Saved: ", args[2], " + macrophage_subtype_summary.csv")
}
