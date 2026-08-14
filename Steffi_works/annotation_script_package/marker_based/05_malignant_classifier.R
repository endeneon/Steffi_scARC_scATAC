# ============================================================================
# 05_malignant_classifier.R — Malignant vs normal hepatocyte classifier
# ============================================================================
#
# Step 5 of the marker-based HCC annotation pipeline.
#
# Original analysis context (see annotation_dictionary.md, "Limitations" #1):
#   - Applied to the 59,105 cells with celltype_broad == "Hepatocyte"
#   - Two AddModuleScore panels on the RNA assay:
#       malignant (9 genes): AFP, GPC3, MDK, TOP2A, IGF2BP1, FOXM1, PEG10,
#                            BIRC5, SPP1
#       mature    (8 genes): ALB, HNF4A, CYP3A4, CYP2E1, TF, APOA1, CPS1,
#                            SERPINA1
#   - NOTE: the initial plan used delta scoring (malignant - mature), but on
#     this data delta collapsed most AFP+/GPC3+ cells into "ambiguous" because
#     HCC tumor cells retain moderate mature-hepatocyte expression
#     (co-expression, not replacement). The delivered classifier therefore
#     uses QUARTILES OF THE MALIGNANT SCORE ALONE on hepatocytes:
#       mod_hcc_malignant1 > Q3  -> malignant_hepatocyte (n=14,776; 25%)
#       mod_hcc_malignant1 < Q1  -> normal_hepatocyte    (n=14,776; 25%)
#       middle 50%               -> ambiguous            (n=29,553)
#   - Non-hepatocytes are labeled "non_hepatocyte"
#   - Both module scores and the delta (hep_delta) are retained in metadata
#   - NOT CNV-validated: an inferCNV/CopyKAT run using clusters
#     2,4,5,6,8,9,11,12,13 as diploid reference is recommended to confirm
#
# Usage:
#   Rscript 05_malignant_classifier.R <input.qs2> <output.qs2>
#   # or source() and call classify_malignant()
suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

malignant_genes <- c("AFP", "GPC3", "MDK", "TOP2A", "IGF2BP1",
                     "FOXM1", "PEG10", "BIRC5", "SPP1")
mature_genes    <- c("ALB", "HNF4A", "CYP3A4", "CYP2E1", "TF",
                     "APOA1", "CPS1", "SERPINA1")

#' Quartile-based malignant classifier for hepatocytes
#'
#' @param obj            Seurat object with celltype_broad column
#' @param hepatocyte_label Value of celltype_broad marking hepatocytes
#' @return Seurat object with mod_hcc_malignant1, mod_mature_hepatocyte1,
#'         hep_delta, and malignant_status columns
classify_malignant <- function(obj, hepatocyte_label = "Hepatocyte") {

  stopifnot("celltype_broad" %in% colnames(obj@meta.data))
  DefaultAssay(obj) <- "RNA"
  # Seurat v5: join split RNA layers and log-normalize before scoring
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- NormalizeData(obj, verbose = FALSE)

  missing <- setdiff(c(malignant_genes, mature_genes), rownames(obj))
  if (length(missing) > 0) {
    stop("Panel genes missing from object: ", paste(missing, collapse = ", "))
  }

  obj <- AddModuleScore(obj, features = list(malignant_genes),
                        name = "mod_hcc_malignant", assay = "RNA")
  obj <- AddModuleScore(obj, features = list(mature_genes),
                        name = "mod_mature_hepatocyte", assay = "RNA")
  obj$hep_delta <- obj$mod_hcc_malignant1 - obj$mod_mature_hepatocyte1

  is_hep <- obj$celltype_broad == hepatocyte_label
  hep_scores <- obj$mod_hcc_malignant1[is_hep]
  q1 <- quantile(hep_scores, 0.25, na.rm = TRUE)
  q3 <- quantile(hep_scores, 0.75, na.rm = TRUE)
  message(sprintf("Hepatocyte malignant-score quartiles: Q1=%.4f Q3=%.4f",
                  q1, q3))

  obj$malignant_status <- "non_hepatocyte"
  obj$malignant_status[is_hep & obj$mod_hcc_malignant1 > q3] <- "malignant_hepatocyte"
  obj$malignant_status[is_hep & obj$mod_hcc_malignant1 < q1] <- "normal_hepatocyte"
  obj$malignant_status[is_hep & obj$mod_hcc_malignant1 >= q1 &
                       obj$mod_hcc_malignant1 <= q3] <- "ambiguous"
  obj$malignant_status <- factor(obj$malignant_status,
                                 levels = c("normal_hepatocyte",
                                            "malignant_hepatocyte",
                                            "ambiguous",
                                            "non_hepatocyte"))

  message("Malignant status counts:")
  print(table(obj$malignant_status))

  obj
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 05_malignant_classifier.R <input.qs2> <output.qs2> [cluster_col]")
  }
  obj <- qs_read(args[1], nthreads = .nthreads)
  cc  <- if (length(args) >= 3) args[3] else "cluster_final"
  obj <- classify_malignant(obj)
  qs_save(obj, args[2], nthreads = .nthreads)

  # Per-cluster malignant-status breakdown
  cl <- obj@meta.data %>%
    count(.data[[cc]], malignant_status) %>%
    group_by(.data[[cc]]) %>%
    mutate(pct = 100 * n / sum(n))
  write.csv(cl, "malignant_status_per_cluster.csv", row.names = FALSE)
  message("Saved: ", args[2], " + malignant_status_per_cluster.csv")
}
