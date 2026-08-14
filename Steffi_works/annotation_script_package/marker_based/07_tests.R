# ============================================================================
# 07_tests.R — Sanity test suite for the marker-based annotation
# ============================================================================
#
# Step 7 of the marker-based HCC annotation pipeline.
#
# Reproduces object-agnostic sanity checks for any clustering/annotation:
#   1. (optional) cell count after QC
#   2. No NA in celltype_fine / celltype_broad; every cluster labeled
#   3. Label consistency: each cluster's assigned lineage equals the argmax of
#      its 15 lineage module scores
#   4. Malignant classifier: ~25/25/50 split among hepatocytes
#   5. Macrophage subtypes: LYVE1 peaks in LYVE1_TRM, TFRC in TFRC_TAM,
#      CD68 depressed in hepatocyte-doublet clusters
#   6. Contamination markers peak in their expected doublet subtype
#
# Usage:
#   Rscript 07_tests.R <annotated.qs2> [cluster_col] [expected_cells]
#   # or source() and call run_tests(obj)
#
# Exit code 0 if all applicable tests pass, 1 otherwise.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

#' Run the sanity test suite
#' @param obj Annotated Seurat object (output of steps 01-06)
#' @param cluster_col Cluster column used for annotation
#' @param expected_cells Expected post-QC cell count (NULL = skip test 1)
#' @return invisible TRUE if all tests pass
run_tests <- function(obj, cluster_col = "cluster_final",
                      expected_cells = NULL) {

  results <- list()
  check <- function(name, expr) {
    ok <- isTRUE(expr)
    results[[name]] <<- ok
    message(sprintf("[%s] %s", if (ok) "PASS" else "FAIL", name))
    ok
  }

  # 1. Cell count
  if (!is.null(expected_cells)) {
    check("cell count matches expected", ncol(obj) == expected_cells)
  }

  # 2. Label coverage
  check("no NA in celltype_fine",  !any(is.na(obj$celltype_fine)))
  check("no NA in celltype_broad", !any(is.na(obj$celltype_broad)))
  check("every cluster labeled",
        cluster_col %in% colnames(obj@meta.data) &&
        length(unique(obj$celltype_fine)) >= 1 &&
        !any(is.na(obj$celltype_fine)))

  # 3. Label consistency: assigned lineage == argmax lineage module per cluster
  mod_cols <- grep("^mod_.*1$", colnames(obj@meta.data), value = TRUE)
  lineage_mods <- setdiff(mod_cols,
                          c("mod_hcc_malignant1", "mod_mature_hepatocyte1"))
  if (length(lineage_mods) >= 15 && cluster_col %in% colnames(obj@meta.data)) {
    cl        <- as.character(obj@meta.data[[cluster_col]])
    cl_levels <- sort(unique(cl))
    set_names <- sub("1$", "", sub("^mod_", "", lineage_mods))
    cl_mean   <- sapply(lineage_mods, function(mc)
      tapply(obj@meta.data[[mc]], cl, mean)[cl_levels])
    cl_mean   <- matrix(cl_mean, nrow = length(cl_levels),
                        dimnames = list(cl_levels, set_names))
    argmax_fine <- set_names[apply(cl_mean, 1, which.max)]
    assigned    <- tapply(as.character(obj$celltype_fine), cl,
                          function(x) x[1])[cl_levels]
    check("cluster labels match argmax lineage module",
          all(argmax_fine == assigned))
  }

  # 4. Malignant classifier split
  if ("malignant_status" %in% colnames(obj@meta.data)) {
    hep <- obj$malignant_status[obj$celltype_broad == "Hepatocyte"]
    frac <- prop.table(table(hep))
    check("malignant ~25% of hepatocytes",
          abs(frac["malignant_hepatocyte"] - 0.25) < 0.02)
    check("normal ~25% of hepatocytes",
          abs(frac["normal_hepatocyte"] - 0.25) < 0.02)
  }

  # 5-6. Macrophage subtype sanity (only if subclustering was run)
  if ("macrophage_subtype" %in% colnames(obj@meta.data) &&
      !all(is.na(obj$macrophage_subtype))) {
    mac <- subset(obj, subset = !is.na(macrophage_subtype))
    DefaultAssay(mac) <- "RNA"
    # Seurat v5: join split RNA layers before reading the data layer
    mac[["RNA"]] <- JoinLayers(mac[["RNA"]])
    dat <- GetAssayData(mac, assay = "RNA", layer = "data")
    st  <- mac$macrophage_subtype
    m   <- function(g) tapply(dat[g, ], st, mean)

    check("LYVE1 peaks in LYVE1_TRM", names(which.max(m("LYVE1"))) == "LYVE1_TRM")
    check("TFRC peaks in TFRC_TAM",   names(which.max(m("TFRC"))) == "TFRC_TAM")
    # AFP marks hepatocyte contamination in either doublet subtype; GPC3 is the
    # cleaner discriminator for the plain hepatocyte doublet
    check("AFP peaks in a hepatocyte doublet",
          names(which.max(m("AFP"))) %in%
            c("Doublet_hepatocyte", "Doublet_hepatocyte_acute_phase"))
    check("CD3D peaks in Doublet_Tcell",
          names(which.max(m("CD3D"))) == "Doublet_Tcell")
    check("VWF peaks in Doublet_endothelial",
          names(which.max(m("VWF"))) == "Doublet_endothelial")
    cd68 <- m("CD68")
    check("CD68 depressed in hepatocyte doublets",
          cd68["Doublet_hepatocyte"] < cd68["LYVE1_TRM"] &&
          cd68["Doublet_hepatocyte_acute_phase"] < cd68["LYVE1_TRM"])
  }

  n_pass <- sum(unlist(results))
  n_tot  <- length(results)
  message(sprintf("\n=== %d/%d tests passed ===", n_pass, n_tot))
  invisible(n_pass == n_tot)
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript 07_tests.R <annotated.qs2> [cluster_col] [expected_cells]")
  obj <- qs_read(args[1], nthreads = .nthreads)
  cc  <- if (length(args) >= 2) args[2] else "cluster_final"
  ec  <- if (length(args) >= 3) as.integer(args[3]) else NULL
  ok  <- run_tests(obj, cluster_col = cc, expected_cells = ec)
  quit(status = if (ok) 0 else 1)
}
