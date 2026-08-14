# ============================================================================
# 01b_choose_resolution.R — Pick a clustering resolution giving 10-12 clusters
# ============================================================================
#
# Step 1b of the marker-based HCC annotation pipeline.
#
# The all-samples integrated object is over-clustered at SCT_snn_res.1 (66
# clusters). This step searches the resolution that yields a target number of
# clusters (default 10-12) and stores it in a stable column `cluster_final`
# used by every downstream step, so the pipeline no longer depends on a
# hard-coded resolution/column.
#
# Search strategy (reuses the existing SNN graph, so it is fast):
#   1. Start at start_res (default 1.0).
#   2. Halve the resolution until the cluster count drops to <= target_max.
#   3. If the count is now < target_min, bisect between the last two
#      resolutions (mean of the bracketing low/high resolutions) until the
#      count lands in [target_min, target_max].
#
# Usage:
#   Rscript 01b_choose_resolution.R <input.qs2> <output.qs2> \
#       [graph_name] [target_min] [target_max] [start_res]
#   # or source() and call choose_resolution()

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
})

.nthreads <- as.integer(Sys.getenv("NTHREADS", "8"))

#' Search for a resolution that yields target_min..target_max clusters
#'
#' @param obj        Seurat object with a precomputed SNN graph
#' @param graph_name Name of the SNN graph to cluster on (default "SCT_snn")
#' @param target_min,target_max Inclusive target cluster-count band
#' @param start_res  Starting resolution (default 1.0)
#' @param max_iter   Safety cap on search iterations
#' @return list(obj = object with `cluster_final`, resolution, k, history)
choose_resolution <- function(obj,
                              graph_name = "SCT_snn",
                              target_min = 10,
                              target_max = 12,
                              start_res  = 1.0,
                              max_iter   = 40) {

  stopifnot(graph_name %in% Graphs(obj))

  history <- data.frame(res = numeric(0), k = integer(0))
  n_clusters_at <- function(res) {
    obj <<- FindClusters(obj, graph.name = graph_name,
                         resolution = res, verbose = FALSE)
    k <- length(unique(Idents(obj)))
    history[nrow(history) + 1L, ] <<- list(res, k)
    message(sprintf("  res=%.6f -> %d clusters", res, k))
    k
  }

  res <- start_res
  k   <- n_clusters_at(res)

  # Phase 1: halve until at/under the upper bound; remember the last
  # "too many" resolution as the high bracket.
  hi_res <- res
  while (k > target_max && nrow(history) < max_iter) {
    hi_res <- res
    res    <- res / 2
    k      <- n_clusters_at(res)
  }

  # Phase 2: if too few, bisect between lo_res (too few) and hi_res (too many)
  lo_res <- res
  while (!(k >= target_min && k <= target_max) && nrow(history) < max_iter) {
    if (k < target_min) {
      lo_res <- res
    } else {                       # k > target_max
      hi_res <- res
    }
    res <- mean(c(lo_res, hi_res))
    k   <- n_clusters_at(res)
  }

  # Choose the best resolution: prefer one inside the band, else the closest
  # to the band midpoint, then re-cluster once so Idents match cluster_final.
  in_band <- history$k >= target_min & history$k <= target_max
  if (any(in_band)) {
    best_res <- history$res[which(in_band)[1]]
  } else {
    mid      <- (target_min + target_max) / 2
    best_res <- history$res[which.min(abs(history$k - mid))]
    warning(sprintf("No resolution gave %d-%d clusters after %d iterations; ",
                    target_min, target_max, nrow(history)),
            sprintf("using res=%.6f (closest: %d clusters).",
                    best_res, history$k[which.min(abs(history$k - mid))]))
  }

  obj <- FindClusters(obj, graph.name = graph_name,
                      resolution = best_res, verbose = FALSE)
  obj$cluster_final <- factor(as.integer(as.character(Idents(obj))))
  k_final <- length(levels(obj$cluster_final))

  message(sprintf("Selected resolution %.6f -> %d clusters (column 'cluster_final')",
                  best_res, k_final))
  list(obj = obj, resolution = best_res, k = k_final, history = history)
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript 01b_choose_resolution.R <input.qs2> <output.qs2> [graph_name] [target_min] [target_max] [start_res]")
  }
  graph_name <- if (length(args) >= 3) args[3] else "SCT_snn"
  target_min <- if (length(args) >= 4) as.integer(args[4]) else 10L
  target_max <- if (length(args) >= 5) as.integer(args[5]) else 12L
  start_res  <- if (length(args) >= 6) as.numeric(args[6]) else 1.0

  obj <- qs_read(args[1], nthreads = .nthreads)
  res <- choose_resolution(obj, graph_name = graph_name,
                           target_min = target_min, target_max = target_max,
                           start_res = start_res)
  qs_save(res$obj, args[2], nthreads = .nthreads)
  write.csv(res$history, "resolution_search_history.csv", row.names = FALSE)
  message("Saved: ", args[2], " + resolution_search_history.csv")
}
