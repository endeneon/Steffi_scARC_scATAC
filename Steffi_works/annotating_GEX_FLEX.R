#! /usr/bin/env Rscript

# Siwei 24 Mar 2026
# init ####
{
  library(BiocParallel)
  library(glmGamPoi)
  library(future)
  library(parallel)
  library(foreach)
  library(doParallel)

  library(snow)
  library(harmony)
  # library(optparse)

  library(dplyr)
  library(data.table)

  library(edgeR)

  # library(scuttle)
  library(Matrix)
  library(matrixStats)

  library(Seurat)
  # library(Signac)
  # library(loomR)
  # library(anndata)
  # library(hdf5r)
  # library(arrow)
  # library(rhdf5)

  library(stringr)

  library(qs2)
  library(fs)

  # library(SingleR)

  library(ggplot2)
  library(gplots)
  library(patchwork)
  library(scCustomize)

  # library(scuttle)

  # library(scMerge)
  # library(scater)
}

# setwd(
#   "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/
#    szhang37/projects/szhang_dev/STEREO_seq/Human_liver/R_liver"
# )
setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)
# determine if R is running in RSTUDIO/VSCode/Positron
if (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode")) {
  print("Running under RStudio/VSCode/Positron IDE, use plan(multisession)")
  session_plan <- "multisession"
} else {
  print("Running under Rscript, use plan(multicore)")
  session_plan <- "multicore"
}

# preload functions ####
get_available_workers <-
  function(x) {
    future::plan(session_plan) # check here!
    return(future::nbrOfFreeWorkers())
  }

# LSF does not expose its core allocation to future::availableCores() by
# default, so it falls back to reporting 1. Read LSB_DJOB_NUMPROC directly,
# falling back to parallelly's detection when not running under LSF.
lsf_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = NA))
available_cores <-
  if (!is.na(lsf_cores) && lsf_cores >= 1) {
    lsf_cores
  } else {
    parallelly::availableCores()
  }
print(available_cores)

# OpenMP-backed code (e.g. Rtsne) is single-process / shared-memory, so it can
# only use cores that live on ONE host. Keep all LSF slots on a single node
# (`#BSUB -R "span[hosts=1]"`), otherwise the -n 40 slots get spread across hosts
# and only the master host's share is actually usable.
# Rtsne() also IGNORES OMP_NUM_THREADS unless you pass num_threads = 0; its
# num_threads argument (default 1) otherwise wins -- which is why the log showed
# "OpenMP is working. 1 threads.". We therefore set an explicit thread count and
# forward it to each RunTSNE(num_threads = omp_threads) call below. Cap it (e.g.
# `min(8L, ...)`) if you prefer fewer threads than the full allocation.
omp_threads <-
  min(
    8L,
    max(
      1L,
      as.integer(available_cores / 2)
    )
  )
Sys.setenv(OMP_NUM_THREADS = omp_threads)

if (session_plan == "multisession") {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        12
      )
    )
} else {
  # ~50 GB/worker peak on the FindIntegrationAnchors step; 16 workers ~ 800 GB,
  # leaving margin under the 1000 GB (20 x 50 GB) LSF reservation.
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        16
      )
    )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  set.seed(42)
  # Seurat's parallelized steps (e.g. IntegrateLayers / FindIntegrationAnchors)
  # call future_lapply() without future.seed = TRUE, so future warns about
  # "unreliable" RNG. We cannot pass future.seed through Seurat, and our
  # stochastic steps are already explicitly seeded (set.seed(42) above, plus
  # seed.use = 42 / random.seed = 42 on each Seurat call), so silence the check.
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = workers_2_use * 20 * 1024^3) # 20 G per thread
  future::plan(
    session_plan, # Do NOT use "multisession" here if use LSF, use "multicore" instead
    workers = workers_2_use
  )
}

# functions ####
## extract raw count matrix from Seurat object
extract_raw_count_mtx <-
  function(
    object,
    assay = "RNA"
  ) {
    object <-
      UpdateSeuratObject(object)
    DefaultAssay(object) <- assay
    dimnames_count_mtx <-
      list(
        gene_id = rownames(object),
        barcodes = colnames(object)
      )
    if (
      is.null(dimnames_count_mtx$gene_id) |
        is.null(dimnames_count_mtx$barcodes)
    ) {
      return(NA)
    } else {
      raw_count_mtx <-
        object[[assay]]@layers$counts
      rownames(raw_count_mtx) <-
        dimnames_count_mtx$gene_id
      colnames(raw_count_mtx) <-
        dimnames_count_mtx$barcodes
      return(raw_count_mtx)
    }
  }

# find PCs
# use Surbhi's function with changes
# this function will return npcs only, no figures
findPcsElbow <-
  function(obj) {
    # credit: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
    # calculate where the principal components start to elbow by taking the larger value of:
    #  1. The point where the principal components only contribute 5% of standard deviation and
    # 		the principal components cumulatively contribute 90% of the standard deviation.
    #  2. The point where the percent change in variation between the consecutive PCs is less than 0.1%
    # 		Determine percent of variation associated with each PC
    stopifnot(is(obj, "Seurat"))

    # helper: percent and cumulative percent of stdev per PC
    computePct <- function(obj) {
      pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
      list(pct = pct, cumu = cumsum(pct))
    }

    # (re)run PCA if there is no reduction yet, or if the cumu > 90 & pct < 5
    # condition is never satisfied simultaneously with the current PCs
    needsPca <- !"pca" %in% Reductions(obj)
    if (!needsPca) {
      p <- computePct(obj)
      needsPca <- length(which(p$cumu > 90 & p$pct < 5)) == 0
    }
    if (needsPca) {
      obj <-
        obj %>%
        NormalizeData(verbose = T) %>%
        FindVariableFeatures(verbose = T) %>%
        ScaleData(verbose = T) %>%
        RunPCA(npcs = 50, verbose = T)
    }

    p <- computePct(obj)
    pct <- p$pct
    cumu <- p$cumu
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <-
      which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <-
      sort(
        which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
        decreasing = T
      )[1] +
      1
    # Minimum of the two calculation; i.e. PCs covering the majority of the variation in the data.
    pcs <-
      min(c(co1, co2), na.rm = TRUE)
    return(ndims = pcs)
  }

merged_liver_obj <-
  qs_read("merged_all_samples_integrated_seurat_obj.qs2", nthreads = 8)

ref_obj <-
  qs_read("seurat_gex_with_predicted_annotation.qs2", nthreads = 8)

############
## Transfer cell-identity labels from ref_obj (predicted.annotation)
## onto merged_liver_obj via Seurat anchor-based label transfer.
############

## Reference: use RNA assay, log-normalize, HVGs, scale, PCA.
## Label transfer needs a normalized `data` layer and a `pca` reduction,
## neither of which ref_obj carries out of the box (only `counts` + IterativeLSI).
DefaultAssay(ref_obj) <- "RNA"
ref_obj <-
  ref_obj %>%
  NormalizeData(verbose = T) %>%
  FindVariableFeatures(verbose = T) %>%
  ScaleData(verbose = T) %>%
  RunPCA(npcs = 50, verbose = T)

## Choose transfer dims from the reference PCA elbow
n_dims <-
  findPcsElbow(ref_obj)

## Query: use RNA assay, join the per-sample split layers, then log-normalize.
## merged_liver_obj defaults to SCT with 36 un-joined counts.*/data.* RNA layers,
## so anchoring on RNA requires joining and normalizing first.
DefaultAssay(merged_liver_obj) <- "RNA"
merged_liver_obj[["RNA"]] <-
  JoinLayers(merged_liver_obj[["RNA"]])
merged_liver_obj <-
  NormalizeData(merged_liver_obj, verbose = T)

## Batch-aware label transfer: anchor and transfer separately within each
## preparation batch (GEX / FLEX / multiome), then recombine. This avoids
## anchoring across chemistries that the integration was built to correct.
query_batches <-
  SplitObject(merged_liver_obj, split.by = "preparation")
gc()

predictions_list <-
  lapply(
    query_batches,
    function(query_sub) {
      anchors_sub <-
        FindTransferAnchors(
          reference = ref_obj,
          query = query_sub,
          dims = 1:n_dims,
          reference.reduction = "pca",
          verbose = T
        )
      pred_sub <-
        TransferData(
          anchorset = anchors_sub,
          refdata = ref_obj$predicted.annotation,
          dims = 1:n_dims,
          verbose = T
        )
      # free the batch's anchors and subset before the next iteration
      rm(anchors_sub, query_sub)
      gc()
      pred_sub
    }
  )

## Stack per-batch predictions and restore original cell order
predictions <-
  do.call(rbind, predictions_list)
predictions <-
  predictions[colnames(merged_liver_obj), ]

## Add the transferred labels (and prediction scores) to query metadata
merged_liver_obj <-
  AddMetaData(merged_liver_obj, metadata = predictions)

merged_liver_obj$annotation <-
  merged_liver_obj$predicted.id

## Per-batch summary of transfer confidence (prediction.score.max)
prediction_score_summary <-
  merged_liver_obj@meta.data %>%
  dplyr::group_by(preparation) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    min_score = min(prediction.score.max),
    median_score = median(prediction.score.max),
    mean_score = mean(prediction.score.max),
    max_score = max(prediction.score.max),
    # fraction of cells transferred with low confidence
    frac_below_0.5 = mean(prediction.score.max < 0.5),
    .groups = "drop"
  )
print(prediction_score_summary)

qs_save(
  merged_liver_obj,
  "merged_all_samples_integrated_seurat_obj_annotated.qs2",
  nthreads = 8
)
