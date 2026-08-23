#! /usr/bin/env Rscript

# init
{
  library(Seurat)
  library(gplots)
  library(ArchR)
  library(future)
  library(stringr)
  # library(pheatmap)

  library(BiocParallel)
  # library(BiocParallel.FutureParam)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(doFuture)
  library(snow)

  library(Matrix)
  library(matrixStats)

  library(qs2)
  library(fs)

  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg38)

  if (
    interactive() &&
      ((Sys.getenv("TERM_PROGRAM") == "vscode") &&
        (Sys.getenv("POSITRON") != "1"))
  ) {
    print("Running under VSCode, load languageserver, showtext, httpgd")
    library(languageserver)
    library(showtext)
    library(httpgd)

    httpgd::hgd()
    options(vsc.use_httpgd = TRUE) # Use httpgd for plotting in VSCode
    httpgd::hgd_view() # Open the httpgd viewer pane in VSCode
    showtext::showtext_auto()
  }
}


# determine if R is running in RSTUDIO/VSCode/Positron
if (
  interactive() &&
    (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode"))
) {
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
print(paste0(
  "Detected ",
  available_cores,
  " available cores (LSF_DJOB_NUMPROC = ",
  lsf_cores,
  ")."
))

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
print(paste0(
  "Detected ",
  omp_threads,
  " OpenMP threads (OMP_NUM_THREADS = ",
  Sys.getenv("OMP_NUM_THREADS"),
  ")."
))

if (session_plan == "multisession") {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        16
      )
    )
  doFuture::registerDoFuture()
} else {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        32
      )
    )
  # multicore_workers <- MulticoreParam(
  #   workers = workers_2_use - 1,                # Number of allocated CPU cores
  #   progressbar = TRUE,         # Show visual progress bars
  #   stop.on.error = TRUE        # Halt execution if a core errors out
  # )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  options(useUCSCChromosomeNames = FALSE)
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
print(paste0(
  "R session plan set to ",
  session_plan,
  " with ",
  workers_2_use,
  " workers."
))

# ArchR settings ####
addArchRThreads(threads = workers_2_use)
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)

# addArchRThreads(threads = 1)
merged_output_dir <- "ArchR_merged_ATAC_multiome_obj"

projMerged <-
  loadArchRProject(
    path = merged_output_dir,
    showLogo = FALSE,
    force = TRUE
  )
projMerged <-
  ArchR::addIterativeLSI(
    ArchRProj = projMerged,
    useMatrix = "PeakMatrix",
    name = "IterativeLSI",
    firstSelection = "top",
    iterations = 4,
    clusterParams = list(
      # single resolution only: a vector makes Seurat::FindClusters return
      # multiple cluster columns, which breaks ArchR's .LSICluster extraction.
      resolution = 2,
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    sampleCellsPre = 25000,
    projectCellsPre = FALSE,
    sampleCellsFinal = 25000,
    seed = 42,
    saveIterations = FALSE,
    force = TRUE
  )


print(paste0(
  "Adding Harmony for ",
  length(getArrowFiles(projMerged)),
  " samples..."
))
print(
  "!! We can only use IterativeLSI (ATAC) reducedDims for Harmony clustering, tSNE, UMAP, and peak calling !!"
)
projMerged <-
  ArchR::addHarmony(
    ArchRProj = projMerged,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )
print("Harmony added and ArchRProject saved successfully.")


# try iteration approach to find the resolution whose cluster count brackets 10.
# --- Iterative resolution search for ~10 clusters -------------------------
# Start at resolution 0.5 and adjust by factors of 2 until the cluster count
# brackets 10. Each attempt is stored under its own str_c("Cluster_", res) name.
#   - if the count is < 10, double the resolution until the count exceeds 10
#   - if the count is > 10, halve the resolution until the count drops below 10
iter_resolution <- 0.5
iter_cluster_counts <- integer(0) # named: cluster_name -> n clusters
iter_max_steps <- 20 # safety cap to avoid an unbounded loop

add_clusters_at_resolution <-
  function(
    proj,
    resolution_value
  ) {
    cluster_name <- str_c("Cluster_", resolution_value)
    print(paste0(
      "Adding Clusters '",
      cluster_name,
      "' at resolution ",
      resolution_value,
      " using IterativeLSI reducedDims..."
    ))
    proj <-
      ArchR::addClusters(
        input = proj,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = cluster_name,
        resolution = resolution_value,
        force = TRUE
      )
    n_clusters <-
      length(unique(getCellColData(proj, select = cluster_name)[, 1]))
    list(proj = proj, cluster_name = cluster_name, n_clusters = n_clusters)
  }

# step 1: cluster at the starting resolution (0.5)
iter_result <- add_clusters_at_resolution(projMerged, iter_resolution)
projMerged <- iter_result$proj
iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters

# step 2/3: adjust the resolution until the cluster count brackets 10
iter_step <- 1
if (iter_result$n_clusters < 10) {
  # too few clusters -> keep doubling the resolution until count > 10
  while (iter_result$n_clusters < 10 && iter_step < iter_max_steps) {
    iter_resolution <- iter_resolution * 2
    iter_result <- add_clusters_at_resolution(projMerged, iter_resolution)
    projMerged <- iter_result$proj
    iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters
    iter_step <- iter_step + 1
  }
} else if (iter_result$n_clusters > 10) {
  # too many clusters -> keep halving the resolution until count < 10
  while (iter_result$n_clusters > 10 && iter_step < iter_max_steps) {
    iter_resolution <- iter_resolution / 2
    iter_result <- add_clusters_at_resolution(projMerged, iter_resolution)
    projMerged <- iter_result$proj
    iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters
    iter_step <- iter_step + 1
  }
}

# warn if the search hit the step cap without reaching ~10 clusters
if (iter_step >= iter_max_steps && iter_result$n_clusters != 10) {
  warning(paste0(
    "Resolution search hit the step cap (",
    iter_max_steps,
    ") without reaching 10 clusters. Last resolution = ",
    iter_resolution,
    " gave ",
    iter_result$n_clusters,
    " clusters."
  ))
}

projMerged <-
  saveArchRProject(
    ArchRProj = projMerged,
    outputDirectory = merged_output_dir,
    load = TRUE,
    overwrite = TRUE
  )
# selected_cells_obj <-
#   qs_read(
#     "annotation_script_package/cluster_based/manual_selection_hepatocytes.qs2",
#     nthreads = 8
#   )

# head(rownames(projMerged@cellColData))
# head(colnames(selected_cells_obj))

# barcodes_w_ident <-
#   colnames(selected_cells_obj)
# barcodes_w_ident <-
#   sub(
#     "^([^_]*__[^_]*)_",
#     "\\1#",
#     barcodes_w_ident
#   )
# head(barcodes_w_ident)

# unique(str_split(
#   rownames(projMerged@cellColData),
#   pattern = "#",
#   simplify = TRUE
# )[, 1])
#  [1] "X__pt14pre"  "X__pt24pre"  "X__pt17pre"  "X__pt16pre"  "X__24t"      "X__pt1_pre"  "X__pt12pre"
#  [8] "X__pt6_pre"  "X__pt5_pre"  "X__pt11_pre" "X__pt9_pre"  "X__13t"      "X__20t"      "X__3t"
# [15] "X__1t"       "X__22t"      "X__19t"      "X__18t"      "X__17t"      "X__15t"      "X__pt22pre"
# [22] "X__pt29pre"  "X__16t"      "X__11t"      "X__pt30pre"  "X__7t"       "X__12t"      "X__pt26pre"
# [29] "X__pt20pre"  "X__pt23pre"

# unique(str_split(
#   barcodes_w_ident,
#   pattern = "#",
#   simplify = TRUE
# )[, 1]) |>
#   str_split(
#     pattern = "_",
#     simplify = TRUE
#   ) |>
#   _[, 1] |>
#   unique()
#  [1] "BC001"      "BC002"      "BC003"      "BC004"      "BC005"      "BC006"      "BC007"
#  [8] "BC008"      "BC009"      "BC010"      "BC011"      "BC012"      "BC013"      "BC014"
# [15] "BC015"      "BC016"      "GEX_pt1pre" "X__11t"     "X__12t"     "X__13t"     "X__15t"
# [22] "X__16t"     "X__17t"     "X__18t"     "X__19t"     "X__1t"      "X__20t"     "X__22t"
# [29] "X__3t"      "X__7t"      "X__pt20pre" "X__pt22pre" "X__pt23pre" "X__pt26pre" "X__pt29pre"
# [36] "X__pt30pre"

# sum(rownames(projMerged@cellColData) %in% barcodes_w_ident)
# # # add transferred_barcode column to cellColData and fill it with NA
# projMerged@cellColData$transferred_barcode <- NA
# projMerged@cellColData$transferred_barcode[
#   rownames(projMerged@cellColData) %in% barcodes_w_ident
# ] <- "Hepatocytes"

# projMerged <-
#   saveArchRProject(
#     ArchRProj = projMerged,
#     outputDirectory = merged_output_dir,
#     load = TRUE,
#     overwrite = TRUE
#   )

cat(
  paste0(
    "Iterative resolution search cluster counts:\n",
    paste0(
      names(iter_cluster_counts),
      ": ",
      iter_cluster_counts,
      collapse = "\n"
    )
  ),
  "\n"
)

print(paste0(
  "Adding UMAP for ",
  length(ArchR::getArrowFiles(projMerged)),
  " samples using Harmony reducedDims..."
))
projMerged <-
  ArchR::addUMAP(
    ArchRProj = projMerged,
    reducedDims = "Harmony",
    name = "UMAP_RNA",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )

print("UMAP added and ArchRProject saved successfully.")

print(paste0(
  "Adding tSNE for ",
  length(ArchR::getArrowFiles(projMerged)),
  " samples using Harmony reducedDims..."
))
projMerged <-
  ArchR::addTSNE(
    ArchRProj = projMerged,
    reducedDims = "Harmony",
    name = "tSNE_RNA",
    perplexity = 30,
    force = TRUE
  )
projMerged <-
  saveArchRProject(
    ArchRProj = projMerged,
    outputDirectory = merged_output_dir,
    load = TRUE,
    overwrite = TRUE
  )
print("tSNE added and ArchRProject saved successfully.")

# print("Peak set copied from projPeakSource and peak matrix added successfully.")
