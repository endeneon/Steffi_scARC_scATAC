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

  if (interactive() &&
    ((Sys.getenv("TERM_PROGRAM") == "vscode") && (Sys.getenv("POSITRON") != "1"))
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
if (interactive() &&
  (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode"))) {
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

setwd("/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works")

projPeakSource <-
  loadArchRProject(
    path = "ArchR_multiome_obj",
    showLogo = FALSE,
    force = TRUE
  )


projMultiome <-
  loadArchRProject(
    path = "ArchR_atac_obj",
    showLogo = FALSE,
    force = TRUE
  )

iter_resolution <- 0.03125
iter_cluster_counts <- integer(0) # named: cluster_name -> n clusters
# iter_max_steps <- 20 # safety cap to avoid an unbounded loop

add_clusters_at_resolution <- function(proj, resolution_value) {
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

# step 1: cluster at the pre-determined resolution (0.03125)
iter_result <- add_clusters_at_resolution(projMultiome, iter_resolution)
projMultiome <- iter_result$proj
iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters

# optimum resolution = the attempt whose cluster count is nearest to 10.
# NOTE: the grid-based block below re-assigns cluster_name_4_peak_calling; keep
# only the approach you actually want to drive peak calling.
cluster_name_4_peak_calling <-
  names(iter_cluster_counts)[which.min(abs(iter_cluster_counts - 10))]
print(paste0(
  "Iterative search selected '",
  cluster_name_4_peak_calling,
  "' (",
  iter_cluster_counts[cluster_name_4_peak_calling],
  " clusters) as nearest to 10 clusters for peak calling."
))
# --------------------------------------------------------------------------
# print("Clusters added and ArchRProject saved successfully.")

print(paste0(
  "Adding UMAP for ",
  length(ArchR::getArrowFiles(projMultiome)),
  " samples using Harmony reducedDims..."
))
projMultiome <-
  ArchR::addUMAP(
    ArchRProj = projMultiome,
    reducedDims = "Harmony",
    name = "UMAP_RNA",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("UMAP added and ArchRProject saved successfully.")

print(paste0(
  "Adding tSNE for ",
  length(ArchR::getArrowFiles(projMultiome)),
  " samples using Harmony reducedDims..."
))
projMultiome <-
  ArchR::addTSNE(
    ArchRProj = projMultiome,
    reducedDims = "Harmony",
    name = "tSNE_RNA",
    perplexity = 30,
    force = TRUE
  )
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("tSNE added and ArchRProject saved successfully.")

# plotEmbedding(
#   ArchRProj = projMultiome,
#   colorBy = "cellColData",
#   name = cluster_name_4_peak_calling,
#   embedding = "UMAP"
# )

# Reuse the peak set already called on the multiome project (projPeakSource)
# instead of re-calling peaks with MACS2 on the ATAC-only project.
projMultiome <-
  ArchR::addPeakSet(
    ArchRProj = projMultiome,
    peakSet = ArchR::getPeakSet(projPeakSource),
    force = TRUE
  )
projMultiome <-
  ArchR::addPeakMatrix(
    ArchRProj = projMultiome,
    force = TRUE
  )
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("Peak set copied from projPeakSource and peak matrix added successfully.")
