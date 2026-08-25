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

  library(ggplot2)
  library(Gviz)

  if (
    interactive() &&
      (Sys.getenv("TERM_PROGRAM") == "vscode") &&
      (Sys.getenv("POSITRON") != "1")
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


setwd(
  "Steffi_works"
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
if (
  interactive() &&
    (Sys.getenv("TERM_PROGRAM") == "vscode")
) {
  print("Running under IDE, use 1 ArchR Thread")
  addArchRThreads(threads = 1)
} else {
  print("Running under Rscript, use all usable ArchR Threads")
  addArchRThreads(threads = workers_2_use)
}
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

projMerged <-
  ArchR::loadArchRProject(path = "ArchR_merged_ATAC_multiome_obj")
colnames(projMerged@cellColData)
unique(projMerged@cellColData$transferred_barcode)

# transfer barcodes with PeakMatrix UMAP embedding
plotEmbedding(
  ArchRProj = projMerged,
  embedding = "UMAP_Peaks",
  colorBy = "cellColData",
  name = "transferred_barcode",
  size = 0.1,
  # sampleCells = 10000,
  highlightCells = NULL,
  rastr = FALSE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10
)

# 2D density map of Hepatocytes in UMAP_Peaks space ####
# UMAP_Peaks embedding coordinates for all cells
umap_df <- ArchR::getEmbedding(
  ArchRProj = projMerged,
  embedding = "UMAP_Peaks",
  returnDF = TRUE
)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# align cellColData and subset to Hepatocytes
cell_meta <- as.data.frame(projMerged@cellColData)
umap_df$transferred_barcode <- cell_meta[
  rownames(umap_df),
  "transferred_barcode"
]

hep_df <- umap_df[
  !is.na(umap_df$transferred_barcode) &
    umap_df$transferred_barcode == "Hepatocytes",
]

lower_density_cutoff <- 0.005
density_plot <-
  ggplot(
    hep_df,
    aes(x = UMAP1, y = UMAP2)
  ) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    breaks = seq(
      lower_density_cutoff,
      0.1,
      by = lower_density_cutoff
    ), # lower density cut-off at 0.008
    alpha = 0.5
  ) +
  scale_fill_viridis_c(option = "magma") +
  geom_point(size = 0.1, alpha = 0.2, color = "grey20") +
  labs(
    title = "UMAP_Peaks 2D density: Hepatocytes",
    x = "UMAP 1",
    y = "UMAP 2",
    fill = "Density"
  ) +
  theme_minimal()

density_plot

projMerged@cellColData$peaks_projected_hepatocytes <- FALSE

# --- Project the Hepatocyte density coverage onto all cells ---------------
# Rebuild the same 2D kernel density that stat_density_2d draws (MASS::kde2d
# with per-axis bandwidth.nrd, matching ggplot's defaults), then flag every
# cell whose UMAP_Peaks position falls inside the plotted coverage, i.e. where
# the fitted density is >= lower_density_cutoff (the plot's lowest break).

# per-axis bandwidth (same as ggplot's stat_density_2d default)
kde_bandwidth <- c(
  MASS::bandwidth.nrd(hep_df$UMAP1),
  MASS::bandwidth.nrd(hep_df$UMAP2)
)

# fit the density on Hepatocyte cells over a grid spanning ALL cells so that
# non-Hepatocyte barcodes inside the region can also be evaluated
kde_fit <- MASS::kde2d(
  x = hep_df$UMAP1,
  y = hep_df$UMAP2,
  h = kde_bandwidth,
  n = 200,
  lims = c(range(umap_df$UMAP1), range(umap_df$UMAP2))
)

# look up the fitted density at each cell via its nearest grid cell
ix <- findInterval(umap_df$UMAP1, kde_fit$x, all.inside = TRUE)
iy <- findInterval(umap_df$UMAP2, kde_fit$y, all.inside = TRUE)
cell_density <- kde_fit$z[cbind(ix, iy)]

# barcodes inside the density coverage -> TRUE
inside_barcodes <- rownames(umap_df)[
  !is.na(cell_density) & cell_density >= lower_density_cutoff
]

coverage_flag <- rownames(projMerged@cellColData) %in% inside_barcodes
projMerged@cellColData$peaks_projected_hepatocytes <- coverage_flag

table(projMerged@cellColData$peaks_projected_hepatocytes)
plotEmbedding(
  ArchRProj = projMerged,
  embedding = "UMAP_Peaks",
  colorBy = "cellColData",
  name = "peaks_projected_hepatocytes",
  size = 0.1,
  # sampleCells = 10000,
  highlightCells = NULL,
  rastr = FALSE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10
)

projMerged <-
  saveArchRProject(
    ArchRProj = projMerged,
    outputDirectory = "ArchR_merged_ATAC_multiome_obj",
    load = TRUE,
    overwrite = TRUE
  )
