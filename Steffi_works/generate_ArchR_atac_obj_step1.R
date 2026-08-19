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
    (Sys.getenv("TERM_PROGRAM") == "vscode") && (Sys.getenv("POSITRON") != "1")
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
addArchRThreads(threads = workers_2_use)
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

# functions ####
create_granges_from_bed <-
  function(bed_file) {
    # The 4th column ("Low Mappability" / "High Signal Region") contains spaces, so
    # read as strictly tab-delimited to keep it as a single field.
    blacklist_df <-
      read.table(
        bed_file,
        header = FALSE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE,
        col.names = c("chrom", "start", "end", "reason")
      )

    blacklist_gr <-
      GenomicRanges::GRanges(
        seqnames = blacklist_df$chrom,
        ranges = IRanges::IRanges(
          start = blacklist_df$start + 1L, # BED 0-based -> GRanges 1-based
          end = blacklist_df$end
        ),
        strand = "*",
        reason = blacklist_df$reason
      )

    # Tag the genome build for downstream compatibility (e.g. ArchR / overlap ops).
    GenomeInfoDb::genome(blacklist_gr) <- "hg38"

    # Sort by chromosome then coordinate for tidy, reproducible output.
    blacklist_gr <- GenomeInfoDb::sortSeqlevels(blacklist_gr)
    blacklist_gr <- sort(blacklist_gr)

    return(blacklist_gr)
  }

setwd("/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works")

# find all fragments.tsv.gz files in the current working directory and its subdirectories
input_bams_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/atac_bams"
fragments_files <-
  list.files(
    path = input_bams_dir,
    pattern = "fragments.tsv.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
# Some sample directories contain an accidental duplicated nested
# "outs/outs" copy of the cellranger-atac output; drop those so each
# sample only contributes a single fragments.tsv.gz file.
fragments_files <-
  fragments_files[!str_detect(fragments_files, pattern = "/outs/outs/")]
names(fragments_files) <-
  paste0(
    "X__",
    str_replace(
      str_replace(
        fragments_files,
        pattern = paste0(input_bams_dir, "/"),
        replacement = ""
      ),
      pattern = "/outs/fragments.tsv.gz",
      replacement = ""
    )
  )

# find all barcodes.tsv files in the current working directory and its subdirectories
csv_barcode_files <-
  list.files(
    path = input_bams_dir,
    pattern = "barcodes.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
csv_barcode_files <-
  csv_barcode_files[str_detect(
    csv_barcode_files,
    pattern = "filtered_peak_bc_matrix"
  ) & !str_detect(csv_barcode_files, pattern = "/outs/outs/")]
names(csv_barcode_files) <-
  paste0(
    "X__",
    str_replace(
      str_replace(
        csv_barcode_files,
        pattern = paste0(input_bams_dir, "/"),
        replacement = ""
      ),
      pattern = "/outs/filtered_peak_bc_matrix/barcodes.tsv",
      replacement = ""
    )
  )

# find all h5 files in the current working directory and its subdirectories
h5_files <-
  list.files(
    path = input_bams_dir,
    pattern = "filtered_peak_bc_matrix.h5$",
    recursive = TRUE,
    full.names = TRUE
  )
h5_files <-
  h5_files[!str_detect(h5_files, pattern = "/outs/outs/")]
names(h5_files) <-
  paste0(
    "X__",
    str_replace(
      str_replace(
        h5_files,
        pattern = paste0(input_bams_dir, "/"),
        replacement = ""
      ),
      pattern = "/outs/filtered_peak_bc_matrix.h5",
      replacement = ""
    )
  )

print(paste0(
  "Found ",
  length(fragments_files),
  " fragments.tsv.gz files, ",
  length(csv_barcode_files),
  " barcodes.tsv files, and ",
  length(h5_files),
  " h5 files."
))

if (
  all(
    names(fragments_files) == names(csv_barcode_files) &
      names(fragments_files) == names(h5_files)
  )
) {
  print(
    "All fragments.tsv.gz files have corresponding barcodes.tsv files and h5 files."
  )
} else {
  stop("Mismatch between fragments.tsv.gz, barcodes.tsv, and h5 files.")
}

# Test: only keep the first 2 samples for testing
# fragments_files <- fragments_files[1:2]
# csv_barcode_files <- csv_barcode_files[1:2]
# h5_files <- h5_files[1:2]

# barcodes.tsv files are plain barcode lists (no header), not 10x
# singlecell.csv files, so ArchR::getValidBarcodes() (which requires
# cell_id/is__cell_barcode columns) does not apply here. Build the
# SimpleList of valid barcodes directly instead.
filtered_barcode_list <- S4Vectors::SimpleList()
# typeof(filtered_barcode_list@listData[[1]]) # should be "character"
# create a new barcode list and use it to replace the "listData" slot of the filtered_barcode_list object
insert_listData <-
  lapply(csv_barcode_files, function(csv_file) {
    df <- read.csv(
      csv_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )
    df <- unname(unlist(df)) # convert the list to a character vector
    return(df)
  })
names(insert_listData) <- as.character(names(csv_barcode_files))
filtered_barcode_list@listData <- insert_listData

# generate the gene and genome annotation objects for hg38
geneAnno <- createGeneAnnotation("hg38")

# generate the genome annotation object for hg38, with the blacklist GRanges object created from the ENCODE hg38 blacklist BED file
blacklist_bed_file <- "hg38-blacklist.v2.bed"
blacklist_gr <- create_granges_from_bed(blacklist_bed_file)
genomeAnno <-
  createGenomeAnnotation("hg38", blacklist = blacklist_gr)

print(paste0(
  "Creating Arrow files for ",
  length(fragments_files),
  " samples..."
))
# choose + create the destination directory for Arrow files so it will not save them in the current working directory
arrow_dir <- "atac_arrow_files"
if (!dir.exists(arrow_dir)) {
  dir.create(arrow_dir, recursive = TRUE)
  print(paste0("Created Arrow files directory: ", arrow_dir))
} else {
  print(paste0("Arrow files directory already exists: ", arrow_dir))
}

# create Arrow files for each sample, using the filtered barcodes and the gene and genome annotation objects
arrow_files <-
  ArchR::createArrowFiles(
    inputFiles = fragments_files,
    sampleNames = names(fragments_files),
    outputNames = file.path(arrow_dir, names(fragments_files)),
    minTSS = 1,
    minFrags = 1000,
    maxFrags = 1e+07,
    excludeChr = c("chrM", "chrX", "chrY"),
    force = TRUE,
    validBarcodes = filtered_barcode_list,
    threads = workers_2_use - 1,
    geneAnnotation = geneAnno,
    genomeAnnotation = genomeAnno,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )
names(arrow_files) <- names(fragments_files)

print(paste0("Creating ArchRProject for ", length(arrow_files), " samples..."))
projMultiome <-
  ArchR::ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = "ArchR_atac_obj",
    copyArrows = TRUE
  )
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("ArchRProject created and saved successfully.")


# filter out doublets from the ArchRProject
projMultiome <-
  addDoubletScores(
    input = projMultiome,
    useMatrix = "TileMatrix",
    k = 10, # Refers to how many similar cells should be used when constructing the artificial doublets
    knnMethod = "LSI", # Refers to the embedding to use for nearest neighbor search
    LSIMethod = 1,
    LSIParams = list(
      dimsToUse = 1:20, # try 5 and 20
      varFeatures = 1000,
      iterations = 5 # try 2 and 5
    ),
    force = TRUE
  )

# calculate doubletenrichment per sample at given cutEnrich threshold (2) and print the results
cutenrich <- 2
doublets_per_sample <-
  tapply(
    X = getCellColData(
      projMultiome,
      select = c("Sample", "DoubletEnrichment")
    )$DoubletEnrichment >=
      cutenrich,
    INDEX = getCellColData(
      projMultiome,
      select = c("Sample", "DoubletEnrichment")
    )$Sample,
    FUN = sum,
    na.rm = TRUE
  ) # doublets per sample at cutEnrich = 1
cat(
  paste0(
    "Doublets per sample at cutEnrich = ",
    cutenrich,
    ":\n",
    names(doublets_per_sample),
    ": ",
    as.integer(doublets_per_sample),
    collapse = "\n"
  ),
  "\n"
)

# calculate doubletscored per sample at given cutscore threshold (100) and print the results
cutscore <- 100
doublets_per_sample <-
  tapply(
    X = getCellColData(
      projMultiome,
      select = c("Sample", "DoubletScore")
    )$DoubletScore >=
      cutscore,
    INDEX = getCellColData(
      projMultiome,
      select = c("Sample", "DoubletScore")
    )$Sample,
    FUN = sum,
    na.rm = TRUE
  ) # doublets per sample at cutscore = 100
cat(
  paste0(
    "Doublets per sample at cutscore = ",
    cutscore,
    ":\n",
    names(doublets_per_sample),
    ": ",
    as.integer(doublets_per_sample),
    collapse = "\n"
  ),
  "\n"
)

# addDoubletScores() leaves DoubletEnrichment/DoubletScore as NA for samples it
# could not score (e.g. too few cells). filterDoublets() ignores those cells, so
# drop them here to keep only cells that received a valid doublet score.
cells_scored <-
  getCellNames(projMultiome)[
    (!is.na(getCellColData(projMultiome, select = "DoubletEnrichment")[, 1])) &
      (!is.na(getCellColData(projMultiome, select = "DoubletScore")[, 1]))
  ]
print(paste0(
  "Removing ",
  nCells(projMultiome) - length(cells_scored),
  " cells with NA doublet scores; ",
  length(cells_scored),
  " scored cells remain."
))
projMultiome <-
  subsetArchRProject(
    ArchRProj = projMultiome,
    cells = cells_scored,
    outputDirectory = "ArchR_atac_obj",
    force = TRUE
  )

apply_filtering_switch <- FALSE # set to TRUE to filter doublets, FALSE to skip filtering
if (apply_filtering_switch) {
  print(paste0(
    "Filtering doublets from ArchRProject for ",
    length(arrow_files),
    " samples..."
  ))
  projMultiome <-
    filterDoublets(
      ArchRProj = projMultiome,
      cutEnrich = cutenrich,
      cutScore = cutscore,
      filterRatio = 1.0
    )
  cat(paste0(
    "Doublets filtered from ArchRProject for ",
    length(arrow_files),
    " samples with remaining cell numbers: \n",
    {
      cells_per_sample <-
        table(getCellColData(projMultiome, select = "Sample")[, 1])
      paste0(
        names(cells_per_sample),
        ": ",
        as.integer(cells_per_sample),
        collapse = "\n"
      )
    },
    "\n"
  ))
} else {
  print("Skipping doublet filtering, keeping all cells in ArchRProject.")
}

print(paste0("Adding Iterative LSI for ", length(arrow_files), " samples..."))
projMultiome <-
  ArchR::addIterativeLSI(
    ArchRProj = projMultiome,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    firstSelection = "top",
    iterations = 2,
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
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("Iterative LSI added and ArchRProject saved successfully.")

print(paste0("Adding Harmony for ", length(arrow_files), " samples..."))
print(
  "!! We can only use IterativeLSI (ATAC) reducedDims for Harmony clustering, tSNE, UMAP, and peak calling !!"
)
projMultiome <-
  ArchR::addHarmony(
    ArchRProj = projMultiome,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )
projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("Harmony added and ArchRProject saved successfully.")

# ! Remember to use the "Harmony" reducedDims for downstream analysis, such as clustering and UMAP visualization.

# Cluster at a range of resolutions, storing each result under its own name
# (str_c("Cluster_", resolution_value)). Afterwards, select the resolution whose
# cluster count is closest to 10 and keep its name in cluster_name_4_peak_calling.
# resolution_values <- c(0.01, 0.02, 0.04, 0.08, 0.16)
# cluster_counts <- integer(length(resolution_values))
# names(cluster_counts) <- str_c("Cluster_", resolution_values)

# for (resolution_value in resolution_values) {
#   cluster_name <- str_c("Cluster_", resolution_value)
#   print(paste0(
#     "Adding Clusters '",
#     cluster_name,
#     "' for ",
#     length(arrow_files),
#     " samples at resolution ",
#     resolution_value,
#     " using Harmony reducedDims..."
#   ))
#   projMultiome <-
#     ArchR::addClusters(
#       input = projMultiome,
#       reducedDims = "Harmony",
#       method = "Seurat",
#       name = cluster_name,
#       resolution = resolution_value,
#       force = TRUE
#     )
#   cluster_counts[cluster_name] <-
#     length(unique(
#       getCellColData(projMultiome, select = cluster_name)[, 1]
#     ))
# }
# projMultiome <-
#   saveArchRProject(
#     ArchRProj = projMultiome,
    #  outputDirectory = "ArchR_atac_obj",
#     load = TRUE,
#     overwrite = TRUE
#   )
# cat(
#   paste0(
#     "Cluster counts per resolution:\n",
#     paste0(names(cluster_counts), ": ", cluster_counts, collapse = "\n")
#   ),
#   "\n"
# )
# # pick the resolution whose cluster count is nearest to 10 (ties -> lowest
# # resolution, i.e. first match) and store its name for downstream peak calling.
# cluster_name_4_peak_calling <-
#   names(cluster_counts)[which.min(abs(cluster_counts - 10))]
# print(paste0(
#   "Selected '",
#   cluster_name_4_peak_calling,
#   "' (",
#   cluster_counts[cluster_name_4_peak_calling],
#   " clusters) as nearest to 10 clusters for peak calling."
# ))

# try iteration approach to find the resolution whose cluster count brackets 10.
# --- Iterative resolution search for ~10 clusters -------------------------
# Start at resolution 0.5 and adjust by factors of 2 until the cluster count
# brackets 10. Each attempt is stored under its own str_c("Cluster_", res) name.
#   - if the count is < 10, double the resolution until the count exceeds 10
#   - if the count is > 10, halve the resolution until the count drops below 10
iter_resolution <- 0.5
iter_cluster_counts <- integer(0) # named: cluster_name -> n clusters
iter_max_steps <- 20 # safety cap to avoid an unbounded loop

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

# step 1: cluster at the starting resolution (0.5)
iter_result <- add_clusters_at_resolution(projMultiome, iter_resolution)
projMultiome <- iter_result$proj
iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters

# step 2/3: adjust the resolution until the cluster count brackets 10
iter_step <- 1
if (iter_result$n_clusters < 10) {
  # too few clusters -> keep doubling the resolution until count > 10
  while (iter_result$n_clusters < 10 && iter_step < iter_max_steps) {
    iter_resolution <- iter_resolution * 2
    iter_result <- add_clusters_at_resolution(projMultiome, iter_resolution)
    projMultiome <- iter_result$proj
    iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters
    iter_step <- iter_step + 1
  }
} else if (iter_result$n_clusters > 10) {
  # too many clusters -> keep halving the resolution until count < 10
  while (iter_result$n_clusters > 10 && iter_step < iter_max_steps) {
    iter_resolution <- iter_resolution / 2
    iter_result <- add_clusters_at_resolution(projMultiome, iter_resolution)
    projMultiome <- iter_result$proj
    iter_cluster_counts[iter_result$cluster_name] <- iter_result$n_clusters
    iter_step <- iter_step + 1
  }
}

projMultiome <-
  saveArchRProject(
    ArchRProj = projMultiome,
    outputDirectory = "ArchR_atac_obj",
    load = TRUE,
    overwrite = TRUE
  )

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

# optimum resolution = the attempt whose cluster count is nearest to 10.
# NOTE: the grid-based block below re-assigns cluster_name_4_peak_calling; keep
# only the approach you actually want to drive peak calling.
# cluster_name_4_peak_calling <-
#   names(iter_cluster_counts)[which.min(abs(iter_cluster_counts - 10))]
# print(paste0(
#   "Iterative search selected '",
#   cluster_name_4_peak_calling,
#   "' (",
#   iter_cluster_counts[cluster_name_4_peak_calling],
#   " clusters) as nearest to 10 clusters for peak calling."
# ))
# # --------------------------------------------------------------------------
# print("Clusters added and ArchRProject saved successfully.")

# print(paste0(
#   "Adding UMAP for ",
#   length(arrow_files),
#   " samples using Harmony reducedDims..."
# ))
# projMultiome <-
#   ArchR::addUMAP(
#     ArchRProj = projMultiome,
#     reducedDims = "Harmony",
#     name = "UMAP_RNA",
#     nNeighbors = 30,
#     minDist = 0.5,
#     metric = "cosine",
#     force = TRUE
#   )
# projMultiome <-
#   saveArchRProject(
#     ArchRProj = projMultiome,
#     outputDirectory = "ArchR_atac_obj",
#     load = TRUE,
#     overwrite = TRUE
#   )
# print("UMAP added and ArchRProject saved successfully.")

# print(paste0(
#   "Adding tSNE for ",
#   length(arrow_files),
#   " samples using Harmony reducedDims..."
# ))
# projMultiome <-
#   ArchR::addTSNE(
#     ArchRProj = projMultiome,
#     reducedDims = "Harmony",
#     name = "tSNE_RNA",
#     perplexity = 30,
#     force = TRUE
#   )
# projMultiome <-
#   saveArchRProject(
#     ArchRProj = projMultiome,
#     outputDirectory = "ArchR_atac_obj",
#     load = TRUE,
#     overwrite = TRUE
#   )
# print("tSNE added and ArchRProject saved successfully.")

# # plotEmbedding(
# #   ArchRProj = projMultiome,
# #   colorBy = "cellColData",
# #   name = cluster_name_4_peak_calling,
# #   embedding = "UMAP"
# # )

# # Prepare for peak calling by adding a pseudo-bulk replicate for each cluster, using the cluster_name_4_peak_calling determined above.
# projMultiome <-
#   ArchR::addGroupCoverages(
#     ArchRProj = projMultiome,
#     groupBy = cluster_name_4_peak_calling,
#     excludeChr = c("chrM", "chrX", "chrY"),
#     force = TRUE
#   )

# pathToMacs2 <- findMacs2()
# projMultiome <-
#   ArchR::addReproduciblePeakSet(
#     ArchRProj = projMultiome,
#     groupBy = cluster_name_4_peak_calling,
#     pathToMacs2 = pathToMacs2,
#     excludeChr = c("chrM", "chrX", "chrY"),
#     force = TRUE
#   )
# projMultiome <-
#   saveArchRProject(
#     ArchRProj = projMultiome,
#     outputDirectory = "ArchR_atac_obj",
#     load = TRUE,
#     overwrite = TRUE
#   )
# print("Reproducible peak set added and ArchRProject saved successfully.")
