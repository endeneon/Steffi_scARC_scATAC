#! /usr/bin/env Rscript

# init
{
  library(Seurat)
  library(gplots)
  library(ArchR)
  library(future)
  library(stringr)
  library(scCustomize)
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

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main"
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
addArchRThreads(threads = workers_2_use)
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

cluster_name_4_peak_calling <- "Cluster_0.03125"

# load ArchR object
projMultiome <-
  loadArchRProject(
    path = "Steffi_works/ArchR_multiome_obj",
    showLogo = FALSE
  )

# extract GeneExpressionMatrix + clustering and build a Seurat object ####

names(projMultiome@embeddings)
# reduced dimensions (IterativeLSI, Harmony, etc.)
names(projMultiome@reducedDims)
# [1] "UMAP_RNA" "tSNE_RNA"
# all available matrices (GeneExpressionMatrix, PeakMatrix, ...)
getAvailableMatrices(projMultiome)
# all per-cell metadata columns (every Cluster_* assignment, sample, QC)
colnames(getCellColData(projMultiome))

# GeneExpressionMatrix comes back as a SummarizedExperiment:
#   assay()   -> genes x cells sparse counts
#   rowData() -> gene metadata (symbol lives in $name, NOT the rownames)
#   colnames()-> ArchR cell names (match projMultiome cellNames)
gem_se <-
  getMatrixFromProject(
    ArchRProj = projMultiome,
    useMatrix = "GeneExpressionMatrix",
    verbose = TRUE,
    logFile = createLogFile("getGeneExpressionMatrix")
  )

# build the count matrix with proper gene + cell names
gex_counts <- assay(gem_se) # dgCMatrix, genes x cells
rownames(gex_counts) <- rowData(gem_se)$name # gene symbols
colnames(gex_counts) <- colnames(gem_se) # ArchR cell names

# ArchR gene rows can carry NA / duplicated symbols -> make them safe
gex_counts <- gex_counts[!is.na(rownames(gex_counts)), ]
rownames(gex_counts) <- make.unique(as.character(rownames(gex_counts)))

# pull all per-cell metadata (every Cluster_* column, sample, QC, etc.)
cell_meta <- as.data.frame(getCellColData(projMultiome))
cell_meta <- cell_meta[colnames(gex_counts), , drop = FALSE]

# create the Seurat object (GeneExpressionMatrix stores raw counts)
seurat_gex <-
  CreateSeuratObject(
    counts = gex_counts,
    assay = "RNA",
    meta.data = cell_meta,
    project = "Multiome_GEX"
  )

# use the peak-calling clustering column as the active identity
Idents(seurat_gex) <- cluster_name_4_peak_calling

# transfer ArchR dimensionality reductions as Seurat DimReduc objects ####

# helper: copy one ArchR embedding into the Seurat object, aligning cells
add_archr_embedding <- function(
  seurat_obj,
  archr_proj,
  embedding_name,
  reduction_key
) {
  if (!embedding_name %in% names(archr_proj@embeddings)) {
    message(paste0("Embedding '", embedding_name, "' not found, skipping."))
    return(seurat_obj)
  }

  emb <- getEmbedding(
    ArchRProj = archr_proj,
    embedding = embedding_name,
    returnDF = TRUE
  )
  emb <- as.matrix(emb)

  # align embedding rows to the Seurat cells and drop any missing ones
  common_cells <- intersect(colnames(seurat_obj), rownames(emb))
  emb <- emb[common_cells, , drop = FALSE]
  colnames(emb) <- paste0(reduction_key, seq_len(ncol(emb)))

  seurat_obj[[embedding_name]] <-
    CreateDimReducObject(
      embeddings = emb,
      key = reduction_key,
      assay = DefaultAssay(seurat_obj)
    )
  seurat_obj
}

# inspect what reductions exist in the ArchR project, then copy them over
print(names(projMultiome@embeddings))

seurat_gex <- add_archr_embedding(
  seurat_obj = seurat_gex,
  archr_proj = projMultiome,
  embedding_name = "UMAP_RNA",
  reduction_key = "UMAPRNA_"
)
seurat_gex <- add_archr_embedding(
  seurat_obj = seurat_gex,
  archr_proj = projMultiome,
  embedding_name = "tSNE_RNA",
  reduction_key = "tSNERNA_"
)

# copy the IterativeLSI (or Harmony) reduction if present
if ("IterativeLSI" %in% names(projMultiome@reducedDims)) {
  lsi <- getReducedDims(projMultiome, reducedDims = "IterativeLSI")
  lsi <- lsi[intersect(colnames(seurat_gex), rownames(lsi)), , drop = FALSE]
  colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))
  seurat_gex[["IterativeLSI"]] <-
    CreateDimReducObject(
      embeddings = lsi,
      key = "LSI_",
      assay = DefaultAssay(seurat_gex)
    )
}

qs_save(seurat_gex,
  file = "Steffi_works/seurat_gex_4_identity_projection.qs2",
  nthreads = 4
)

seurat_gex <- UpdateSeuratObject(seurat_gex)
seurat_gex$orig.ident <- 
str_split(
    seurat_gex$orig.ident, "#", simplify = TRUE)[, 1]

library(SingleR)
# [22] "DoubletEnrichment" "Cluster_0.5"       "Cluster_0.25"     
# [25] "Cluster_0.125"     "Cluster_0.0625"    "Cluster_0.03125" 
unique(seurat_gex$Cluster_0.0625)
pseudobulk_seurat_mtx <-
AggregateExpression(
  seurat_gex,
  group.by = "Cluster_0.125",
  assays = DefaultAssay(seurat_gex),
  slot = "counts"
)
pseudobulk_seurat_mtx <-
pseudobulk_seurat_mtx$RNA
pseudobulk_seurat_mtx <-
edgeR::cpm(pseudobulk_seurat_mtx, log = TRUE)

ref_human_liver_atlas <-
Read10X(data.dir = "Steffi_works/liver_atlas_ref/rawData_human/countTable_human",
gene.column = 1)
ref_human_liver_annot <-
read.csv("Steffi_works/liver_atlas_ref/annot_humanAll.csv")
sum(ref_human_liver_annot$cell %in% colnames(ref_human_liver_atlas))
length(colnames(ref_human_liver_atlas) %in% ref_human_liver_annot$cell)
ref_human_liver_atlas <- ref_human_liver_atlas[, colnames(ref_human_liver_atlas) %in% ref_human_liver_annot$cell]
ref_human_liver_atlas <-
ref_human_liver_atlas[, order(colnames(ref_human_liver_atlas))]
ref_human_liver_annot <- ref_human_liver_annot[order(ref_human_liver_annot$cell), ]
all(colnames(ref_human_liver_atlas) == ref_human_liver_annot$cell)
ref_human_liver_atlas <- CreateSeuratObject(counts = ref_human_liver_atlas, project = "HumanLiverAtlas", meta.data = ref_human_liver_annot)
qs_save(ref_human_liver_atlas,
  file = "Steffi_works/ref_human_liver_atlas.qs2",
  nthreads = 4
)
pseudobulk_ref_mtx <-
AggregateExpression(
  ref_human_liver_atlas,
  group.by = "annot",
  assays = DefaultAssay(ref_human_liver_atlas),
  slot = "counts"
)
pseudobulk_ref_mtx <-
pseudobulk_ref_mtx$RNA
pseudobulk_ref_mtx <-
edgeR::cpm(pseudobulk_ref_mtx, log = TRUE)

predicted_annot <-
SingleR(
  test = pseudobulk_seurat_mtx,
  ref = pseudobulk_ref_mtx,
  labels = colnames(pseudobulk_ref_mtx),
  genes = "sd"
)
predicted_annot$pruned.labels
predicted_annot$Cluster_index <- rownames(predicted_annot)
seurat_gex$predicted.annotation <-
predicted_annot$pruned.labels[match(seurat_gex$Cluster_0.125, predicted_annot$Cluster_index)]
table(seurat_gex$predicted.annotation)
qs_save(seurat_gex,
  file = "Steffi_works/seurat_gex_with_predicted_annotation.qs2",
  nthreads = 4
)
DimPlot_scCustom(seurat_gex, 
group.by = "predicted.annotation",
colors_use = "glasbey")

DimPlot(seurat_gex, group.by = "predicted.annotation")
