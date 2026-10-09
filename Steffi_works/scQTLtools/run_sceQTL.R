#! /usr/bin/env Rscript

# use conda env r_45_python_312 ####

# init
{
  library(Seurat)

  library(stringr)
  # library(pheatmap)

  # library(BiocParallel)
  # library(BiocParallel.FutureParam)
  # library(parallel)

  library(Matrix)
  library(matrixStats)

  library(qs2)
  library(ggplot2)

  library(scQTLtools)

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


# determine if R is running in RSTUDIO/VSCode/Positron and assign session directories
if (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode")) {
  print("Running under RStudio/VSCode/Positron IDE, use plan(multisession)")
  session_plan <- "multisession"
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  print("Running under Rscript, use plan(multicore)")
  session_plan <- "multicore"
  setwd(
    "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/scQTLtools"
  )
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

cat(c("Starting analysis at", as.character(Sys.time()), "\n"))
# # ArchR settings ####
# if (
#   interactive() &&
#     (Sys.getenv("TERM_PROGRAM") == "vscode")
# ) {
#   print("Running under IDE, use 1 ArchR Thread")
#   addArchRThreads(threads = 1)
# } else {
#   print("Running under Rscript, use all usable ArchR Threads")
#   addArchRThreads(threads = workers_2_use)
# }
# addArchRGenome("hg38")

# print(paste0(
#   "ArchR threads set to ",
#   getArchRThreads(),
#   " and genome set to ",
#   getArchRGenome()
# ))
print("All settings initialized successfully.")

# run mode ####
# Usage:
#   Rscript run_sceQTL.R prep          -> split data per cell type, write cell_types.txt
#   Rscript run_sceQTL.R run [index]   -> process one cell type (index defaults to $LSB_JOBINDEX)
cli_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(cli_args) >= 1) cli_args[1] else "prep"
prep_dir <- "per_celltype_input"
cell_type_file <- file.path(prep_dir, "cell_types.txt")

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

if (run_mode == "prep") {
  # prep step: runs once, shared by all array tasks ####
  dir.create(prep_dir, showWarnings = FALSE)

  seurat_multiome_obj <-
    qs2::qs_read(
      "../annotation_script_package/marker_based/annotated_seurat_marker_based.qs2",
      nthreads = 8
    )
  df_genotype_raw <-
    read.table(
      "Macrophage_ASoC_genotypes_extended_GT_only.tsv",
      header = TRUE,
      sep = "\t"
    )
  SeuratObject::Idents(seurat_multiome_obj) <- "preparation"
  seurat_subsetted_obj <-
    subset(
      seurat_multiome_obj,
      cells = SeuratObject::WhichCells(
        seurat_multiome_obj,
        idents = "multiome"
      )
    )
  rm(seurat_multiome_obj)
  gc()

  SeuratObject::Idents(seurat_subsetted_obj) <- "celltype_broad"

  # Genotype recoding does not depend on cell type, so do it once here.
  df_genotype_multiome <-
    df_genotype_raw[, stringr::str_detect(
      colnames(df_genotype_raw),
      "multiome"
    )]
  rownames(df_genotype_multiome) <-
    as.vector(stringr::str_c(
      df_genotype_raw[, 1],
      df_genotype_raw[, 2],
      df_genotype_raw[, 3],
      df_genotype_raw[, 4],
      sep = ":"
    ))
  colnames(df_genotype_multiome) <-
    stringr::str_replace(
      colnames(df_genotype_multiome),
      "multiome_",
      "X__"
    )
  colnames(df_genotype_multiome) <-
    stringr::str_replace(
      colnames(df_genotype_multiome),
      "_GT",
      ""
    )

  # Recode genotype strings to integer codes; anything unmatched (incl. NA) -> 0.
  df_genotype_multiome[] <- lapply(
    df_genotype_multiome,
    function(col) {
      dplyr::case_when(
        col == "0/0" ~ 1L,
        col %in% c("0/1", "1/0") ~ 3L,
        col == "1/1" ~ 2L,
        .default = 0L
      )
    }
  )
  qs2::qs_save(
    df_genotype_multiome,
    file = file.path(prep_dir, "genotype_multiome_recoded.qs2"),
    nthreads = 8
  )

  # Skip cell types too small for a reliable per-cell ZINB fit.
  # Override with env var MIN_CELLS_PER_TYPE.
  min_cells_per_type <- as.integer(Sys.getenv("MIN_CELLS_PER_TYPE", "100"))
  cell_type_counts <- table(as.character(seurat_subsetted_obj$celltype_broad))
  skipped_types <- names(cell_type_counts)[
    cell_type_counts < min_cells_per_type
  ]
  if (length(skipped_types) > 0) {
    cat(
      "Skipping cell types with <",
      min_cells_per_type,
      "cells:",
      paste0(
        skipped_types,
        " (",
        cell_type_counts[skipped_types],
        ")",
        collapse = ", "
      ),
      "\n"
    )
  }
  cell_type_list <-
    as.character(unique(seurat_subsetted_obj$celltype_broad))
  cell_type_list <- setdiff(cell_type_list, skipped_types)

  for (each_cell_type in cell_type_list) {
    cat("Saving subset for cell type:", each_cell_type, "\n")
    qs2::qs_save(
      subset(
        seurat_subsetted_obj,
        cells = SeuratObject::WhichCells(
          seurat_subsetted_obj,
          idents = each_cell_type
        )
      ),
      file = file.path(
        prep_dir,
        paste0("seurat_", safe_name(each_cell_type), ".qs2")
      ),
      nthreads = 8
    )
  }

  # One line per array task; line i is processed by LSB_JOBINDEX = i.
  writeLines(cell_type_list, cell_type_file)
  cat(c("Prep finished at", as.character(Sys.time()), "\n"))
} else if (run_mode == "run") {
  # per-cell-type step: one LSF array task each ####
  task_index <- as.integer(
    if (length(cli_args) >= 2) cli_args[2] else Sys.getenv("LSB_JOBINDEX", NA)
  )
  cell_type_list <- readLines(cell_type_file)
  if (
    is.na(task_index) || task_index < 1 || task_index > length(cell_type_list)
  ) {
    stop("Invalid task index: ", task_index)
  }
  each_cell_type <- cell_type_list[task_index]

  cat(
    c(
      "Running eQTL analysis for cell type:",
      each_cell_type,
      "(task",
      task_index,
      ") at",
      as.character(Sys.time()),
      "\n"
    )
  )
  seurat_celltype_obj <-
    qs2::qs_read(
      file.path(prep_dir, paste0("seurat_", safe_name(each_cell_type), ".qs2")),
      nthreads = 8
    )
  df_genotype_multiome <-
    qs2::qs_read(
      file.path(prep_dir, "genotype_multiome_recoded.qs2"),
      nthreads = 8
    )

  df_full_genotype <-
    as.data.frame(matrix(
      0,
      nrow = nrow(df_genotype_multiome),
      ncol = ncol(seurat_celltype_obj)
    ))
  rownames(df_full_genotype) <- rownames(df_genotype_multiome)
  colnames(df_full_genotype) <- colnames(seurat_celltype_obj)

  # Fill each cell column from its sample's genotype column. The sample prefix is
  # everything before the final "_<barcode>" segment of the cell name. Cells whose
  # sample has no genotype column (e.g. X__20t) are dropped from df_full_genotype.
  cell_sample <- sub("_[^_]+$", "", colnames(df_full_genotype))
  keep <- cell_sample %in% colnames(df_genotype_multiome)
  df_full_genotype <- df_full_genotype[, keep, drop = FALSE]
  df_full_genotype[] <- df_genotype_multiome[, cell_sample[keep], drop = FALSE]

  snp_matrix_raw <-
    as.matrix(df_full_genotype)

  cat(
    "Starting eQTL analysis with scQTLtools...",
    as.character(Sys.time()),
    "\n"
  )
  eqtl_seurat <-
    scQTLtools::createQTLObject(
      snpMatrix = snp_matrix_raw,
      genedata = seurat_celltype_obj,
      biClassify = FALSE,
      species = 'human',
      group = "category"
    )

  cat("Normalizing gene expression data...\n")
  eqtl_seurat <-
    scQTLtools::normalizeGene(
      eQTLObject = eqtl_seurat,
      method = "logNormalize"
    )
  cat("Filtering genes and SNPs...\n")
  eqtl_seurat_filtered <- scQTLtools::filterGeneSNP(
    eQTLObject = eqtl_seurat,
    snpNumOfCellsPercent = 2,
    expressionMin = 0,
    expressionNumOfCellsPercent = 2
  )
  # Offline gene locations ####
  # The firewall blocks Ensembl/biomaRt, which callQTL() uses (via createGeneLoc)
  # when upstream/downstream are set. Replace createGeneLoc with a local EnsDb
  # lookup. Chromosomes are returned with a "chr" prefix so they match the SNP IDs
  # ("chr1:pos:ref:alt"), which scQTLtools parses into chr_name = "chr1".
  # SNP positions need no lookup because the IDs are not rsIDs.
  createGeneLoc_offline <-
    function(geneList, gene_mart = NULL, geneDataset = NULL, OrgDb = NULL) {
      geneList <- unique(geneList)
      id_col <- if (grepl("^ENSG", geneList[[1]][1])) "gene_id" else "gene_name"
      genes_gr <-
        ensembldb::genes(
          EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
          filter = AnnotationFilter::SeqNameFilter(c(1:22)),
          columns = c("gene_id", "gene_name")
        )
      genes_gr <- genes_gr[S4Vectors::mcols(genes_gr)[[id_col]] %in% geneList]
      gene_loc <-
        data.frame(
          gene = S4Vectors::mcols(genes_gr)[[id_col]],
          chromosome_name = paste0(
            "chr",
            as.character(GenomicRanges::seqnames(genes_gr))
          ),
          start_position = GenomicRanges::start(genes_gr),
          end_position = GenomicRanges::end(genes_gr)
        )
      # duplicated gene symbols: keep the first locus
      gene_loc <- gene_loc[!duplicated(gene_loc$gene), ]
      colnames(gene_loc)[1] <-
        if (id_col == "gene_id") "ensembl_gene_id" else "external_gene_name"
      rownames(gene_loc) <- NULL
      gene_loc
    }
  assignInNamespace("createGeneLoc", createGeneLoc_offline, ns = "scQTLtools")

  source("callQTL_parallel.R")
  eqtl_output <-
    callQTL_parallel(
      eQTLObject = eqtl_seurat_filtered,
      gene_ids = NULL,
      downstream = -1e6,
      upstream = 1e6,
      gene_mart = NULL,
      snp_mart = NULL,
      pAdjustMethod = "BH",
      useModel = "zinb",
      pAdjustThreshold = 0.05,
      logfcThreshold = 0.1,
      n_cores = workers_2_use
    )
  cat(
    "eQTL analysis for cell type:",
    each_cell_type,
    "completed. Saving results...\n"
  )
  qs2::qs_save(
    eqtl_output,
    file = paste0(
      "seurat_QTL_output_",
      safe_name(each_cell_type),
      "_multiome.qs2"
    ),
    nthreads = 8
  )
  cat(c("Run finished at", as.character(Sys.time()), "\n"))
} else {
  stop("Unknown mode '", run_mode, "'. Use 'prep' or 'run [index]'.")
}
