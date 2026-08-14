# ============================================================================
# 01_fetch_reference.R — Download and build the Lu et al. 2022 HCC reference
# ============================================================================
#
# Step 1 of the reference-based annotation pipeline.
#
# Reference dataset:
#   Lu Y, Yang A, Quan C, et al. "A single-cell atlas of the multicellular
#   ecosystem of primary and metastatic hepatocellular carcinoma."
#   Nature Communications 2022;13:4594. doi:10.1038/s41467-022-32283-3
#   GEO: GSE149614 — 71,915 cells, 10 patients, 4 tissue sites
#   (non-tumor liver / primary tumor / PVTT / metastatic lymph node).
#
# Files (verified 2026-08-12):
#   GSE149614_HCC.metadata.updated.txt.gz      478 KB   cell annotations
#     columns: Cell, sample, res.3, site, patient, stage, virus, celltype
#     celltype levels: T/NK (25,591), Hepatocyte (20,782), Myeloid (15,947),
#                      B (3,685), Endothelial (3,644), Fibroblast (2,266)
#   GSE149614_HCC.scRNAseq.S71915.count.txt.gz 158 MB   genes x cells counts
#     (tab-separated; first row = cell barcodes matching metadata$Cell)
#
# Usage:
#   Rscript 01_fetch_reference.R <out_dir>
#   # or source() and call fetch_reference() / build_reference()
setwd("/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/annotation_script_package/reference_based")
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleR)
  library(SingleCellExperiment)
  library(data.table)
})

# Threads granted by the scheduler (exported as NTHREADS in the bsub script).
nthreads <- as.integer(Sys.getenv("NTHREADS", unset = "8"))
if (is.na(nthreads) || nthreads < 1L) nthreads <- 8L
setDTthreads(nthreads)

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl"
META_URL  <- file.path(GEO_BASE, "GSE149614_HCC.metadata.updated.txt.gz")
COUNT_URL <- file.path(GEO_BASE, "GSE149614_HCC.scRNAseq.S71915.count.txt.gz")

#' Download the reference files from GEO
#'
#' @param out_dir Directory to store downloaded files
#' @return list(meta_path, count_path)
fetch_reference <- function(out_dir = "reference_data") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  meta_path  <- file.path(out_dir, "GSE149614_HCC.metadata.updated.txt.gz")
  count_path <- file.path(out_dir, "GSE149614_HCC.scRNAseq.S71915.count.txt.gz")

  if (!file.exists(meta_path)) {
    message("Downloading metadata (478 KB) ...")
    download.file(META_URL, meta_path, mode = "wb")
  }
  if (!file.exists(count_path)) {
    message("Downloading count matrix (158 MB) — this takes a few minutes ...")
    download.file(COUNT_URL, count_path, mode = "wb")
  }
  list(meta_path = meta_path, count_path = count_path)
}

#' Build a SingleR reference from the downloaded files
#'
#' @param meta_path  Path to metadata txt.gz
#' @param count_path Path to count matrix txt.gz
#' @param aggregate  If TRUE, aggregate to per-celltype pseudobulk profiles
#'                   (much faster SingleR; recommended for >50k-cell references)
#' @param save_rds   Optional path to save the reference as .rds
#' @return SingleR reference (SingleCellExperiment or aggregated matrix)
build_reference <- function(meta_path, count_path, aggregate = TRUE,
                            save_rds = file.path(dirname(meta_path),
                                                 "gse149614_singler_ref.rds")) {

  message("Reading metadata ...")
  meta <- read.delim(gzfile(meta_path), stringsAsFactors = FALSE)
  stopifnot(all(c("Cell", "celltype") %in% colnames(meta)))

  message("Reading count matrix (genes x cells) ...")
  # fread (multi-threaded) is far faster than read.delim for the 158 MB matrix;
  # decompress via a gzip pipe so it works without R.utils.
  # The file uses an 'R-style' ragged header: line 1 lists only the cell
  # barcodes (no name for the gene-id column), so it has one fewer field than
  # the data rows. fread would discard those names on the mismatch, so read the
  # header separately and assign names manually.
  barcodes <- scan(text = readLines(gzfile(count_path), n = 1L),
                   what = character(), sep = "\t", quiet = TRUE)
  dt <- fread(cmd = paste("gzip -cd", shQuote(count_path)),
              sep = "\t", header = FALSE, skip = 1L, nThread = nthreads)
  gene_ids <- dt[[1]]
  dt[[1]] <- NULL
  counts <- as(as.matrix(dt), "dgCMatrix")
  rownames(counts) <- gene_ids
  # Drop a leading label field if the header also names the gene column.
  if (length(barcodes) == ncol(counts) + 1L) barcodes <- barcodes[-1L]
  stopifnot(length(barcodes) == ncol(counts))
  colnames(counts) <- barcodes
  rm(dt); gc()

  # Align cells between counts and metadata
  common <- intersect(colnames(counts), meta$Cell)
  message(sprintf("Cells: %d in counts, %d in metadata, %d matched",
                  ncol(counts), nrow(meta), length(common)))
  counts <- counts[, common]
  labels <- meta$celltype[match(common, meta$Cell)]

  # Log-normalize (SingleR expects log-normalized expression)
  lib <- colSums(counts)
  norm <- log1p(t(t(counts) / lib * 1e4))

  sce <- SingleCellExperiment(assays = list(logcounts = norm))
  sce$label <- labels

  if (aggregate) {
    message("Aggregating to per-celltype pseudobulk reference ...")
    ref <- aggregateReference(sce, labels = sce$label,
                              assay.type = "logcounts",
                              BPPARAM = BiocParallel::MulticoreParam(nthreads))
  } else {
    ref <- sce
  }

  if (!is.null(save_rds)) {
    saveRDS(ref, save_rds)
    message("Saved reference: ", save_rds)
  }
  ref
}

# --- Run as script -----------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "reference_data"
  paths <- fetch_reference(out_dir)
  invisible(build_reference(paths$meta_path, paths$count_path))
}
