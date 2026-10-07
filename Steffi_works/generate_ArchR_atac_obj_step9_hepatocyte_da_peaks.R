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
  # library(Gviz)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(org.Hs.eg.db)

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


# decide whether to add the GWAS track above the coverage tracks
annot_vcf_EAS <- "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/HCC_GWAS/GRCh38/hcc_ea_011123.vcf.gz"
annot_vcf_EUR <- "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/HCC_GWAS/GRCh38/hcc_eur_200324.vcf.gz"
use_GWAS_track <- TRUE

projMerged <-
  ArchR::loadArchRProject(
    path = "ArchR_merged_ATAC_multiome_obj_multiple_cell_types"
  )

df_raw <-
  read.table(
    "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    # annotation fields contain apostrophes; disable quote handling so
    # read.table() doesn't hit "EOF within quoted string" and drop rows
    quote = ""
  )

df_raw <-
  df_raw <- df_raw[!duplicated(df_raw$variantID), ]

df_2_plot <-
  df_raw[, c(
    "seqnames",
    "start",
    "variantID",
    "SYMBOL"
  )]


# count the Hepatocyte-specific DA intervals
# ---- Hepatocyte-specificity filter on +/- 250 bp SNP windows (RPGC) ---------
# window splits into exactly 10 x 50 bp bins). Fragments are assigned to bins
# by midpoint and scaled per group as in plot_snp_tracks():
#   RPGC = count * egs / (n_frags_grp * frag_len)
# Self-contained: helpers defined further below are not used here.
out_dir <- "DA_SNP_by_cell_type"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
target_group <- "Hepatocyte"
fold_cutoff <- 2.5
min_frags_window <- 10 # min raw fragment count in window, summed over groups
win_up <- 249
win_down <- 249
filt_bin_size <- 50

df_windows <- data.frame(
  variantID = df_2_plot$variantID,
  seqnames = as.character(df_2_plot$seqnames),
  snp_pos = as.integer(df_2_plot$start),
  win_start = as.integer(df_2_plot$start) - win_up,
  win_end = as.integer(df_2_plot$start) + win_down
)

filt_cells <- rownames(projMerged@cellColData)
filt_groups <- as.character(
  projMerged@cellColData$projected_barcodes_multiple_cell_types
)
filt_keep <- !is.na(filt_groups)
filt_cells <- filt_cells[filt_keep]
filt_groups <- filt_groups[filt_keep]
filt_levels <- sort(unique(filt_groups))
stopifnot(target_group %in% filt_levels)

# one row per 50 bp bin per window
n_bins_win <- (win_up + win_down + 1) %/% filt_bin_size
bins_df <- data.frame(
  win_idx = rep(seq_len(nrow(df_windows)), each = n_bins_win),
  seqnames = rep(df_windows$seqnames, each = n_bins_win),
  bin_start = rep(df_windows$win_start, each = n_bins_win) +
    rep((seq_len(n_bins_win) - 1L) * filt_bin_size, nrow(df_windows))
)
bins_gr <- GenomicRanges::GRanges(
  bins_df$seqnames,
  IRanges::IRanges(start = bins_df$bin_start, width = filt_bin_size)
)
filt_chrs <- unique(bins_df$seqnames)
filt_arrows <- ArchR::getArrowFiles(projMerged)

# Count fragments for one Arrow file: read once per chromosome, count all
# windows at once. Each worker opens its own Arrow file (HDF5 handles cannot be
# shared across processes) and returns a bins x groups count matrix.
count_arrow_bins <- function(af) {
  counts <- matrix(
    0,
    nrow = length(bins_gr),
    ncol = length(filt_levels),
    dimnames = list(NULL, filt_levels)
  )
  cells_af <- intersect(filt_cells, ArchR:::.availableCells(af, "TileMatrix"))
  if (!length(cells_af)) {
    return(counts)
  }
  for (chr in filt_chrs) {
    fr <- tryCatch(
      ArchR:::.getFragsFromArrow(
        af,
        chr = chr,
        out = "GRanges",
        cellNames = cells_af
      ),
      error = function(e) NULL
    )
    if (is.null(fr) || !length(fr)) {
      next
    }
    # fragment midpoint; floor() keeps x.5 midpoints in the same bin as
    # findInterval() does in plot_snp_tracks()
    mids <- floor(GenomicRanges::start(fr) + (GenomicRanges::width(fr) - 1) / 2)
    mid_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(mids, width = 1))
    hits <- GenomicRanges::findOverlaps(mid_gr, bins_gr)
    if (!length(hits)) {
      next
    }
    frag_grp <- filt_groups[match(
      as.character(S4Vectors::mcols(fr)$RG)[S4Vectors::queryHits(hits)],
      filt_cells
    )]
    ok <- !is.na(frag_grp)
    tab <- table(
      factor(S4Vectors::subjectHits(hits)[ok], levels = seq_along(bins_gr)),
      factor(frag_grp[ok], levels = filt_levels)
    )
    counts <- counts + unclass(tab)
  }
  counts
}

use_par_filt <- foreach::getDoParRegistered() &&
  foreach::getDoParWorkers() > 1 &&
  length(filt_arrows) > 1
if (use_par_filt) {
  `%dopar%` <- foreach::`%dopar%`
  bin_counts_list <- foreach::foreach(
    af = filt_arrows,
    .packages = c("ArchR", "GenomicRanges", "S4Vectors", "IRanges"),
    # export only what the worker needs; keeps projMerged out of the globals
    .export = c(
      "count_arrow_bins",
      "bins_gr",
      "filt_cells",
      "filt_groups",
      "filt_levels",
      "filt_chrs"
    ),
    .noexport = "projMerged",
    .errorhandling = "stop"
  ) %dopar%
    suppressPackageStartupMessages(suppressMessages(count_arrow_bins(af)))
} else {
  bin_counts_list <- lapply(filt_arrows, count_arrow_bins)
}
bin_counts <- Reduce(`+`, bin_counts_list)
rm(bin_counts_list)

# RPGC scale factors per group
filt_ga <- ArchR::getGenomeAnnotation(projMerged)
filt_egs <- sum(as.numeric(GenomicRanges::width(filt_ga$chromSizes)))
if (!is.null(filt_ga$blacklist) && length(filt_ga$blacklist)) {
  filt_egs <- filt_egs -
    sum(as.numeric(GenomicRanges::width(
      GenomicRanges::reduce(filt_ga$blacklist)
    )))
}

# mean fragment length estimated once (most frequent SNP chromosome), reused
set.seed(5813)
frag_len_chr <- names(sort(table(df_windows$seqnames), decreasing = TRUE))[1]
frag_sample <- unlist(lapply(filt_levels, function(lv) {
  cl <- filt_cells[filt_groups == lv]
  if (length(cl) > 300) sample(cl, 300) else cl
}))
# per-Arrow width sum and count (not raw widths) to keep worker returns small;
# sum/n pooled across Arrow files gives the same mean as pooling the widths
width_stats_arrow <- function(af) {
  cl <- intersect(frag_sample, ArchR:::.availableCells(af, "TileMatrix"))
  if (!length(cl)) {
    return(c(sum = 0, n = 0))
  }
  fr <- tryCatch(
    ArchR:::.getFragsFromArrow(
      af,
      chr = frag_len_chr,
      out = "GRanges",
      cellNames = cl
    ),
    error = function(e) NULL
  )
  if (is.null(fr) || !length(fr)) {
    return(c(sum = 0, n = 0))
  }
  w <- as.numeric(GenomicRanges::width(fr))
  c(sum = sum(w), n = length(w))
}

if (use_par_filt) {
  width_stats <- foreach::foreach(
    af = filt_arrows,
    .packages = c("ArchR", "GenomicRanges"),
    .export = c("width_stats_arrow", "frag_sample", "frag_len_chr"),
    .noexport = "projMerged",
    .errorhandling = "stop"
  ) %dopar%
    suppressPackageStartupMessages(suppressMessages(width_stats_arrow(af)))
} else {
  width_stats <- lapply(filt_arrows, width_stats_arrow)
}
width_stats <- Reduce(`+`, width_stats)
filt_frag_len <- if (width_stats[["n"]] > 0) {
  width_stats[["sum"]] / width_stats[["n"]]
} else {
  warning("Could not estimate fragment length; falling back to 100 bp.")
  100
}

filt_n_frags <- vapply(
  filt_levels,
  function(lv) {
    sum(as.numeric(
      projMerged@cellColData[filt_cells[filt_groups == lv], "nFrags"]
    ))
  },
  numeric(1)
)
filt_scale <- filt_egs / (filt_n_frags * filt_frag_len)
bin_rpgc <- sweep(bin_counts, 2, filt_scale[colnames(bin_counts)], `*`)

# per-window sums over the 10 bins
win_rpgc <- rowsum(bin_rpgc, bins_df$win_idx, reorder = TRUE)
win_frags <- rowSums(rowsum(bin_counts, bins_df$win_idx, reorder = TRUE))

other_levels <- setdiff(filt_levels, target_group)
df_windows$total_frags <- win_frags
df_windows$target_rpgc <- win_rpgc[, target_group]
df_windows$others_mean_rpgc <- rowMeans(win_rpgc[, other_levels, drop = FALSE])
df_windows$fold_vs_others <- df_windows$target_rpgc /
  df_windows$others_mean_rpgc
df_windows$keep <- df_windows$target_rpgc > 0 &
  df_windows$total_frags >= min_frags_window &
  df_windows$target_rpgc >= fold_cutoff * df_windows$others_mean_rpgc

colnames(win_rpgc) <- paste0("rpgc_", colnames(win_rpgc))
write.table(
  cbind(df_windows, win_rpgc),
  file.path(out_dir, "SNP_window_RPGC_hepatocyte_filter.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(paste0(
  "Hepatocyte filter kept ",
  sum(df_windows$keep),
  " of ",
  nrow(df_windows),
  " SNPs."
))
# carry fold_vs_others (same values as the TSV) for the panel titles
df_2_plot$fold_vs_others <- df_windows$fold_vs_others
df_2_plot <- df_2_plot[df_windows$keep, , drop = FALSE]


df_2_plot$annot <-
  stringr::str_c(
    df_2_plot$SYMBOL,
    ", ",
    df_2_plot$variantID,
    " (",
    df_2_plot$seqnames,
    ":",
    df_2_plot$start,
    ") +/- 5 kb, 50 bp bins (RPGC)",
    sep = ""
  )


# ---- helpers (inlined from plot_gviz_SPI1_macrophage.R) ----------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

# like %||% but also catches NA lookups from named-vector indexing
`%NA%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

.cytoband_cache <- new.env(parent = emptyenv())

# Local cytoband table for IdeogramTrack. Gviz otherwise queries UCSC at plot
# time, which fails on compute nodes without outbound network access. Refresh
# with:
#   curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz
#   mv cytoBandIdeo.txt.gz cytoBandIdeo_hg38.txt.gz
load_cytobands <-
  function(
    genome = "hg38",
    file = paste0("cytoBandIdeo_", genome, ".txt.gz")
  ) {
    key <- paste(genome, file, sep = "_")
    if (!is.null(.cytoband_cache[[key]])) {
      return(.cytoband_cache[[key]])
    }
    if (!file.exists(file)) {
      return(NULL)
    }
    bands <- utils::read.table(
      file,
      sep = "\t",
      header = FALSE,
      col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
      stringsAsFactors = FALSE
    )
    .cytoband_cache[[key]] <- bands
    bands
  }

# Effective genome size for RPGC: ArchR chromSizes minus blacklisted bp.
.effective_genome_size <-
  function(ref_ArchR_obj) {
    ga <- ArchR::getGenomeAnnotation(ref_ArchR_obj)
    total <- sum(as.numeric(GenomicRanges::width(ga$chromSizes)))
    bl <- ga$blacklist
    if (!is.null(bl) && length(bl)) {
      total <- total -
        sum(as.numeric(GenomicRanges::width(GenomicRanges::reduce(bl))))
    }
    total
  }

# Iterate over Arrow files in parallel when a foreach backend is registered,
# falling back to a serial lapply otherwise. Arrow/HDF5 handles cannot be
# shared across processes, so each worker opens its own file.
.arrow_lapply <-
  function(arrow_files, FUN, ...) {
    use_par <- requireNamespace("foreach", quietly = TRUE) &&
      foreach::getDoParRegistered() &&
      foreach::getDoParWorkers() > 1 &&
      length(arrow_files) > 1
    if (!use_par) {
      return(lapply(arrow_files, FUN, ...))
    }
    af <- NULL # silence R CMD check note for the foreach variable
    `%dopar%` <- foreach::`%dopar%`
    # each worker attaches ArchR, which prints a banner; keep it quiet
    wrapped <- function(x, ...) {
      suppressPackageStartupMessages(suppressMessages(FUN(x, ...)))
    }
    foreach::foreach(
      af = arrow_files,
      .packages = c("ArchR", "Matrix"),
      .errorhandling = "stop"
    ) %dopar%
      wrapped(af, ...)
  }

# Mean fragment width, estimated from a random sample of cells per group.
# Reading every fragment genome-wide is not tractable, but the mean stabilises
# quickly, so a few hundred cells per group is enough for a scale factor.
.estimate_frag_length <-
  function(
    arrow_files,
    cells_use,
    groups,
    grp_levels,
    chr,
    n_sample = 300L,
    seed = 5813L
  ) {
    set.seed(seed)
    sampled <- unlist(lapply(grp_levels, function(lv) {
      cl <- cells_use[groups == lv]
      if (length(cl) > n_sample) sample(cl, n_sample) else cl
    }))
    widths <- .arrow_lapply(arrow_files, function(af) {
      cl <- intersect(sampled, ArchR:::.availableCells(af, "TileMatrix"))
      if (!length(cl)) {
        return(numeric(0))
      }
      fr <- tryCatch(
        ArchR:::.getFragsFromArrow(
          af,
          chr = chr,
          out = "GRanges",
          cellNames = cl
        ),
        error = function(e) NULL
      )
      if (is.null(fr) || !length(fr)) {
        return(numeric(0))
      }
      as.numeric(GenomicRanges::width(fr))
    })
    widths <- unlist(widths)
    if (!length(widths)) {
      warning(
        "Could not estimate fragment length from the Arrow files; ",
        "falling back to 100 bp."
      )
      return(100)
    }
    mean(widths)
  }

# Pull GWAS records in [from, to] from a tabix-indexed VCF and return
# CHROM/POS/-log10(PVAL). scanTabix reads only the region (the full VCFs are
# hundreds of MB); vcfR parses the returned lines and extracts INFO/PVAL.
.read_gwas_region <- function(vcf_file, chr, from, to, ancestry) {
  param <- GenomicRanges::GRanges(chr, IRanges::IRanges(from, to))
  lines <- tryCatch(
    Rsamtools::scanTabix(vcf_file, param = param)[[1]],
    error = function(e) character(0)
  )
  if (!length(lines)) {
    return(NULL)
  }
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(
    c(Rsamtools::headerTabix(vcf_file)$header, lines),
    tmp
  )
  vcf <- vcfR::read.vcfR(tmp, verbose = FALSE)
  pval <- as.numeric(vcfR::extract.info(vcf, element = "PVAL"))
  out <- data.frame(
    CHROM = vcfR::getCHROM(vcf),
    POS = vcfR::getPOS(vcf),
    mlog10p = -log10(pval),
    ancestry = ancestry
  )
  out[is.finite(out$mlog10p), , drop = FALSE]
}

plot_snp_tracks <- function(
  snp_chr,
  snp_pos,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = sprintf(
    "%s:%s +/- %s bp, %s bp bins (RPGC)",
    snp_chr,
    format(snp_pos, big.mark = ","),
    format(gviz_window, big.mark = ","),
    gviz_bin_size
  )
) {
  from <- snp_pos - gviz_window
  to <- snp_pos + gviz_window

  # ---- cells and groups -------------------------------------------------------
  cells_use <- rownames(projMerged@cellColData)
  groups <- as.character(
    projMerged@cellColData$projected_barcodes_multiple_cell_types
  )
  keep <- !is.na(groups)
  cells_use <- cells_use[keep]
  groups <- groups[keep]
  grp_levels <- sort(unique(groups))

  # ---- 50 bp region-local pileups from fragments ------------------------------
  bin_starts <- seq(
    floor(from / gviz_bin_size) * gviz_bin_size,
    to,
    by = gviz_bin_size
  )
  region_gr <- GenomicRanges::GRanges(snp_chr, IRanges::IRanges(from, to))

  arrow_files <- ArchR::getArrowFiles(projMerged)
  mat_grp <- matrix(
    0,
    nrow = length(bin_starts),
    ncol = length(grp_levels),
    dimnames = list(NULL, grp_levels)
  )
  for (af in arrow_files) {
    cells_af <- intersect(cells_use, ArchR:::.availableCells(af, "TileMatrix"))
    if (!length(cells_af)) {
      next
    }
    fr <- tryCatch(
      ArchR:::.getFragsFromArrow(
        af,
        chr = snp_chr,
        out = "GRanges",
        cellNames = cells_af
      ),
      error = function(e) NULL
    )
    if (is.null(fr) || !length(fr)) {
      next
    }
    fr <- IRanges::subsetByOverlaps(fr, region_gr)
    if (!length(fr)) {
      next
    }
    # count each fragment in the bin holding its midpoint
    mids <- GenomicRanges::start(fr) + (GenomicRanges::width(fr) - 1) / 2
    bin_idx <- findInterval(mids, bin_starts)
    cell_ids <- as.character(S4Vectors::mcols(fr)$RG)
    frag_grp <- groups[match(cell_ids, cells_use)]
    ok <- bin_idx >= 1 &
      bin_idx <= length(bin_starts) &
      !is.na(frag_grp)
    if (!any(ok)) {
      next
    }
    tab <- table(
      factor(bin_idx[ok], levels = seq_along(bin_starts)),
      factor(frag_grp[ok], levels = grp_levels)
    )
    mat_grp <- mat_grp + unclass(tab)
  }

  # ---- RPGC normalization per cell type ---------------------------------------
  egs <- .effective_genome_size(projMerged)
  frag_len <- .estimate_frag_length(
    arrow_files,
    cells_use,
    groups,
    grp_levels,
    chr = snp_chr
  )
  n_frags_grp <- vapply(
    grp_levels,
    function(lv) {
      sum(as.numeric(projMerged@cellColData[cells_use[groups == lv], "nFrags"]))
    },
    numeric(1)
  )
  scale_grp <- egs / (n_frags_grp * frag_len)
  for (lv in grp_levels) {
    mat_grp[, lv] <- mat_grp[, lv] * scale_grp[[lv]]
  }

  tile_gr <- GenomicRanges::GRanges(
    seqnames = snp_chr,
    ranges = IRanges::IRanges(start = bin_starts, width = gviz_bin_size)
  )

  # ---- tracks -----------------------------------------------------------------
  bands <- load_cytobands()
  ideo_track <- tryCatch(
    if (is.null(bands)) {
      Gviz::IdeogramTrack(genome = "hg38", chromosome = snp_chr)
    } else {
      Gviz::IdeogramTrack(
        genome = "hg38",
        chromosome = snp_chr,
        bands = bands[bands$chrom == snp_chr, , drop = FALSE]
      )
    },
    error = function(e) NULL
  )
  axis_track <- Gviz::GenomeAxisTrack(
    range = region_gr,
    add53 = TRUE,
    add35 = TRUE
  )

  pal <- setNames(
    rep_len(
      RColorBrewer::brewer.pal(max(3L, min(8L, length(grp_levels))), "Dark2"),
      length(grp_levels)
    ),
    grp_levels
  )

  # shared y-limit across all cell-type tracks, with headroom
  y_max <- max(mat_grp, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }
  ylim_use <- c(0, y_max * 1.05)

  # one vertical track per cell type
  cov_tracks <- lapply(grp_levels, function(lv) {
    Gviz::DataTrack(
      range = tile_gr,
      data = matrix(mat_grp[, lv], nrow = 1),
      genome = "hg38",
      chromosome = snp_chr,
      name = lv,
      type = "polygon",
      fill.mountain = rep(pal[[lv]], 2),
      col.mountain = NA,
      ylim = ylim_use,
      showAxis = TRUE,
      cex.axis = 0.6
    )
  })

  # peak annotation from the ArchR_hepato project's peak set
  # (ArchR_hepatocytes_subset has no Save-ArchR-Project.rds and cannot be loaded)
  projPeakAnnot <-
    ArchR::loadArchRProject(path = "ArchR_hepato")
  peak_track <- NULL
  ps <- tryCatch(ArchR::getPeakSet(projPeakAnnot), error = function(e) NULL)
  if (!is.null(ps) && length(ps)) {
    ps <- IRanges::subsetByOverlaps(ps, region_gr)
  }
  if (!is.null(ps) && length(ps)) {
    peak_col <- c(
      Promoter = "#E7298A",
      Distal = "#7570B3",
      Exonic = "#66A61E",
      Intronic = "#E6AB02"
    )
    ptype <- as.character(S4Vectors::mcols(ps)$peakType)
    if (is.null(ptype)) {
      ptype <- rep("Peak", length(ps))
    }
    ptype[is.na(ptype)] <- "Other"
    type_cols <- vapply(
      unique(ptype),
      function(tt) unname(peak_col[tt]) %NA% "grey50",
      character(1)
    )
    peak_track <- Gviz::AnnotationTrack(
      range = GenomicRanges::GRanges(snp_chr, IRanges::ranges(ps)),
      genome = "hg38",
      chromosome = snp_chr,
      name = "Peaks",
      feature = ptype,
      stacking = "dense",
      col = NA,
      showFeatureId = FALSE
    )
    Gviz::displayPars(peak_track) <- as.list(type_cols)
  }

  # squished gene track at the bottom
  gene_track <- Gviz::GeneRegionTrack(
    TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    chromosome = snp_chr,
    start = from,
    end = to,
    genome = "hg38",
    name = "Genes",
    stacking = "squish",
    collapseTranscripts = "longest"
  )
  if (length(Gviz::gene(gene_track))) {
    sym <- suppressMessages(tryCatch(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = sub("\\..*$", "", Gviz::gene(gene_track)),
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      ),
      error = function(e) NULL
    ))
    if (!is.null(sym)) {
      sym[is.na(sym)] <- Gviz::gene(gene_track)[is.na(sym)]
      Gviz::symbol(gene_track) <- unname(sym)
    }
  }

  # optional GWAS -log10(PVAL) dot track (EAS + EUR), placed at the top
  gwas_track <- NULL
  if (isTRUE(use_GWAS_track)) {
    gwas_df <- rbind(
      .read_gwas_region(annot_vcf_EAS, snp_chr, from, to, "EAS"),
      .read_gwas_region(annot_vcf_EUR, snp_chr, from, to, "EUR")
    )
    if (!is.null(gwas_df) && nrow(gwas_df)) {
      pos_u <- sort(unique(gwas_df$POS))
      # one row per ancestry; NA where a position is absent from that VCF
      gwas_mat <- vapply(
        c("EAS", "EUR"),
        function(anc) {
          sub_df <- gwas_df[gwas_df$ancestry == anc, ]
          sub_df$mlog10p[match(pos_u, sub_df$POS)]
        },
        numeric(length(pos_u))
      )
      # fix the y-axis to 0-5 unless a point exceeds 5, then autoscale up to it
      gwas_ymax <- suppressWarnings(max(gwas_mat, na.rm = TRUE))
      if (!is.finite(gwas_ymax) || gwas_ymax <= 5) {
        ylim_gwas <- c(0, 5)
      } else {
        ylim_gwas <- c(0, gwas_ymax)
      }
      gwas_points <- Gviz::DataTrack(
        range = GenomicRanges::GRanges(
          snp_chr,
          IRanges::IRanges(start = pos_u, width = 1)
        ),
        data = t(gwas_mat),
        groups = factor(c("EAS", "EUR"), levels = c("EAS", "EUR")),
        genome = "hg38",
        chromosome = snp_chr,
        name = "HCC GWAS\n-log10(P)",
        type = "p",
        col = c("darkred", "darkblue"),
        alpha = 0.8,
        legend = FALSE,
        cex = 0.6,
        cex.axis = 0.6,
        ylim = ylim_gwas,
        # y = 0 lower-margin line (single baseline; vector baselines error in
        # Gviz, so the significance line is a separate overlay track below)
        baseline = 0,
        col.baseline = "grey70",
        lty.baseline = "solid",
        lwd.baseline = 0.8
      )
      # genome-wide significance line (P = 5e-8, -log10 = 7.3) as its own flat
      # DataTrack, overlaid on the points so it shares the same y-axis.
      sig_line <- Gviz::DataTrack(
        range = GenomicRanges::GRanges(
          snp_chr,
          IRanges::IRanges(start = c(from, to), width = 1)
        ),
        data = matrix(rep(-log10(5e-8), 2), nrow = 1),
        genome = "hg38",
        chromosome = snp_chr,
        type = "l",
        col = "grey50",
        lty = "dashed",
        lwd = 0.8,
        ylim = ylim_gwas
      )
      gwas_track <- Gviz::OverlayTrack(
        trackList = list(gwas_points, sig_line),
        name = "HCC GWAS\n-log10(P)"
      )
    }
  }

  # highlight the SNP position across the data tracks
  data_tracks <- c(
    if (!is.null(gwas_track)) list(gwas_track),
    cov_tracks,
    if (!is.null(peak_track)) list(peak_track),
    list(gene_track)
  )
  hl <- Gviz::HighlightTrack(
    trackList = data_tracks,
    start = snp_pos,
    end = snp_pos,
    chromosome = snp_chr,
    col = "red",
    fill = "#FFE9E9",
    inBackground = TRUE
  )

  track_sizes <- c(
    0.6, # ideogram
    0.8, # axis
    if (!is.null(gwas_track)) 3, # GWAS dots
    rep(3, length(cov_tracks)), # one coverage track per cell type
    if (!is.null(peak_track)) 1, # peaks
    2 # genes (squished)
  )

  Gviz::plotTracks(
    c(
      Filter(Negate(is.null), list(ideo_track, axis_track)),
      list(hl)
    ),
    from = from,
    to = to,
    chromosome = snp_chr,
    transcriptAnnotation = "symbol",
    shape = "arrow",
    sizes = track_sizes[
      seq_len(length(track_sizes) - is.null(ideo_track))
    ],
    background.title = "white",
    col.title = "black",
    fontcolor.title = "black",
    col.axis = "black",
    col.border.title = "transparent",
    cex.title = 0.7,
    main = main,
    cex.main = 1.1,
    col.main = "black"
  )
}


# # rs150652488
# plot_snp_tracks(
#   snp_chr = "chr19",
#   snp_pos = 45479071,
#   gviz_window = 5000,
#   gviz_bin_size = 50,
#   main = "ERCC1, rs150652488 (chr19:45479071) +/- 5 kb, 50 bp bins (RPGC)"
# )

# plot_snp_tracks(
#   snp_chr = "chr19",
#   snp_pos = 45479071,
#   gviz_window = 1000,
#   gviz_bin_size = 50,
#   main = "ERCC1, rs150652488 (chr19:45479071) +/- 1 kb, 50 bp bins (RPGC)"
# )

# ---- batch SNP track panels -> 3x3 multipage PDF ----------------------------

# widen the window when the GWAS track is shown so more variants are visible
gviz_window <- if (isTRUE(use_GWAS_track)) 50000 else 5000
gviz_bin_size <- 50

# Capture each row's plotTracks() output as a grid grob (plotTracks draws to a
# device, so grid.grabExpr() records it without opening one). Built in parallel
# over the foreach/doFuture backend registered above; each worker re-derives its
# own Arrow handles inside plot_snp_tracks() (HDF5 handles can't cross processes).
if (!foreach::getDoParRegistered()) {
  doFuture::registerDoFuture()
  future::plan(session_plan, workers = workers_2_use)
}

i <- NULL # silence R CMD check note for the foreach variable
`%dopar%` <- foreach::`%dopar%`
gviz_grobs <- foreach::foreach(
  i = seq_len(nrow(df_2_plot)),
  .packages = c(
    "Gviz",
    "grid",
    "gridExtra",
    "ArchR",
    "GenomicRanges",
    "Rsamtools",
    "vcfR"
  ),
  .errorhandling = "pass"
) %dopar%
  {
    cat("Processing row ", i, "\n")
    row <- df_2_plot[i, ]
    main_i <- sprintf(
      "%s, %s (%s:%s) +/- %s bp, %s bp bins (RPGC), \nHepatocyte fold vs others = %s",
      row$SYMBOL,
      row$variantID,
      row$seqnames,
      format(row$start, big.mark = ","),
      format(gviz_window, big.mark = ","),
      gviz_bin_size,
      # Inf when all other groups are 0
      formatC(row$fold_vs_others, format = "f", digits = 2)
    )
    panel <- grid::grid.grabExpr(
      plot_snp_tracks(
        snp_chr = row$seqnames,
        snp_pos = row$start,
        gviz_window = gviz_window,
        gviz_bin_size = gviz_bin_size,
        main = main_i
      ),
      # Gviz reuses grob names across its per-track panels; without wrapping,
      # grid.grabExpr() silently overwrites all but one coverage track (the grab
      # is "not faithful"), collapsing the 6 cell-type tracks to a single panel.
      wrap.grobs = TRUE
    )
    # attach a short label (SYMBOL, variantID, chr:pos) as a per-panel subtitle
    subtitle_i <- sprintf(
      "%s, %s (%s:%s)",
      row$SYMBOL,
      row$variantID,
      row$seqnames,
      format(row$start, big.mark = ",")
    )
    gridExtra::arrangeGrob(
      panel,
      bottom = grid::textGrob(subtitle_i, gp = grid::gpar(fontsize = 5))
    )
  }

# With .errorhandling = "pass", a failed SNP returns its error condition instead
# of a grob; drop those so one bad SNP doesn't abort the batch, and report which.
failed <- !vapply(gviz_grobs, grid::is.grob, logical(1))
if (any(failed)) {
  warning(
    sum(failed),
    " of ",
    length(gviz_grobs),
    " panels failed and were dropped: rows ",
    paste(which(failed), collapse = ", ")
  )
  for (j in which(failed)) {
    message(
      "  row ",
      j,
      " (",
      df_2_plot$variantID[j],
      "): ",
      conditionMessage(gviz_grobs[[j]])
    )
  }
  gviz_grobs <- gviz_grobs[!failed]
}

out_dir <- "DA_SNP_by_cell_type"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2x2 grid per landscape page; larger page keeps the 6 stacked tracks legible.
pdf(
  file.path(out_dir, "DA_SNP_by_cell_type_panels_hepatocyte.pdf"),
  width = 14,
  height = 10.5
)
per_page <- 4
n_pages <- ceiling(length(gviz_grobs) / per_page)
for (pg in seq_len(n_pages)) {
  idx <- seq((pg - 1) * per_page + 1, min(pg * per_page, length(gviz_grobs)))
  gridExtra::grid.arrange(
    grobs = gviz_grobs[idx],
    nrow = 2,
    ncol = 2,
    bottom = grid::textGrob(
      sprintf("Page %d of %d", pg, n_pages),
      x = 0.98,
      hjust = 1,
      gp = grid::gpar(fontsize = 8, col = "grey40")
    )
  )
}
invisible(dev.off())

# # rs150652488
# plot_snp_tracks(
#   snp_chr = "chr19",
#   snp_pos = 45479071,
#   gviz_window = 5000,
#   gviz_bin_size = 50,
#   main = "ERCC1, rs150652488 (chr19:45479071) +/- 5 kb, 50 bp bins (RPGC)"
# )

# plot_snp_tracks(
#   snp_chr = "chr19",
#   snp_pos = 45479071,
#   gviz_window = 1000,
#   gviz_bin_size = 50,
#   main = "ERCC1, rs150652488 (chr19:45479071) +/- 1 kb, 50 bp bins (RPGC)"
# )
