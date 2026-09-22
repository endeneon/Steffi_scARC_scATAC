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

selected_cells_obj <-
  qs_read(
    "annotation_script_package/reference_based/annotated_singler_reference_based.qs2",
    nthreads = 8
  )
selected_cells_obj$simplified_singler_label <-
  str_split(
    selected_cells_obj$singler_label_cell,
    pattern = "\\.",
    simplify = TRUE
  )[, 1]
qs_save(
  selected_cells_obj,
  "annotation_script_package/reference_based/annotated_singler_reference_based.qs2",
  nthreads = 8
)

# use cluster-based
selected_cells_obj <-
  qs_read(
    "annotation_script_package/marker_based/step5.qs2",
    nthreads = 8
  )

df_barcodes_w_ident <-
  data.frame(
    barcodes_w_ident = colnames(selected_cells_obj),
    celltype_broad = selected_cells_obj$celltype_broad
  )
df_barcodes_w_ident$barcodes_w_ident_transfer <-
  sub(
    "^([^_]*__[^_]*)_",
    "\\1#",
    df_barcodes_w_ident$barcodes_w_ident
  )

sum(
  rownames(projMerged@cellColData) %in%
    df_barcodes_w_ident$barcodes_w_ident_transfer
)

# add transferred_barcode column to cellColData and fill the rest with NA
projMerged@cellColData$transferred_celltype_multiple <- NA
idx <- match(
  rownames(projMerged@cellColData),
  df_barcodes_w_ident$barcodes_w_ident_transfer
)
projMerged@cellColData$transferred_celltype_multiple <- as.character(
  df_barcodes_w_ident$celltype_broad
)[idx]
# sum(!is.na(idx))

# transfer barcodes with PeakMatrix UMAP embedding
# plotEmbedding(
#   ArchRProj = projMerged,
#   embedding = "UMAP_Peaks",
#   colorBy = "cellColData",
#   name = "transferred_barcode",
#   size = 0.1,
#   # sampleCells = 10000,
#   highlightCells = NULL,
#   rastr = FALSE,
#   quantCut = c(0.01, 0.99),
#   discreteSet = NULL,
#   continuousSet = NULL,
#   randomize = TRUE,
#   keepAxis = FALSE,
#   baseSize = 10
# )

# name the extended barcode idents to "projected_barcodes_multiple_cell_types"
projMerged@cellColData$projected_barcodes_multiple_cell_types <- NA

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
umap_df$transferred_celltype_multiple <- cell_meta[
  rownames(umap_df),
  "transferred_celltype_multiple"
]

unique_cell_types <-
  unique(umap_df$transferred_celltype_multiple)[
    !is.na(unique(umap_df$transferred_celltype_multiple))
  ]

# --- Project each cell type's density onto all cells ----------------------
# For every cell type, rebuild the same 2D kernel density that
# stat_density_2d draws (MASS::kde2d with per-axis bandwidth.nrd, matching
# ggplot's defaults) on that type's labeled cells, over a grid spanning ALL
# cells, and evaluate the fitted density at every cell's UMAP_Peaks position.

# nearest grid cell indices are the same for every cell type (shared grid)
# skip cell types with too few cells or degenerate bandwidths (kde2d errors)
min_cells_for_kde <- 25
density_mat <- sapply(unique_cell_types, function(ct) {
  ct_df <- umap_df[
    !is.na(umap_df$transferred_celltype_multiple) &
      umap_df$transferred_celltype_multiple == ct,
  ]
  bw <- c(
    MASS::bandwidth.nrd(ct_df$UMAP1),
    MASS::bandwidth.nrd(ct_df$UMAP2)
  )
  if (nrow(ct_df) < min_cells_for_kde || any(!is.finite(bw)) || any(bw <= 0)) {
    warning(
      "Skipping cell type '",
      ct,
      "': too few cells (",
      nrow(ct_df),
      ") or degenerate bandwidth for kde2d"
    )
    return(rep(NA_real_, nrow(umap_df)))
  }
  kde_fit <- MASS::kde2d(
    x = ct_df$UMAP1,
    y = ct_df$UMAP2,
    h = bw,
    n = 200,
    lims = c(range(umap_df$UMAP1), range(umap_df$UMAP2))
  )
  ix <- findInterval(umap_df$UMAP1, kde_fit$x, all.inside = TRUE)
  iy <- findInterval(umap_df$UMAP2, kde_fit$y, all.inside = TRUE)
  kde_fit$z[cbind(ix, iy)]
})
rownames(density_mat) <- rownames(umap_df)

# drop skipped cell types so downstream winner/cutoff logic sees no NAs
kept <- colSums(is.na(density_mat)) < nrow(density_mat)
density_mat <- density_mat[, kept, drop = FALSE]
unique_cell_types <- unique_cell_types[kept]

# per-cell winner: cell type with the highest fitted density
winner_idx <- max.col(density_mat, ties.method = "first")
winner_type <- unique_cell_types[winner_idx]
winner_density <- density_mat[cbind(seq_len(nrow(density_mat)), winner_idx)]

# per-cell lower cut-off:
# (a) the 2nd-smallest of this cell's per-type densities
second_lowest <- apply(density_mat, 1, function(v) sort(v)[2])
# (b) the 5% quantile (across all cells) of the density of this cell's
#     lowest-density cell type
lowest_type_idx <- apply(density_mat, 1, which.min)
type_q05 <- apply(density_mat, 2, quantile, probs = 0.05, na.rm = TRUE)
q05_of_lowest_type <- type_q05[lowest_type_idx]
# take whichever is lower
cutoff <- pmin(second_lowest, q05_of_lowest_type)

# assign the winning cell type only where its density clears the cut-off
projected <- ifelse(winner_density >= cutoff, winner_type, NA)
projMerged@cellColData$projected_barcodes_multiple_cell_types <-
  projected[match(rownames(projMerged@cellColData), rownames(umap_df))]

table(
  projMerged@cellColData$projected_barcodes_multiple_cell_types,
  useNA = "ifany"
)

plotEmbedding(
  ArchRProj = projMerged,
  embedding = "UMAP_Peaks",
  colorBy = "cellColData",
  name = "projected_barcodes_multiple_cell_types",
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
) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

plotEmbedding(
  ArchRProj = projMerged,
  embedding = "UMAP_Tiles",
  colorBy = "cellColData",
  name = "projected_barcodes_multiple_cell_types",
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
) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

projMerged <-
  saveArchRProject(
    ArchRProj = projMerged,
    outputDirectory = "ArchR_merged_ATAC_multiome_obj_multiple_cell_types",
    load = TRUE,
    overwrite = TRUE
  )

# Gviz pileups around rs2298881, one track per projected cell type ####
# Adapted from plot_gviz_SPI1_macrophage.R, but instead of overlaying all
# pileups in one panel, each level of projected_barcodes_multiple_cell_types
# gets its own DataTrack, stacked vertically. Coverage is re-binned region-
# locally at 50 bp from the fragments (nothing written back to Arrow files),
# RPGC-normalized per cell type, and drawn on a shared y-axis.
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

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

snp_chr <- "chr19"
# snp_pos <- 45423658 # rs2298881
snp_pos <- 45479071 #rs150652488

gviz_window <- 5000
gviz_bin_size <- 50
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

# peak annotation from the ArchR_macrophages project's peak set
projMacrophages <-
  ArchR::loadArchRProject(path = "ArchR_macrophages")
peak_track <- NULL
ps <- tryCatch(ArchR::getPeakSet(projMacrophages), error = function(e) NULL)
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

# highlight the SNP position across the data tracks
data_tracks <- c(
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
  # main = "rs2298881 (chr19:45423658) +/- 5 kb, 50 bp bins (RPGC)",
  main = "rs150652488 (chr19:45479071) +/- 5 kb, 50 bp bins (RPGC)",
  cex.main = 1.1,
  col.main = "black"
)
