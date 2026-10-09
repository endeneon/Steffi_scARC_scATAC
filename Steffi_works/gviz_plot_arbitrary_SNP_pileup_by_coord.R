#! /usr/bin/env Rscript

# extract macrophages ####
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
  library(biomaRt)

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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# setwd(
#   "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
# )
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

# load data ####

projMerged <-
  ArchR::loadArchRProject(
    path = "ArchR_merged_ATAC_multiome_obj_multiple_cell_types"
  )

projMacrophages <-
  ArchR::loadArchRProject(path = "ArchR_macrophages")

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

plot_snp_tracks <- function(
  snp_chr,
  snp_pos,
  gviz_window = 5000,
  gviz_bin_size = 50,
  ymax = NULL,
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

  # shared y-limit across all cell-type tracks: user-supplied or auto with headroom
  if (is.null(ymax)) {
    y_max <- max(mat_grp, na.rm = TRUE)
    if (!is.finite(y_max) || y_max <= 0) {
      y_max <- 1
    }
    ylim_use <- c(0, y_max * 1.05)
  } else {
    stopifnot(
      is.numeric(ymax),
      length(ymax) == 1,
      is.finite(ymax),
      ymax > 0
    )
    ylim_use <- c(0, ymax)
  }

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
  # projMacrophages <-
  #   ArchR::loadArchRProject(path = "ArchR_macrophages")
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
    main = main,
    cex.main = 1.1,
    col.main = "black"
  )
}

# rs10792832
plot_snp_tracks(
  snp_chr = "chr11",
  snp_pos = 86156833,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "PICALM, rs10792832 (chr11:86156833) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs150652488
plot_snp_tracks(
  snp_chr = "chr19",
  snp_pos = 45423658,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "ERCC1, rs2298881 (chr19:45423658) +/- 5 kb, 50 bp bins (RPGC)"
)

plot_snp_tracks(
  snp_chr = "chr19",
  snp_pos = 45479071,
  gviz_window = 1000,
  gviz_bin_size = 50,
  main = "ERCC1, rs150652488 (chr19:45479071) +/- 1 kb, 50 bp bins (RPGC)"
)

# rs2977306
plot_snp_tracks(
  snp_chr = "chr1",
  snp_pos = 17237472,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "PADI4, rs2977306 (chr1:17237472) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs3902981
plot_snp_tracks(
  snp_chr = "chr18",
  snp_pos = 12658192,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "PSMG2, rs3902981 (chr18:12658192) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs145809697
plot_snp_tracks(
  snp_chr = "chr6",
  snp_pos = 3751890,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "PXDC1, rs145809697 (chr6:3751890) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs1697138
plot_snp_tracks(
  snp_chr = "chr5",
  snp_pos = 67215959,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "CD180, rs1697138 (chr5:67215959) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs12648696
plot_snp_tracks(
  snp_chr = "chr4",
  snp_pos = 102625739,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "MANBA, rs12648696 (chr4:102625739) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs2027349
plot_snp_tracks(
  snp_chr = "chr1",
  snp_pos = 150067621,
  gviz_window = 5000,
  gviz_bin_size = 50,
  main = "VPS45, rs2027349 (chr1:150067621) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs10792832
plot_snp_tracks(
  snp_chr = "chr11",
  snp_pos = 86156833,
  gviz_window = 50000,
  gviz_bin_size = 50,
  main = "PICALM, rs10792832 (chr11:86156833) +/- 5 kb, 50 bp bins (RPGC)"
)

# rs3754884
plot_snp_tracks(
  snp_chr = "chr2",
  snp_pos = 98508913,
  gviz_window = 10000,
  gviz_bin_size = 50,
  ymax = 50,
  main = "INPP4A, rs3754884 (chr2:98508913) +/- 10 kb, 50 bp bins (RPGC)"
)
