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

####
# projMerged <-
#   ArchR::loadArchRProject(path = "ArchR_merged_ATAC_multiome_obj")

# confident_macrophage_barcodes <-
#   qs_read(
#     "confident_macrophage_barcodes.qs2",
#     nthreads = 4
#   )

# projMacrophages <-
#   subsetArchRProject(
#     ArchRProj = projMerged,
#     cells = confident_macrophage_barcodes,
#     outputDirectory = "ArchR_macrophages",
#     dropCells = TRUE,
#     force = TRUE
#   )

# projMacrophages@cellColData$category <- "Primary"
# projMacrophages@cellColData$category[str_detect(
#   as.character(projMacrophages@cellColData$Sample),
#   pattern = "pre$",
#   negate = TRUE
# )] <- "Resistant"
# projMacrophages <-
#   saveArchRProject(
#     ArchRProj = projMacrophages,
#     outputDirectory = "ArchR_macrophages",
#     load = TRUE,
#     overwrite = TRUE
#   )

ref_ArchR_obj <-
  ArchR::loadArchRProject(path = "ArchR_macrophages")

df_sig_snp_list <-
  read.table(
    "SPI1_coord.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )

# plot_gviz_pileups_by_category() ####
# Build a Gviz track list (ideogram / axis / overlaid horizon coverage split by
# a cellColData column / gene models) for a slide window around a locus, using
# the TileMatrix stored in an ArchRProject's Arrow files.
#
# The tile counts for a given (region, slot, barcode-set) are cached in
# `.tile_cache` so that repeated calls on the same window do not re-read the
# Arrow files.

.tile_cache <- new.env(parent = emptyenv())
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

`%||%` <- function(x, y) if (is.null(x)) y else x

# like %||% but also catches NA lookups from named-vector indexing
`%NA%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# Draw a horizontal colour key. Gviz builds its plots with grid, so the legend
# is drawn into a grid viewport rather than with graphics::legend().
#' Draw the colour key(s). Each key is a block laid out side by side, with the
#' block title to the left of a grid of swatch/label pairs. Blocks with many
#' entries wrap into several columns so all blocks end up the same height.
#'
#' @param legend_rows Rows per block. Defaults to the smallest block's entry
#'   count (min 2), so e.g. a 4-entry key wraps to 2x2 beside a 2-entry key.
#' @return Invisibly, the number of rows used (for height calculations).
.draw_track_legend <-
  function(
    pal,
    alpha = 1,
    title = NULL,
    cex = 0.8,
    peak_pal = NULL,
    legend_rows = NULL,
    draw = TRUE
  ) {
    if (is.null(pal) || !length(pal)) {
      return(invisible(0L))
    }

    blocks <- list(
      list(
        title = title,
        labels = names(pal),
        fills = grDevices::adjustcolor(unname(pal), alpha.f = alpha)
      )
    )
    if (!is.null(peak_pal) && length(peak_pal)) {
      blocks[[length(blocks) + 1L]] <- list(
        title = "peakType",
        labels = names(peak_pal),
        fills = unname(unlist(peak_pal))
      )
    }

    sizes <- vapply(blocks, function(b) length(b$labels), integer(1))
    n_rows <- legend_rows %||% max(2L, min(sizes))
    n_rows <- max(1L, min(n_rows, max(sizes)))
    if (!draw) {
      return(invisible(n_rows))
    }

    # All geometry is resolved to inches up front. Measuring with stringWidth()
    # while drawing at a different cex is what previously let swatches overlap
    # their labels, so widths are measured at the cex actually used.
    txt_w <- function(s, bold = FALSE) {
      vapply(
        s,
        function(x) {
          grid::convertWidth(
            grid::stringWidth(x),
            "inches",
            valueOnly = TRUE
          ) *
            cex
        },
        numeric(1)
      )
    }
    line_in <- grid::convertHeight(
      grid::unit(1, "lines"),
      "inches",
      valueOnly = TRUE
    )
    vp_w <- grid::convertWidth(
      grid::unit(1, "npc"),
      "inches",
      valueOnly = TRUE
    )

    sw_in <- 0.7 * line_in * cex # swatch side
    pad_in <- 0.30 * line_in # swatch -> label
    colgap_in <- 0.75 * line_in # between wrapped columns
    titlegap_in <- 0.5 * line_in # title -> first column
    row_in <- 1.4 * line_in * cex

    for (bi in seq_along(blocks)) {
      b <- blocks[[bi]]
      k <- length(b$labels)
      n_col <- ceiling(k / n_rows)
      col_of <- ceiling(seq_len(k) / n_rows)
      row_of <- (seq_len(k) - 1L) %% n_rows + 1L

      lab_in <- txt_w(b$labels)
      col_w <- vapply(
        seq_len(n_col),
        function(cc) sw_in + pad_in + max(lab_in[col_of == cc]),
        numeric(1)
      )
      title_in <- if (is.null(b$title)) 0 else txt_w(b$title)
      total_in <- title_in +
        titlegap_in +
        sum(col_w) +
        colgap_in * (n_col - 1L)

      # left edge of this block, centred in its horizontal slice
      x0 <- (bi - 0.5) / length(blocks) * vp_w - total_in / 2

      if (!is.null(b$title)) {
        grid::grid.text(
          b$title,
          x = grid::unit(x0, "inches"),
          y = grid::unit(0.5, "npc"),
          just = c("left", "centre"),
          gp = grid::gpar(cex = cex, fontface = "bold", col = "black")
        )
      }

      rows_used <- min(k, n_rows)
      for (i in seq_len(k)) {
        cc <- col_of[i]
        x_sw <- x0 + title_in + titlegap_in
        if (cc > 1L) {
          x_sw <- x_sw + sum(col_w[seq_len(cc - 1L)]) + colgap_in * (cc - 1L)
        }
        yy <- grid::unit(0.5, "npc") +
          grid::unit((rows_used - 1) / 2 * row_in, "inches") -
          grid::unit((row_of[i] - 1) * row_in, "inches")

        grid::grid.rect(
          x = grid::unit(x_sw, "inches"),
          y = yy,
          width = grid::unit(sw_in, "inches"),
          height = grid::unit(sw_in, "inches"),
          just = c("left", "centre"),
          gp = grid::gpar(fill = b$fills[i], col = "grey30", lwd = 0.5)
        )
        grid::grid.text(
          b$labels[i],
          x = grid::unit(x_sw + sw_in + pad_in, "inches"),
          y = yy,
          just = c("left", "centre"),
          gp = grid::gpar(cex = cex, col = "black")
        )
      }
    }
    invisible(n_rows)
  }

# Shared by the main function and replot_gviz_tracks(): draw the tracks in the
# upper part of the device and the colour key in a strip underneath.
.plot_tracks_with_legend <-
  function(
    tracks,
    from,
    to,
    chr,
    pal,
    alpha,
    slot,
    legend = TRUE,
    legend_height = NULL,
    legend_rows = NULL,
    sizes = NULL,
    peak_pal = NULL,
    main_title = NULL,
    ...
  ) {
    show_legend <- isTRUE(legend) && !is.null(pal) && length(pal)

    # Size the strip from the wrapped legend so nothing clips: rows x line
    # height, plus padding, expressed as a fraction of the device height.
    if (show_legend && is.null(legend_height)) {
      n_rows <- .draw_track_legend(
        pal,
        alpha = alpha,
        title = slot,
        peak_pal = peak_pal,
        legend_rows = legend_rows,
        draw = FALSE
      )
      # n_rows of line-height plus generous padding; the extra allowance keeps
      # the bottom row clear of the device edge.
      strip_in <- grid::convertHeight(
        grid::unit(n_rows * 1.4 + 1.6, "lines"),
        "inches",
        valueOnly = TRUE
      )
      dev_in <- grid::convertHeight(
        grid::unit(1, "npc"),
        "inches",
        valueOnly = TRUE
      )
      legend_height <- min(0.4, strip_in / dev_in)
    }

    # Track title panels: black text and axis on a white strip. Gviz defaults
    # to white-on-grey, which is hard to read and does not print well.
    title_pars <- list(
      background.title = "white",
      col.title = "black",
      fontcolor.title = "black",
      col.axis = "black",
      col.border.title = "transparent",
      cex.title = 0.8
    )
    # plotTracks() draws no title when `main` is absent, so only pass it when
    # supplied rather than handing it a NULL.
    if (!is.null(main_title)) {
      title_pars$main <- main_title
      title_pars$cex.main <- 1.1
      title_pars$fontface.main <- "bold"
      title_pars$col.main <- "black"
    }

    if (!show_legend) {
      return(invisible(do.call(
        Gviz::plotTracks,
        c(
          list(
            tracks,
            from = from,
            to = to,
            chromosome = chr,
            transcriptAnnotation = "symbol",
            shape = "arrow",
            sizes = sizes
          ),
          title_pars,
          list(...)
        )
      )))
    }

    # Start a fresh page, but only if the current one has already been drawn
    # on. Calling grid.newpage() unconditionally emits a blank leading page on
    # a freshly opened device (e.g. pdf(); replot_gviz_tracks(); dev.off()).
    if (length(grid::grid.ls(print = FALSE)$name)) {
      grid::grid.newpage()
    } else {
      # clear any viewports left behind without advancing the page
      grid::upViewport(0L)
    }
    grid::pushViewport(grid::viewport(
      y = grid::unit(legend_height, "npc"),
      height = grid::unit(1 - legend_height, "npc"),
      just = "bottom"
    ))
    do.call(
      Gviz::plotTracks,
      c(
        list(
          tracks,
          from = from,
          to = to,
          chromosome = chr,
          transcriptAnnotation = "symbol",
          shape = "arrow",
          add = TRUE,
          sizes = sizes
        ),
        title_pars,
        list(...)
      )
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(
      y = grid::unit(0, "npc"),
      height = grid::unit(legend_height, "npc"),
      just = "bottom"
    ))
    .draw_track_legend(
      pal,
      alpha = alpha,
      title = slot,
      peak_pal = peak_pal,
      legend_rows = legend_rows
    )
    grid::popViewport()
    invisible(NULL)
  }

#' Re-draw a track list returned by plot_gviz_pileups_by_category(), reusing the
#' window, palette and alpha stored in its attributes.
replot_gviz_tracks <-
  function(
    tracks,
    legend = TRUE,
    legend_height = NULL,
    legend_rows = NULL,
    main_title = NULL,
    ...
  ) {
    .plot_tracks_with_legend(
      tracks,
      from = attr(tracks, "from"),
      to = attr(tracks, "to"),
      chr = attr(tracks, "chromosome"),
      pal = attr(tracks, "palette"),
      alpha = attr(tracks, "alpha") %||% 1,
      slot = attr(tracks, "slot"),
      legend = legend,
      legend_height = legend_height,
      legend_rows = legend_rows,
      sizes = attr(tracks, "sizes"),
      peak_pal = attr(tracks, "peak_palette"),
      main_title = main_title %||% attr(tracks, "main_title"),
      ...
    )
    invisible(tracks)
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

#'   `start - window` to `end + window`.
#' @param window Flanking bp added on each side.
#' @param ref_ArchR_obj An ArchRProject with a TileMatrix.
#' @param slot Name of a `cellColData` column used to group / colour cells.
#' @param barcodes Optional character vector of cell names to restrict to;
#'   `NULL` uses all cells in the project.
#' @param bin_size Width of the coverage bins in bp. 500 (the default) reads the
#'   project's existing TileMatrix directly. Any other value re-bins the
#'   fragments overlapping the plotted window on the fly; this is region-local,
#'   so nothing is written back into the Arrow files.
#' @param normalize If TRUE, apply deepTools-style RPGC ("1x genomic coverage")
#'   scaling per group: counts are multiplied by
#'   `effective_genome_size / (n_fragments * mean_fragment_length)`, so a value
#'   of 1 corresponds to one-fold genome coverage. Effective genome size
#'   defaults to the project's chromSizes minus blacklisted bp.
#' @param genome_size Effective genome size in bp for RPGC. Defaults to the
#'   project's chromSizes total minus the blacklist.
#' @param frag_len_n Cells sampled per group when estimating mean fragment
#'   length for the RPGC denominator.
#' @param show_peaks If TRUE, draw an AnnotationTrack of the ArchRProject's
#'   PeakMatrix intervals (`getPeakSet()`) overlapping the window, coloured by
#'   the peak set's `peakType` metadata (Promoter / Distal / Exonic / Intronic).
#' @param peak_col Named vector of fill colours for `peakType`. Types absent
#'   from this vector fall back to grey.
#' @param highlight If TRUE, mark the queried locus (`start`-`end`, without the
#'   flanking window) across all tracks with a HighlightTrack.
#' @param highlight_col,highlight_fill Border and fill for that marker. The fill
#'   is drawn behind the tracks, so keep it pale.
#' @param track_sizes Relative heights for coverage:peaks:genes. Defaults to
#'   `c(7, 1, 2)`; the ideogram and axis are given small fixed shares on top.
#' @param cov_type Gviz DataTrack type for the coverage panel. Defaults to
#'   "polygon" (a filled area plot), which supports a y-axis. Note that
#'   "horizon" ignores `showAxis` entirely, so no RPGC scale can be drawn with
#'   it; with a single fill colour per group the two look the same anyway.
#' @param show_axis If TRUE, draw a y-axis on the coverage panel. All groups
#'   share `ylim`, so the axis is enabled on the first overlaid track only
#'   (Gviz draws a single axis per OverlayTrack).
#' @param ylim Optional numeric length-2 y-axis range shared by all coverage
#'   tracks. Defaults to `c(0, max * (1 + ylim_pad))` over the plotted bins,
#'   which keeps groups on a common scale and stops tall peaks being clipped
#'   at the top of the panel. Pass an explicit range to compare loci.
#' @param ylim_pad Fractional headroom added above the observed maximum when
#'   `ylim` is not supplied.
#' @param alpha Opacity in [0, 1] applied to the overlaid coverage horizon
#'   fills. Lower values make overlaps easier to see through; higher values
#'   make individual groups easier to distinguish.
#' @param legend If TRUE, draw a colour key for the `slot` groups in a strip
#'   below the tracks.
#' @param legend_height Fraction of the device height reserved for that strip.
#' @param label_upstream If TRUE (default), place each transcript's gene symbol
#'   at its upstream end -- left of plus-strand genes, right of minus-strand
#'   ones -- rather than always to the left as Gviz does, and clamp the choice
#'   to whichever side keeps the label inside the plotted window. This stops
#'   labels being clipped at either edge.
#' @param txdb,orgdb Annotation packages used to build the GeneRegionTrack.
#' @param main_title Main title for the plot. A single character string passed
#'   to `plotTracks(main = ...)`; `NULL` (default) draws no title. It is stored
#'   on the returned track list, so `replot_gviz_tracks()` reuses it.
#' @param cytoband_file Path to a local UCSC `cytoBandIdeo` table (gz or plain)
#'   used to draw the IdeogramTrack offline. If missing, Gviz falls back to
#'   querying UCSC.
#' @param force If TRUE, ignore the cache and re-extract the TileMatrix.
#' @return A list of Gviz tracks (with the plot window recorded in the
#'   attributes "from", "to", "chromosome"), replottable with
#'   `Gviz::plotTracks()`.
plot_gviz_pileups_by_category <-
  function(
    chr,
    start,
    end,
    window = 50000,
    ref_ArchR_obj,
    main_title = NULL,
    slot = "category",
    barcodes = NULL,
    bin_size = 500,
    normalize = TRUE,
    genome_size = NULL,
    frag_len_n = 300L,
    alpha = 0.5,
    ylim = NULL,
    ylim_pad = 0.05,
    show_axis = TRUE,
    cov_type = c("polygon", "horizon"),
    show_peaks = TRUE,
    peak_col = c(
      Promoter = "#E7298A",
      Distal = "#7570B3",
      Exonic = "#66A61E",
      Intronic = "#E6AB02"
    ),
    highlight = TRUE,
    highlight_col = "red",
    highlight_fill = "#FFE9E9",
    track_sizes = c(7, 1, 2),
    legend = TRUE,
    legend_height = NULL,
    legend_rows = NULL,
    label_upstream = TRUE,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    orgdb = org.Hs.eg.db::org.Hs.eg.db,
    genome = "hg38",
    cytoband_file = paste0("cytoBandIdeo_", genome, ".txt.gz"),
    force = FALSE,
    plot = TRUE
  ) {
    stopifnot(
      inherits(ref_ArchR_obj, "ArchRProject"),
      length(chr) == 1,
      is.numeric(start),
      is.numeric(end),
      is.numeric(window),
      is.numeric(bin_size),
      length(bin_size) == 1,
      bin_size > 0,
      is.numeric(alpha),
      length(alpha) == 1,
      alpha >= 0,
      alpha <= 1,
      is.null(ylim) || (is.numeric(ylim) && length(ylim) == 2),
      is.numeric(ylim_pad),
      length(ylim_pad) == 1,
      ylim_pad >= 0
    )
    cov_type <- match.arg(cov_type)
    stopifnot(is.numeric(track_sizes), length(track_sizes) == 3)
    stopifnot(
      is.null(main_title) ||
        (is.character(main_title) && length(main_title) == 1)
    )
    if (!slot %in% colnames(ref_ArchR_obj@cellColData)) {
      stop("`slot` '", slot, "' not found in cellColData.")
    }
    if (!length(ArchR::getArrowFiles(ref_ArchR_obj))) {
      stop("No Arrow files found in this ArchRProject.")
    }

    from <- max(1, floor(min(start, end) - window))
    to <- ceiling(max(start, end) + window)

    # ---- cells to use -------------------------------------------------------
    all_cells <- rownames(ref_ArchR_obj@cellColData)
    if (is.null(barcodes)) {
      cells_use <- all_cells
    } else {
      cells_use <- intersect(all_cells, as.character(barcodes))
      if (!length(cells_use)) {
        stop("None of the supplied `barcodes` are present in ref_ArchR_obj.")
      }
      if (length(cells_use) < length(unique(barcodes))) {
        message(
          length(unique(barcodes)) - length(cells_use),
          " supplied barcodes were not found in the project and were dropped."
        )
      }
    }

    groups <- as.character(ref_ArchR_obj@cellColData[cells_use, slot])
    keep <- !is.na(groups)
    cells_use <- cells_use[keep]
    groups <- groups[keep]
    grp_levels <- sort(unique(groups))

    # ---- tile counts (cached) ----------------------------------------------
    cache_key <- paste(
      chr,
      from,
      to,
      slot,
      bin_size,
      normalize,
      length(cells_use),
      # cheap fingerprint of the cell set
      if (length(cells_use)) {
        paste(range(cells_use), collapse = "|")
      } else {
        ""
      },
      sep = "_"
    )

    if (!force && !is.null(.tile_cache[[cache_key]])) {
      mat_grp <- .tile_cache[[cache_key]]$mat
      tile_gr <- .tile_cache[[cache_key]]$gr
    } else {
      arrow_files <- ArchR::getArrowFiles(ref_ArchR_obj)
      bin_size <- as.integer(bin_size)

      if (bin_size == 500L) {
        # ---- read the existing 500 bp TileMatrix ----------------------------
        if (!"TileMatrix" %in% ArchR::getAvailableMatrices(ref_ArchR_obj)) {
          stop("No TileMatrix found in this ArchRProject.")
        }
        feat <- ArchR:::.getFeatureDF(arrow_files, "TileMatrix")
        feat <- feat[as.character(feat$seqnames) == chr, , drop = FALSE]
        if (!nrow(feat)) {
          stop("Chromosome '", chr, "' not present in the TileMatrix.")
        }
        tile_size <- as.integer(diff(sort(head(feat$start, 2))))
        feat <- feat[
          feat$start + tile_size > from & feat$start <= to,
          ,
          drop = FALSE
        ]
        if (!nrow(feat)) {
          stop("No TileMatrix bins overlap the requested window.")
        }
        bin_starts <- feat$start

        per_arrow <- .arrow_lapply(arrow_files, function(af) {
          cells_af <- intersect(
            cells_use,
            ArchR:::.availableCells(af, "TileMatrix")
          )
          if (!length(cells_af)) {
            return(NULL)
          }
          m <- tryCatch(
            ArchR:::.getMatFromArrow(
              ArrowFile = af,
              featureDF = feat,
              binarize = FALSE,
              useMatrix = "TileMatrix",
              cellNames = cells_af
            ),
            # some Arrow files store TileMatrix binarized on disk
            error = function(e) {
              ArchR:::.getMatFromArrow(
                ArrowFile = af,
                featureDF = feat,
                binarize = TRUE,
                useMatrix = "TileMatrix",
                cellNames = cells_af
              )
            }
          )
          list(mat = m, cells = colnames(m))
        })
      } else {
        # ---- re-bin fragments in the window on the fly ----------------------
        # Deliberately region-local: a genome-wide addFeatureMatrix at small
        # bin sizes would write tens of millions of features per cell back into
        # the Arrow files for what is only ever plotted one window at a time.
        tile_size <- bin_size
        bin_starts <- seq(
          floor(from / bin_size) * bin_size,
          to,
          by = bin_size
        )
        region <- GenomicRanges::GRanges(
          chr,
          IRanges::IRanges(from, to)
        )

        per_arrow <- .arrow_lapply(arrow_files, function(af) {
          cells_af <- intersect(
            cells_use,
            ArchR:::.availableCells(af, "TileMatrix")
          )
          if (!length(cells_af)) {
            return(NULL)
          }
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
            return(NULL)
          }
          fr <- IRanges::subsetByOverlaps(fr, region)
          if (!length(fr)) {
            return(NULL)
          }
          # count each fragment in the bin holding its midpoint
          mids <- GenomicRanges::start(fr) +
            (GenomicRanges::width(fr) - 1) / 2
          bin_idx <- findInterval(mids, bin_starts)
          cell_ids <- as.character(S4Vectors::mcols(fr)$RG)
          ok <- bin_idx >= 1 & bin_idx <= length(bin_starts) & !is.na(cell_ids)
          if (!any(ok)) {
            return(NULL)
          }
          cf <- factor(cell_ids[ok], levels = cells_af)
          m <- Matrix::sparseMatrix(
            i = bin_idx[ok],
            j = as.integer(cf),
            x = 1,
            dims = c(length(bin_starts), length(cells_af)),
            dimnames = list(NULL, cells_af)
          )
          list(mat = m, cells = cells_af)
        })
      }

      # ---- accumulate per-group sums -------------------------------------
      n_bins <- length(bin_starts)
      mat_grp <- matrix(
        0,
        nrow = n_bins,
        ncol = length(grp_levels),
        dimnames = list(NULL, grp_levels)
      )
      n_cells_grp <- setNames(
        integer(length(grp_levels)),
        grp_levels
      )

      for (res in per_arrow) {
        if (is.null(res)) {
          next
        }
        g <- groups[match(res$cells, cells_use)]
        for (lv in grp_levels) {
          idx <- which(g == lv)
          if (!length(idx)) {
            next
          }
          mat_grp[, lv] <- mat_grp[, lv] +
            Matrix::rowSums(res$mat[, idx, drop = FALSE])
          n_cells_grp[lv] <- n_cells_grp[lv] + length(idx)
        }
      }

      if (normalize) {
        # deepTools RPGC ("1x genomic coverage"):
        #   scale = effective_genome_size / (n_fragments * mean_frag_length)
        egs <- genome_size %||% .effective_genome_size(ref_ArchR_obj)
        frag_len <- .estimate_frag_length(
          arrow_files,
          cells_use,
          groups,
          grp_levels,
          chr = chr,
          n_sample = frag_len_n
        )
        n_frags_grp <- vapply(
          grp_levels,
          function(lv) {
            sum(as.numeric(
              ref_ArchR_obj@cellColData[cells_use[groups == lv], "nFrags"]
            ))
          },
          numeric(1)
        )
        scale_grp <- egs / (n_frags_grp * frag_len)
        message(
          "RPGC scaling: effective genome size ",
          format(egs, big.mark = ","),
          " bp, mean fragment length ",
          round(frag_len, 1),
          " bp; scale factors ",
          paste(
            paste0(grp_levels, "=", signif(scale_grp, 3)),
            collapse = ", "
          )
        )
        for (lv in grp_levels) {
          mat_grp[, lv] <- mat_grp[, lv] * scale_grp[[lv]]
        }
      }

      tile_gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(
          start = bin_starts,
          width = tile_size
        )
      )
      .tile_cache[[cache_key]] <- list(mat = mat_grp, gr = tile_gr)
    }

    # ---- tracks -------------------------------------------------------------
    # Use the locally cached cytoband table when available so no UCSC query is
    # made; fall back to Gviz's online lookup only if the file is missing.
    bands <- load_cytobands(genome = genome, file = cytoband_file)
    ideo_track <- tryCatch(
      if (is.null(bands)) {
        Gviz::IdeogramTrack(genome = genome, chromosome = chr)
      } else {
        Gviz::IdeogramTrack(
          genome = genome,
          chromosome = chr,
          bands = bands[bands$chrom == chr, , drop = FALSE]
        )
      },
      error = function(e) {
        message(
          "IdeogramTrack unavailable (",
          conditionMessage(e),
          "); skipped."
        )
        NULL
      }
    )

    axis_track <- Gviz::GenomeAxisTrack(
      range = GenomicRanges::GRanges(
        chr,
        IRanges::IRanges(from, to)
      ),
      add53 = TRUE,
      add35 = TRUE
    )

    # "Dark2" carries 8 colours; recycle if there are more groups than that
    pal <- setNames(
      rep_len(
        RColorBrewer::brewer.pal(
          max(3L, min(8L, length(grp_levels))),
          "Dark2"
        ),
        length(grp_levels)
      ),
      grp_levels
    )

    # Shared y-limits across groups. Gviz otherwise derives ylim per DataTrack,
    # so overlaid horizon tracks end up on different scales AND the tallest
    # peaks sit flush against the panel top; pad the maximum to leave headroom.
    y_max <- suppressWarnings(max(mat_grp, na.rm = TRUE))
    if (!is.finite(y_max) || y_max <= 0) {
      y_max <- 1
    }
    ylim_use <- ylim %||% c(0, y_max * (1 + ylim_pad))

    # "horizon.scale" is the width of ONE horizon band, not the whole range:
    # values above it are folded into darker bands, which with a single fill
    # colour renders as a flat-topped, clipped peak. Set it to the full padded
    # range so the tallest peak stays inside the first band.
    horizon_scale <- diff(ylim_use)

    cov_tracks <- lapply(seq_along(grp_levels), function(i) {
      lv <- grp_levels[[i]]
      fill_i <- grDevices::adjustcolor(pal[[lv]], alpha.f = alpha)
      Gviz::DataTrack(
        range = tile_gr,
        data = matrix(mat_grp[, lv], nrow = 1),
        genome = genome,
        chromosome = chr,
        name = if (normalize) "Coverage (RPGC)" else "Coverage",
        type = cov_type,
        # horizon fills
        fill.horizon = rep(fill_i, 6),
        col.horizon = NA,
        horizon.origin = 0,
        horizon.scale = horizon_scale,
        # polygon ("mountain") fills
        fill.mountain = rep(fill_i, 2),
        col.mountain = NA,
        ylim = ylim_use,
        # groups share ylim, so one axis suffices; Gviz draws the axis of the
        # first track in an OverlayTrack
        showAxis = isTRUE(show_axis) && i == 1L,
        cex.axis = 0.7,
        alpha = alpha,
        alpha.title = 1,
        legend = TRUE,
        groups = lv
      )
    })
    cov_overlay <- Gviz::OverlayTrack(
      trackList = cov_tracks,
      name = paste0("Coverage (", slot, ")")
    )
    # OverlayTrack does not forward these, so set them on the composite too
    Gviz::displayPars(cov_overlay) <- list(
      showAxis = isTRUE(show_axis),
      ylim = ylim_use
    )

    # ---- called peaks from the PeakMatrix -----------------------------------
    peak_track <- NULL
    peak_palette <- NULL
    if (isTRUE(show_peaks)) {
      ps <- tryCatch(
        ArchR::getPeakSet(ref_ArchR_obj),
        error = function(e) NULL
      )
      if (is.null(ps) || !length(ps)) {
        message("No peak set found in this ArchRProject; peak track skipped.")
      } else {
        ps <- IRanges::subsetByOverlaps(
          ps,
          GenomicRanges::GRanges(chr, IRanges::IRanges(from, to))
        )
        if (!length(ps)) {
          message("No called peaks overlap the window; peak track skipped.")
        } else {
          ptype <- as.character(S4Vectors::mcols(ps)$peakType)
          if (is.null(ptype)) {
            ptype <- rep("Peak", length(ps))
          }
          ptype[is.na(ptype)] <- "Other"
          types <- unique(ptype)
          # Gviz colours AnnotationTrack features by matching feature() values
          # against same-named displayPars entries.
          type_cols <- vapply(
            types,
            function(tt) unname(peak_col[tt]) %NA% "grey50",
            character(1)
          )
          peak_track <- Gviz::AnnotationTrack(
            range = GenomicRanges::GRanges(
              seqnames = chr,
              ranges = IRanges::ranges(ps)
            ),
            genome = genome,
            chromosome = chr,
            name = "Peaks",
            feature = ptype,
            stacking = "dense",
            col = NA,
            showFeatureId = FALSE
          )
          Gviz::displayPars(peak_track) <- as.list(type_cols)
          peak_palette <- type_cols[order(names(type_cols))]
        }
      }
    }

    gene_track <- Gviz::GeneRegionTrack(
      txdb,
      chromosome = chr,
      start = from,
      end = to,
      genome = genome,
      name = "Genes",
      stacking = "squish",
      # cap stacks at 1 row/gene so dense loci don't overflow the device
      collapseTranscripts = "longest"
    )
    # map Entrez gene ids -> HGNC symbols for transcriptAnnotation = "symbol"
    if (length(Gviz::gene(gene_track))) {
      sym <- suppressMessages(tryCatch(
        AnnotationDbi::mapIds(
          orgdb,
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

    # Gviz justifies group labels for the whole track at once (`just.group`),
    # always relative to the transcript's LEFT end. A transcript that starts
    # left of the window therefore has its label drawn off-screen and clipped
    # (e.g. minus-strand KLF10 here); the mirror case clips on the right.
    #
    # Decide each transcript's side by whether the label actually FITS inside
    # the window, preferring the upstream end, then bucket transcripts into a
    # left-justified and a right-justified track and overlay the two.
    if (isTRUE(label_upstream)) {
      gr_all <- gene_track@range
      mc <- S4Vectors::mcols(gr_all)
      if (length(gr_all) && !is.null(mc$transcript)) {
        tx <- as.character(mc$transcript)
        labs <- if (!is.null(mc$symbol)) {
          as.character(mc$symbol)
        } else {
          tx
        }

        # per-transcript extent and label
        tx_start <- tapply(GenomicRanges::start(gr_all), tx, min)
        tx_end <- tapply(GenomicRanges::end(gr_all), tx, max)
        tx_lab <- tapply(labs, tx, function(x) x[[1]])
        tx_strand <- tapply(
          as.character(GenomicRanges::strand(gr_all)),
          tx,
          function(x) x[[1]]
        )

        # estimate label width in bp: character width at the group cex,
        # scaled from the plotted window and the panel's width in inches
        cex_grp <- suppressWarnings(as.numeric(
          Gviz::displayPars(gene_track)[["cex.group"]]
        ))
        if (!is.finite(cex_grp)) {
          cex_grp <- 0.6
        }
        panel_in <- max(
          1,
          grid::convertWidth(
            grid::unit(1, "npc"),
            "inches",
            valueOnly = TRUE
          )
        )
        bp_per_in <- (to - from) / panel_in
        lab_bp <- vapply(
          as.character(tx_lab),
          function(s) {
            grid::convertWidth(
              grid::stringWidth(s),
              "inches",
              valueOnly = TRUE
            ) *
              cex_grp *
              bp_per_in
          },
          numeric(1)
        )

        # a left label occupies [start - width, start]; a right label
        # occupies [end, end + width]
        fits_left <- (tx_start - lab_bp) >= from
        fits_right <- (tx_end + lab_bp) <= to
        prefer_right <- tx_strand == "-"

        side <- ifelse(
          prefer_right & fits_right,
          "right",
          ifelse(
            !prefer_right & fits_left,
            "left",
            # preferred side clips: fall back to whichever fits, else draw the
            # label over the transcript body ("above"), which is always inside
            # the window -- this covers transcripts spanning the whole view.
            ifelse(
              fits_left,
              "left",
              ifelse(fits_right, "right", "above")
            )
          )
        )
        names(side) <- names(tx_start)

        tx_side <- side[tx]
        make_side <- function(which_side, just) {
          sub_gr <- gr_all[tx_side == which_side]
          if (!length(sub_gr)) {
            return(NULL)
          }
          g <- Gviz::GeneRegionTrack(
            sub_gr,
            chromosome = chr,
            genome = genome,
            name = "Genes",
            stacking = "squish"
          )
          Gviz::displayPars(g) <- list(just.group = just)
          g
        }
        sides <- Filter(
          Negate(is.null),
          list(
            make_side("left", "left"),
            make_side("right", "right"),
            make_side("above", "above")
          )
        )
        if (length(sides) > 1L) {
          gene_track <- Gviz::OverlayTrack(
            trackList = sides,
            name = "Genes"
          )
        } else if (length(sides) == 1L) {
          gene_track <- sides[[1]]
        }
      }
    }

    # Keep sizes aligned with whichever tracks actually exist. The ideogram and
    # axis get small fixed shares; coverage:peaks:genes follow track_sizes.
    parts <- list(
      list(track = ideo_track, size = 0.6),
      list(track = axis_track, size = 0.8),
      list(track = cov_overlay, size = track_sizes[[1]]),
      list(track = peak_track, size = track_sizes[[2]]),
      list(track = gene_track, size = track_sizes[[3]])
    )
    parts <- Filter(function(p) !is.null(p$track), parts)
    tracks <- lapply(parts, `[[`, "track")
    track_size_vec <- vapply(parts, `[[`, numeric(1), "size")

    # Mark the queried locus across the data tracks. The ideogram and axis are
    # left outside the highlight so it does not span the karyotype cartoon.
    if (isTRUE(highlight)) {
      is_data <- vapply(
        tracks,
        function(tr) {
          !inherits(tr, "IdeogramTrack") && !inherits(tr, "GenomeAxisTrack")
        },
        logical(1)
      )
      if (any(is_data)) {
        hl <- Gviz::HighlightTrack(
          trackList = tracks[is_data],
          start = min(start, end),
          end = max(start, end),
          chromosome = chr,
          col = highlight_col,
          fill = highlight_fill,
          inBackground = TRUE
        )
        tracks <- c(tracks[!is_data], list(hl))
        # plotTracks() expands a HighlightTrack into its children when it
        # validates `sizes`, so keep one entry per wrapped track, in order.
        track_size_vec <- c(
          track_size_vec[!is_data],
          track_size_vec[is_data]
        )
      }
    }

    attr(tracks, "chromosome") <- chr
    attr(tracks, "from") <- from
    attr(tracks, "to") <- to
    attr(tracks, "palette") <- pal
    attr(tracks, "alpha") <- alpha
    attr(tracks, "slot") <- slot
    attr(tracks, "sizes") <- track_size_vec
    attr(tracks, "peak_palette") <- peak_palette
    attr(tracks, "main_title") <- main_title

    if (plot) {
      .plot_tracks_with_legend(
        tracks,
        from = from,
        to = to,
        chr = chr,
        pal = pal,
        alpha = alpha,
        slot = slot,
        legend = legend,
        legend_height = legend_height,
        legend_rows = legend_rows,
        sizes = track_size_vec,
        peak_pal = peak_palette,
        main_title = main_title
      )
    }

    invisible(tracks)
  }

# test plot
# tr <- plot_gviz_pileups_by_category(
#   chr = df_sig_snp_list$seqnames[1],
#   start = df_sig_snp_list$start[1],
#   end = df_sig_snp_list$end[1],
#   bin_size = 50,
#   window = 1000,
#   ref_ArchR_obj = projHepatocytes,
#   main_title = paste0(
#     df_sig_snp_list$variantID[1],
#     " (",
#     df_sig_snp_list$SYMBOL[1],
#     ")",
#     ",",
#     df_sig_snp_list$annotation[1]
#   ),
#   slot = "category",
#   alpha = 0.85
# )

writeout_dir <- "gviz_macrophage_SNP_pileups"
if (!dir.exists(writeout_dir)) {
  dir.create(writeout_dir, recursive = TRUE)
}

for (i in seq_along(df_sig_snp_list$seqnames)) {
  snp_chr <- df_sig_snp_list$seqnames[i]
  snp_start <- df_sig_snp_list$start[i]
  snp_end <- df_sig_snp_list$end[i]
  snp_id <- df_sig_snp_list$variantID[i]
  gene_symbol <- df_sig_snp_list$SYMBOL[i]
  gene_annotation <- df_sig_snp_list$annotation[i]
  main_title <-
    paste0(
      snp_id,
      " (",
      gene_symbol,
      ")",
      ",",
      gene_annotation
    )

  tryCatch(
    {
      tr <- plot_gviz_pileups_by_category(
        chr = snp_chr,
        start = snp_start,
        end = snp_end,
        bin_size = 50,
        window = 2000,
        ref_ArchR_obj = ref_ArchR_obj,
        slot = "category",
        alpha = 0.85,
        main_title = main_title
      )

      pdf(
        file.path(
          writeout_dir,
          paste0(snp_id, "_", gene_symbol, "_gviz_plot.pdf")
        ),
        width = 11,
        height = 8.5
      )
      replot_gviz_tracks(tr)
      dev.off()
    },
    error = function(e) {
      # close any partially-opened device so the next iteration isn't broken
      if (dev.cur() != 1) {
        dev.off()
      }
      message("Skipped ", snp_id, " (", gene_symbol, "): ", conditionMessage(e))
    }
  )
}
