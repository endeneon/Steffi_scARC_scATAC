#! /usr/bin/env Rscript

# MPI worker for ONE part of the hepatocyte ASoC SNP set.
#
# Launched (once per job-array element) as:
#   mpirun -n <n_workers + 1> --bind-to none \
#     Rscript plot_gviz_pileups_mpi_worker.R --part-index <k> --threads <t>
#
# Two nested levels of parallelism:
#   * OUTER  - doMPI ranks spread this part's SNPs across <n_workers> workers.
#              (SPMD launch: rank 0 = master, ranks 1..N = compute workers.)
#   * INNER  - within each SNP, a node-local fork pool of <t> cores lets ArchR
#              read the ~30 Arrow files in parallel (that is what makes
#              addArchRThreads(t) with t > 1 actually pay off; with t = 1 each
#              SNP reads its Arrow files serially and is slow).
# With the default 10 workers x 6 threads that is 60 cores, which fits on one
# 60-core node (the job script pins it there with span[hosts=1]).
#
# Because we launch SPMD, every rank runs the code BEFORE startMPIcluster(),
# including the source() below -- so every worker already has the helper
# functions and the loaded `projHepatocytes` in its own memory. We therefore do
# NOT ship that (large) ArchRProject over MPI (`.noexport = "projHepatocytes"`);
# each worker uses its local copy.
#
# The master (rank 0) collects the per-SNP Gviz track lists and draws this
# part's panels, 4 per landscape US-Letter page (11 in wide x 8.5 in tall), into
# gviz_hepatocyte_SNP_pileups/parts/chunk_<k>.pdf. The guardian later merges the
# per-part chunk PDFs into the master PDF.

# ---- parse args ------------------------------------------------------------
# Minimal flag parser (no optparse dependency). Recognised:
#   --part-index <int>  (falls back to $LSB_JOBINDEX, then 1)
#   --threads <int>     ArchR / inner-fork threads per worker (default 6)
#   --n-parts <int>     total parts, for messaging only (optional)
.get_flag <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, a)
  if (!is.na(hit) && hit < length(a)) a[[hit + 1L]] else default
}

part_index <- as.integer(
  .get_flag("--part-index", Sys.getenv("LSB_JOBINDEX", "1"))
)
if (!is.finite(part_index) || part_index < 1L) {
  part_index <- 1L
}
threads_per_worker <- as.integer(.get_flag("--threads", "6"))
if (!is.finite(threads_per_worker) || threads_per_worker < 1L) {
  threads_per_worker <- 1L
}

# ---- load the shared functions + projHepatocytes (lib-only) ----------------
# gviz.pipeline.lib_only = TRUE makes plot_gviz_pileups_by_category.R stop after
# defining its functions and loading `projHepatocytes`, i.e. it does NOT run its
# own (serial / PSOCK) driver + PDF block. This runs on EVERY MPI rank.
options(gviz.pipeline.lib_only = TRUE)
source("plot_gviz_pileups_by_category.R", local = FALSE)

suppressPackageStartupMessages({
  library(doMPI)
  library(foreach)
})

# ---- this part's SNPs ------------------------------------------------------
writeout_dir <- "gviz_hepatocyte_SNP_pileups"
parts_dir <- file.path(writeout_dir, "parts")
part_tsv <- file.path(parts_dir, sprintf("part_%02d.tsv", part_index))
chunk_pdf <- file.path(parts_dir, sprintf("chunk_%02d.pdf", part_index))
# Completion signature for resume. Only written after the PDF is fully closed on
# disk (see the end of this script), so its mere presence certifies the chunk
# finished. Removed up-front so an interrupted run never leaves a stale marker.
chunk_done <- file.path(parts_dir, sprintf("chunk_%02d.done", part_index))
suppressWarnings(if (file.exists(chunk_done)) file.remove(chunk_done))

if (!file.exists(part_tsv)) {
  stop("Part file not found: ", part_tsv, " (did the split step run?)")
}
# md5 of the part file, recorded in the signature so a later run re-does this
# chunk if the underlying part changed (the array script / guardian compare it).
part_md5 <- unname(tools::md5sum(part_tsv))
# quote = "" / comment.char = "": the motif / TF columns contain apostrophes
# and quotes that would otherwise merge rows (matches the split step's read).
df_part <- read.table(
  part_tsv,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
n_snps <- nrow(df_part)

# ---- start the MPI cluster from the ranks mpirun launched ------------------
# No `count`: doMPI turns the existing mpirun ranks into the cluster (rank 0
# master, the rest workers). Worker ranks block inside startMPIcluster() and
# never reach the code below; only the master continues.
cl <- doMPI::startMPIcluster()
doMPI::registerDoMPI(cl)
n_workers <- doMPI::clusterSize(cl)

message(sprintf(
  "[part %02d] %d SNPs across %d MPI worker(s), %d ArchR thread(s) each.",
  part_index,
  n_snps,
  n_workers,
  threads_per_worker
))

# packages every worker must have attached to build a track list. (The helper
# functions and projHepatocytes are already present on each worker from the
# source() above; projHepatocytes is explicitly NOT re-shipped.)
export_pkgs <- c(
  "ArchR",
  "Gviz",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "Matrix",
  "matrixStats",
  "RColorBrewer",
  "AnnotationDbi",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "BSgenome.Hsapiens.UCSC.hg38",
  "doParallel"
)

# ---- build every SNP's track list in parallel ------------------------------
tr_list <-
  foreach(
    i = seq_len(n_snps),
    .errorhandling = "pass",
    .noexport = "projHepatocytes",
    .packages = export_pkgs
  ) %dopar%
  {
    # fresh evaluation context on the worker: re-assert plotting options and
    # ArchR settings, then register a NODE-LOCAL fork pool so .arrow_lapply()
    # reads the Arrow files `threads_per_worker`-at-a-time (this is the inner
    # level that makes threads > 1 matter). registerDoParallel(cores = t) uses
    # forking on Linux, so it stays on this rank's node.
    options(useUCSCChromosomeNames = FALSE)
    ArchR::addArchRThreads(threads = threads_per_worker, force = TRUE)
    ArchR::addArchRGenome("hg38")
    doParallel::registerDoParallel(cores = threads_per_worker)

    main_title <-
      paste0(
        df_part$variantID[i],
        " (",
        df_part$SYMBOL[i],
        "),",
        df_part$annotation[i]
      )

    tryCatch(
      plot_gviz_pileups_by_category(
        chr = df_part$seqnames[i],
        start = df_part$start[i],
        end = df_part$end[i],
        bin_size = 50,
        window = 2000,
        ref_ArchR_obj = projHepatocytes,
        slot = "category",
        alpha = 0.85,
        main_title = main_title,
        plot = FALSE
      ),
      error = function(e) {
        message(
          "Skipped ",
          df_part$variantID[i],
          " (",
          df_part$SYMBOL[i],
          "): ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }

doMPI::closeCluster(cl)

# ---- keep only valid track lists -------------------------------------------
# Real track lists carry a "chromosome" attr; drop NULLs (failed builds) and any
# error objects foreach returned with .errorhandling = "pass".
tr_list <-
  Filter(
    function(tr) !is.null(attr(tr, "chromosome", exact = TRUE)),
    tr_list
  )

n_ok <- length(tr_list)
message(sprintf(
  "[part %02d] built %d / %d track list(s); writing %s",
  part_index,
  n_ok,
  n_snps,
  chunk_pdf
))

# ---- render this part: 4 panels (2x2) per landscape US-Letter page ---------
# Landscape Letter = 11 in wide x 8.5 in tall.
panels_per_page <- 4L
n_pages <- max(1L, ceiling(n_ok / panels_per_page))

pdf(chunk_pdf, width = 11, height = 8.5)
on.exit(grDevices::dev.off(), add = TRUE)

if (!n_ok) {
  grid::grid.newpage()
  grid::grid.text(sprintf(
    "part %02d: no track lists could be built",
    part_index
  ))
} else {
  for (pg in seq_len(n_pages)) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))

    idx <-
      seq(
        (pg - 1L) * panels_per_page + 1L,
        min(pg * panels_per_page, n_ok)
      )

    for (k in seq_along(idx)) {
      row <- ((k - 1L) %/% 2L) + 1L
      col <- ((k - 1L) %% 2L) + 1L
      grid::pushViewport(grid::viewport(
        layout.pos.row = row,
        layout.pos.col = col
      ))
      # new_page = FALSE keeps the panel inside the current grid cell instead of
      # advancing the device to a fresh page.
      replot_gviz_tracks(tr_list[[idx[k]]], new_page = FALSE)
      grid::popViewport()
    }

    grid::popViewport()
  }
}

# dev.off() runs via on.exit(); mpi.quit() (from Rmpi) shuts the MPI runtime
# down cleanly and exits R -- doMPI does not re-export it.
grDevices::dev.off()
on.exit() # clear the handler we already ran
message(sprintf("[part %02d] done: %s", part_index, chunk_pdf))

# ---- completion signature (resume) -----------------------------------------
# Written only now, after the device is closed and the PDF exists non-empty, so
# the signature can never precede a finished chunk. The array script / guardian
# skip any part whose signature is present, whose PDF is non-empty, and whose
# recorded part_md5 still matches the current part file.
if (file.exists(chunk_pdf) && file.info(chunk_pdf)$size > 0) {
  writeLines(
    c(
      sprintf("part_index=%d", part_index),
      sprintf("n_snps=%d", n_snps),
      sprintf("n_built=%d", n_ok),
      sprintf("part_md5=%s", part_md5),
      sprintf("pdf=%s", basename(chunk_pdf)),
      sprintf("completed_at=%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
    ),
    chunk_done
  )
  message(sprintf("[part %02d] wrote signature %s", part_index, chunk_done))
} else {
  message(sprintf(
    "[part %02d] PDF missing/empty; NOT writing signature (will re-run).",
    part_index
  ))
}

Rmpi::mpi.quit(save = "no")
