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

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)
projMacrophages <-
  ArchR::loadArchRProject(path = "ArchR_macrophages")
snp_annotation_file <-
  "sig_ASoC_by_celltype/sig_ASoC_in_Macrophage_annotated_PICALM.tsv"

df_snp_raw <-
  read.table(
    snp_annotation_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
df_snp_raw <-
  df_snp_raw[!is.na(df_snp_raw$ENSEMBL), ]

df_2_plot <-
  df_snp_raw[, c("seqnames", "start", "variantID", "ENSEMBL", "SYMBOL")]
colnames(df_2_plot) <-
  c("chr", "pos", "snp_id", "ensembl_id", "gene_symbol")

# RPGC of +/- 50 bp SNP windows by sample ####
# RPGC = mean per-bp fragment coverage / (N_frags * mean_frag_len / G_eff)
# N_frags and mean_frag_len are computed per sample over all fragments from
# the sample's cells in projMacrophages.
{
  flank_bp <- 50
  effective_genome_size <- 2701495761

  gr_windows <-
    GRanges(
      seqnames = df_2_plot$chr,
      ranges = IRanges(
        start = df_2_plot$pos - flank_bp,
        end = df_2_plot$pos + flank_bp
      )
    )

  cell_meta <- as.data.frame(getCellColData(projMacrophages))
  arrow_files <- getArrowFiles(projMacrophages)
  sample_names <- sort(unique(cell_meta$Sample))
  all_chrs <- getChromSizes(projMacrophages)
  all_chrs <- as.character(seqnames(all_chrs))

  rpgc_list <-
    lapply(sample_names, function(i_sample) {
      i_cells <- rownames(cell_meta)[cell_meta$Sample == i_sample]
      n_frags <- 0
      sum_frag_len <- 0
      window_bp_cov <- numeric(length(gr_windows))

      for (i_chr in all_chrs) {
        frags <-
          getFragmentsFromArrow(
            ArrowFile = arrow_files[[i_sample]],
            chr = i_chr,
            cellNames = i_cells,
            verbose = FALSE
          )
        if (length(frags) == 0) {
          next
        }
        n_frags <- n_frags + length(frags)
        sum_frag_len <- sum_frag_len + sum(as.numeric(width(frags)))

        i_win <- which(as.character(seqnames(gr_windows)) == i_chr)
        if (length(i_win) == 0) {
          next
        }
        hits <- findOverlaps(gr_windows[i_win], frags, ignore.strand = TRUE)
        if (length(hits) == 0) {
          next
        }
        overlap_bp <-
          width(pintersect(
            gr_windows[i_win][queryHits(hits)],
            frags[subjectHits(hits)],
            ignore.strand = TRUE
          ))
        window_bp_cov[i_win] <-
          window_bp_cov[i_win] +
          as.numeric(tapply(
            X = overlap_bp,
            INDEX = factor(queryHits(hits), levels = seq_along(i_win)),
            FUN = sum,
            default = 0
          ))
      }

      mean_cov <- window_bp_cov / width(gr_windows)
      scale_factor <- n_frags * (sum_frag_len / n_frags) / effective_genome_size
      print(paste0(
        i_sample,
        ": ",
        length(i_cells),
        " cells, ",
        n_frags,
        " fragments, mean length ",
        round(sum_frag_len / n_frags, 1),
        " bp"
      ))
      mean_cov / scale_factor
    })

  df_rpgc <- as.data.frame(do.call(cbind, rpgc_list))
  colnames(df_rpgc) <- sample_names

  df_2_plot_rpgc <- cbind(df_2_plot, df_rpgc)

  # per-sample cell counts, row order matches the RPGC columns
  df_sample_stats <-
    data.frame(
      sample = sample_names,
      n_cells = as.integer(table(cell_meta$Sample)[sample_names]),
      n_frags = as.numeric(tapply(cell_meta$nFrags, cell_meta$Sample, sum)[
        sample_names
      ])
    )
  write.table(
    df_sample_stats,
    file = sub("\\.tsv$", "_RPGC_sample_stats.tsv", snp_annotation_file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    df_2_plot_rpgc,
    file = sub("\\.tsv$", "_RPGC_by_sample.tsv", snp_annotation_file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
