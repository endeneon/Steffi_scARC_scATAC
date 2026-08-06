#! /usr/bin/env Rscript

# Siwei 24 Mar 2026
# init ####
{
  library(BiocParallel)
  library(glmGamPoi)
  library(future)
  library(parallel)
  library(foreach)
  library(doParallel)

  library(snow)
  library(harmony)
  # library(optparse)

  library(dplyr)
  library(data.table)

  library(edgeR)

  # library(scuttle)
  library(Matrix)
  library(matrixStats)

  library(Seurat)
  library(Signac)
  # library(loomR)
  # library(anndata)
  library(hdf5r)
  library(arrow)
  # library(rhdf5)

  library(stringr)

  library(qs2)
  library(fs)

  # library(SingleR)

  library(ggplot2)
  library(gplots)
  library(patchwork)
  library(scCustomize)

  # library(scuttle)

  # library(scMerge)
  # library(scater)
}

# setwd(
#   "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/
#    szhang37/projects/szhang_dev/STEREO_seq/Human_liver/R_liver"
# )
setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)
# determine if R is running in RSTUDIO/VSCode/Positron
# Gate on interactive(): under Rscript/LSF this is FALSE, so a leaked
# TERM_PROGRAM=vscode from the submit shell won't force memory-duplicating
# multisession (each worker copies the full object -> OOM under LSF).
if (interactive() &&
    (Sys.getenv("RSTUDIO") == "1" || Sys.getenv("TERM_PROGRAM") == "vscode")) {
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
print(available_cores)

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

if (session_plan == "multisession") {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        16
      )
    )
} else {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        32
      )
    )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
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

# plot the genes by sample in pseudobulk ####
merged_integrated_liver_obj <-
  qs_read(
    "merged_integrated_seurat_obj.qs2",
    nthreads = 8
  )

DefaultAssay(merged_integrated_liver_obj) <- "RNA"
merged_integrated_liver_obj <- NormalizeData(
  merged_integrated_liver_obj,
  verbose = TRUE
)
merged_integrated_liver_obj <- ScaleData(
  merged_integrated_liver_obj,
  verbose = TRUE
)

qs_save(
  merged_integrated_liver_obj,
  "merged_integrated_seurat_obj.qs2",
  nthreads = 8
)

pseudobulk_liver <-
  AggregateExpression(
    merged_integrated_liver_obj,
    group.by = c("orig.ident"),
    assays = "RNA",
    slot = "data"
  )

pseudobulk_liver <- as.data.frame(pseudobulk_liver)
head(rownames(pseudobulk_liver))

genes_2_plot <-
  read.table(
    "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated_HCC_crossref_DisGeNET_yes.tsv",
    header = T,
    sep = "\t"
  )
genes_2_plot <- genes_2_plot$SYMBOL
genes_2_plot <- genes_2_plot[!duplicated(genes_2_plot)]

df_2_plot <-
  pseudobulk_liver[rownames(pseudobulk_liver) %in% genes_2_plot, ]
df_2_plot <-
  as.data.frame(t(df_2_plot))
rownames(df_2_plot)
colnames(df_2_plot)

df_long <-
  df_2_plot |>
  tibble::rownames_to_column("sample_name") |>
  data.table::as.data.table() |>
  data.table::melt(
    id.vars = "sample_name",
    variable.name = "gene_name",
    value.name = "exp_level"
  )

head(df_long)
df_long$category <- "Primary"
df_long$category[str_detect(df_long$sample_name, "pre")] <- "Resistant"
unique(df_long$sample_name)

p_box <-
  ggplot(
    df_long,
    aes(
      x = gene_name,
      y = exp_level,
      fill = category
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.size = 0,
    alpha = 0.5
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.05,
      dodge.width = 0.7
    ),
    size = 1,
    alpha = 1
  ) +
  # dashed horizontal line at the mean of each box
  stat_summary(
    fun = mean,
    geom = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    position = position_dodge(width = 0.8),
    width = 1,
    # linetype = "dashed",
    colour = "darkred"
  ) +
  scale_fill_manual(values = c("orange4", "darkblue")) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315, hjust = 0))

ggsave(
  "pseudobulk_liver_gene_boxplot.pdf",
  plot = p_box,
  width = 10,
  height = 6
)
