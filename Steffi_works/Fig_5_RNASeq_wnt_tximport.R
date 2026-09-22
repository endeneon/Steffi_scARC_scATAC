#! /usr/bin/env Rscript

# init ####
{
  library(BiocParallel)
  library(parallel)
  library(future)

  library(readr)

  library(tximport)
  library(edgeR)

  library(Matrix)

  library(dplyr)
  library(reshape2)
  library(stringr)
  library(magrittr)
  library(qs2)

  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)

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


# use this conda env
# /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/standalone_conda_envs/r45_py312_scARC

# set working dir
working_dir <- "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
setwd(working_dir)

# determine if R is running in RSTUDIO or POSITRON, and set the future plan accordingly
# Manual set
# Sys.setenv(VSCODE = "1")
if (
  Sys.getenv("RSTUDIO") == "1" ||
    Sys.getenv("POSITRON") == "1" ||
    Sys.getenv("VSCODE") == "1" ||
    Sys.getenv("TERM_PROGRAM") == "vscode"
) {
  print("Running under IDE, use plan(multisession)")
  session_plan <- "multisession"
} else {
  print("Running under Rscript, use plan(multicore)")
  session_plan <- "multicore"
}

# preload functions for future ####
get_available_workers <-
  function(x) {
    future::plan(session_plan) # check here!
    return(future::nbrOfFreeWorkers())
  }

## nthreads
if (session_plan == "multisession") {
  workers_2_use <-
    min(
      get_available_workers(1) - 1,
      16
    )
} else {
  workers_2_use <-
    min(
      get_available_workers(1) - 1,
      20
    )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  set.seed(42)
  options(future.globals.maxSize = workers_2_use * 20 * 1024^3) # 40 G per thread
  future::plan(
    session_plan, # Do NOT use "multisession" here if submitting LSF jobs, use "multicore" instead
    workers = workers_2_use
  )
}

# load raw tables ####

input_files_list <-
  list.files(
    path = "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/Steffi_sync/RNAseq_29Jul2026/outputBAMs",
    pattern = "*\\.RSEM\\.isoforms\\.results$",
    recursive = T,
    full.names = T
  )
names(input_files_list) <-
  basename(input_files_list) |>
  str_remove("\\.RSEM\\.isoforms\\.results$") |>
  str_remove("X__")

# RSEM isoform files carry both transcript_id and gene_id; build tx2gene so
# tximport can summarize to the gene level (otherwise summarizeFail() errors).
tx2gene <-
  read.table(
    input_files_list[[1]],
    header = TRUE,
    sep = "\t"
  )[, c("transcript_id", "gene_id")]

df_raw <-
  tximport(
    files = input_files_list,
    type = "rsem",
    tx2gene = tx2gene,
    txIn = TRUE,
    txOut = FALSE,
    countsFromAbundance = "lengthScaledTPM"
  )

df_lookup_table <-
  read.table(
    "RNAseq_set_3/outputBAMs/RNAseq-29Jul2026_RSEM_gene_count.2026-08-01_16-27-05.txt",
    header = TRUE,
    sep = "\t"
  )

df_raw_rsem_counts <-
  df_raw$counts
df_raw_rsem_counts <-
  df_raw_rsem_counts[
    rowSums(df_raw_rsem_counts) > 0,
  ]

# Map rownames from Ensembl geneID to geneSymbol, aggregating (summing) counts
# across the 120 symbols shared by multiple gene IDs.
sym <-
  df_lookup_table$geneSymbol[
    match(rownames(df_raw_rsem_counts), df_lookup_table$geneID)
  ]
df_raw_rsem_counts <-
  rowsum(df_raw_rsem_counts, group = sym)


df_meta <-
  read.table(
    "RNAseq_set_3/input_config/meta.txt",
    header = TRUE,
    sep = "\t"
  )

df_DGE <-
  DGEList(
    counts = df_raw_rsem_counts,
    samples = df_meta$ID,
    group = df_meta$GROUP
  )

df_DGE <-
  calcNormFactors(df_DGE, method = "TMM")

df_cpm <-
  as.data.frame(cpm(
    df_DGE,
    normalized.lib.sizes = TRUE,
    log = FALSE
  ))

genes_2_plot <-
  c(
    "APC2",
    "WNT5A",
    "WNT8B",
    "EN1",
    "LMX1A",
    "LMX1B",
    "AXIN2",
    "LEF1",
    "DKK1",
    "MSX1",
    "SNAI1",
    "SNAI2"
  )

genes_2_plot <-
  c(
    "TH",
    "SLC6A3",
    "SLC18A2",
    "AADC",
    "LMX1A",
    "LMX1B",
    "EN1",
    "FOXA2",
    "OTX2",
    "SHH",
    "NR4A2",
    "NURR1",
    "ASCL1",
    "NEUROG2"
  )


cpm_writeout <-
  df_cpm[
    rownames(df_cpm) %in% genes_2_plot,
  ]
gene_list <- rownames(cpm_writeout)
df_cpm <-
  cpm_writeout[, c(1:ncol(cpm_writeout))]
df_cpm <-
  as.data.frame(t(df_cpm))
colnames(df_cpm) <- gene_list

df_grouping <-
  as.data.frame(str_split(df_meta$GROUP, pattern = "_", simplify = TRUE))
colnames(df_grouping) <- c("Condition", "Treatment")
df_grouping$Group <- df_meta$GROUP

df_notrim <- cbind(df_grouping, df_cpm)
df_trim <-
  df_notrim[
    !rownames(df_notrim) %in% c("a18", "a22", "a31", "a34"),
  ]

df_2_plot <-
  reshape2::melt(
    df_trim,
    id.vars = c("Condition", "Treatment", "Group")
  )

colnames(df_2_plot) <-
  c("Condition", "Treatment", "Group", "Gene", "CPM")
df_2_plot$Group <-
  factor(
    df_2_plot$Group,
    levels = c(
      "ctrl_noChir",
      "KO2_noChir",
      "KO3_noChir",
      "ctrl_Chir",
      "KO2_Chir",
      "KO3_Chir"
    )
  )

df_2_plot |>
  summarise(
    CPM_mean = mean(CPM),
    CPM_sem = sd(CPM) / sqrt(n()),
    .by = c(Condition, Treatment, Group, Gene)
  ) |>
  ggplot(
    aes(
      x = Group,
      y = CPM_mean,
      fill = Condition,
      alpha = Treatment,
      shape = Condition
    )
  ) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_errorbar(
    aes(ymin = CPM_mean - CPM_sem, ymax = CPM_mean + CPM_sem),
    position = position_dodge(width = 0.9),
    width = 0.3,
    alpha = 1
  ) +
  geom_point(
    data = df_2_plot,
    aes(y = CPM),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9
    ),
    size = 0.5,
    color = "black"
  ) +
  scale_shape_manual(values = c("ctrl" = 1, "KO2" = 2, "KO3" = 5)) +
  scale_fill_manual(
    values = c("ctrl" = "steelblue", "KO2" = "orange", "KO3" = "darkred")
  ) +
  scale_alpha_manual(values = c("noChir" = 0.5, "Chir" = 1)) +
  scale_y_continuous(
    name = "Mean CPM",
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic() +
  facet_wrap(~Gene, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 315, hjust = 0))
