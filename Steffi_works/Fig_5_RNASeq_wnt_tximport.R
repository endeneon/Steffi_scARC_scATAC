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

### Wnt effectors ####
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

### Dopaminergic neuron markers ####
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

df_2_plot$Treatment <-
  factor(
    df_2_plot$Treatment,
    levels = c("noChir", "Chir")
  )
df_2_plot$Condition <-
  factor(
    df_2_plot$Condition,
    levels = c("ctrl", "KO2", "KO3")
  )

## w/o significance #####

df_2_plot |>
  summarise(
    CPM_mean = mean(CPM),
    CPM_sem = sd(CPM) / sqrt(n()),
    .by = c(Condition, Treatment, Group, Gene)
  ) |>
  ggplot(
    aes(
      x = Treatment,
      y = CPM_mean,
      fill = Treatment,
      alpha = Condition,
      shape = Condition,
      group = Condition
    )
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("noChir" = "orange", "Chir" = "pink")
  ) +
  geom_errorbar(
    aes(ymin = CPM_mean - CPM_sem, ymax = CPM_mean + CPM_sem),
    position = position_dodge(width = 0.7),
    width = 0.3,
    alpha = 1
  ) +
  geom_point(
    data = df_2_plot,
    aes(y = CPM),
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.7
    ),
    size = 1,
    color = "black",
    alpha = 1
  ) +
  scale_shape_manual(values = c("ctrl" = 1, "KO2" = 2, "KO3" = 5)) +
  scale_alpha_manual(values = c("ctrl" = 0, "KO2" = 0.5, "KO3" = 1)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(
    name = "Mean CPM",
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic() +
  facet_wrap(~Gene, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

## w/ significance #####
{
  # significance stars from p value
  sig_labeller <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  }

  # all pairwise two-sample t-tests among the 6 Groups, within each Gene
  compute_pairwise_t <- function(dat) {
    groups <- levels(factor(dat$Group))
    res <- data.frame(
      group1 = character(),
      group2 = character(),
      p = numeric()
    )
    for (pair in combn(groups, 2, simplify = FALSE)) {
      v1 <- dat$CPM[dat$Group == pair[1]]
      v2 <- dat$CPM[dat$Group == pair[2]]
      if (length(v1) < 2 || length(v2) < 2) {
        next
      }
      p <- tryCatch(t.test(v1, v2)$p.value, error = function(e) NA_real_)
      res <- rbind(res, data.frame(group1 = pair[1], group2 = pair[2], p = p))
    }
    res
  }

  # numeric x position of each dodged column on the discrete Treatment axis
  dodge_w <- 0.7
  cond_levels <- levels(df_2_plot$Condition)
  treat_levels <- levels(df_2_plot$Treatment)
  n_cond <- length(cond_levels)
  cond_offset <- (seq_len(n_cond) - (n_cond + 1) / 2) * dodge_w / n_cond

  group_pos <-
    expand.grid(
      Condition = cond_levels,
      Treatment = treat_levels,
      stringsAsFactors = FALSE
    ) |>
    mutate(
      Group = paste(Condition, Treatment, sep = "_"),
      x = match(Treatment, treat_levels) +
        cond_offset[match(Condition, cond_levels)]
    )

  # per-gene top of the data, used to stack the brackets
  y_base <-
    df_2_plot |>
    summarise(y_top = max(CPM, na.rm = TRUE), .by = Gene)

  stat_df <-
    df_2_plot |>
    group_by(Gene) |>
    group_modify(~ compute_pairwise_t(.x)) |>
    ungroup() |>
    filter(!is.na(p), p < 0.05) |>
    mutate(
      label = sig_labeller(p),
      xmin = group_pos$x[match(group1, group_pos$Group)],
      xmax = group_pos$x[match(group2, group_pos$Group)],
      xleft = pmin(xmin, xmax),
      xright = pmax(xmin, xmax),
      xmid = (xleft + xright) / 2
    ) |>
    left_join(y_base, join_by(Gene)) |>
    arrange(Gene, xleft, xright) |>
    mutate(
      y.position = y_top * (1.05 + 0.09 * (row_number() - 1)),
      .by = Gene
    )

  df_2_plot |>
    summarise(
      CPM_mean = mean(CPM),
      CPM_sem = sd(CPM) / sqrt(n()),
      .by = c(Condition, Treatment, Group, Gene)
    ) |>
    ggplot(
      aes(
        x = Treatment,
        y = CPM_mean,
        fill = Treatment,
        alpha = Condition,
        shape = Condition,
        group = Condition
      )
    ) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_col(
      position = position_dodge(width = 0.7),
      width = 0.65,
      color = "black"
    ) +
    scale_fill_manual(
      values = c("noChir" = "orange", "Chir" = "pink")
    ) +
    geom_errorbar(
      aes(ymin = CPM_mean - CPM_sem, ymax = CPM_mean + CPM_sem),
      position = position_dodge(width = 0.7),
      width = 0.3,
      alpha = 1
    ) +
    geom_point(
      data = df_2_plot,
      aes(y = CPM),
      position = position_jitterdodge(
        jitter.width = 0.25,
        dodge.width = 0.7
      ),
      size = 1,
      color = "black",
      alpha = 1
    ) +
    geom_segment(
      data = stat_df,
      aes(x = xleft, xend = xright, y = y.position, yend = y.position),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_text(
      data = stat_df,
      aes(x = xmid, y = y.position, label = label),
      inherit.aes = FALSE,
      vjust = -0.1,
      size = 3
    ) +
    scale_shape_manual(values = c("ctrl" = 1, "KO2" = 2, "KO3" = 5)) +
    scale_alpha_manual(values = c("ctrl" = 0, "KO2" = 0.5, "KO3" = 1)) +
    scale_x_discrete(expand = expansion(add = 0.4)) +
    scale_y_continuous(
      name = "Mean CPM",
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_classic() +
    facet_wrap(~Gene, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

###
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
      fill = Treatment,
      alpha = Treatment,
      shape = Condition
    )
  ) +
  geom_col(
    position = position_dodge(width = 1),
    width = 0.8
  ) +
  geom_errorbar(
    aes(ymin = CPM_mean - CPM_sem, ymax = CPM_mean + CPM_sem),
    position = position_dodge(width = 1),
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
    values = c("noChir" = "orange", "Chir" = "pink")
  ) +
  scale_alpha_manual(values = c("noChir" = 0.5, "Chir" = 1)) +
  scale_y_continuous(
    name = "Mean CPM",
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic() +
  facet_wrap(~Gene, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 315, hjust = 0))
