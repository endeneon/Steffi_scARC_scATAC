#! /usr/bin/env Rscript

# init ####
{
  library(BiocParallel)
  library(parallel)
  library(future)

  library(readr)

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

# load data tables ####
df_wnt <-
  read.table(
    "RNASeq_tables/fpkm_by_condition_and_treatment_wnt.txt",
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )

df_wnt <-
  df_wnt[-c(4, 5, 9, 15, 19, 21), ]

df_2_plot <-
  reshape2::melt(
    df_wnt,
    id.vars = c("Condition", "Treatment", "group")
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
