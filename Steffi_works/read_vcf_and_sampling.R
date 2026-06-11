#! /usr/bin/env Rscript

# init ####

{
  library(BiocParallel)
  library(parallel)
  library(future)

  library(doParallel)
  library(foreach)

  library(readr)
  library(vcfR)
  library(Rsamtools)

  library(Matrix)

  library(dplyr)
  library(stringr)
  library(magrittr)
  library(qs2)

  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)

  if (Sys.getenv("VSCODE") == "1") {
    library(languageserver)
    library(showtext)
    library(httpgd)

    showtext::showtext_auto()
  }
}

# set working dir
working_dir = "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
setwd(working_dir)

# determine if R is running in RSTUDIO or POSITRON, and set the future plan accordingly
if (Sys.getenv("RSTUDIO") == "1" | Sys.getenv("POSITRON") == "1") {
  print("Running under RStudio IDE, use plan(multisession)")
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

# main functions ####

# Count mapped, deduplicated reads for each BAM in df_bam_list ####
# - isUnmappedQuery = FALSE keeps only mapped reads
# - isDuplicate = FALSE drops reads flagged as PCR/optical duplicates
# Assumes df_bam_list and workers_2_use already exist in the session.
count_param <-
  Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE,
      isDuplicate = FALSE
    )
  )


ASoC_output_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/R_ASoC_analysis"
data_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/data_dir"

df_tsv_list <-
  list.files(
    path = ASoC_output_dir,
    pattern = "*.tsv.gz$",
    recursive = T,
    full.names = T
  )
df_tsv_list <-
  sort(df_tsv_list)

df_bam_list <-
  list.files(
    path = ASoC_output_dir,
    pattern = "*.bam$",
    recursive = T,
    full.names = T
  )
df_bam_list <-
  sort(df_bam_list)

# Verify index files up front (fast). Accept both naming conventions:
#   file.bam.bai  and  file.bai
bai_present <-
  file.exists(paste0(df_bam_list, ".bai")) |
  file.exists(sub("\\.bam$", ".bai", df_bam_list))

if (any(!bai_present)) {
  warning(
    "Missing index for ",
    sum(!bai_present),
    " BAM(s); indexing only those:\n",
    paste(df_bam_list[!bai_present], collapse = "\n")
  )
  invisible(lapply(df_bam_list[!bai_present], indexBam))
}

# cluster size = workers_2_use - 1, but never below 1
cluster_size <-
  max(
    workers_2_use - 1,
    1
  )

# Create the cluster once and reuse it for both foreach loops below;
# it is stopped after the final loop.
cl <-
  makeCluster(cluster_size)
doParallel::registerDoParallel(cl)

bam_counts <-
  foreach::foreach(
    bam_file = df_bam_list,
    .combine = rbind,
    .packages = "Rsamtools"
  ) %dopar%
  {
    res <-
      Rsamtools::countBam(
        bam_file,
        param = count_param
      )

    data.frame(
      bam_file = bam_file,
      mapped_dedup_records = res$records,
      stringsAsFactors = FALSE
    )
  }

rownames(bam_counts) <-
  str_split_i(
    bam_counts$bam_file,
    pattern = "/",
    i = -1
  ) %>%
  str_remove(pattern = ".bam$") %>%
  str_remove(pattern = "_bwa.*$")
bam_counts$sample_id <-
  rownames(bam_counts)

tsv_list_4_match <-
  data.frame(
    tsv_file = df_tsv_list,
    sample_id = str_split_i(
      df_tsv_list,
      pattern = "/",
      i = -1
    ) %>%
      str_remove(pattern = ".tsv.gz$") %>%
      str_remove(pattern = "_ASE.*$")
  )

all(tsv_list_4_match$sample_id %in% bam_counts$sample_id)

main_tsv_bam_table <-
  merge(
    tsv_list_4_match,
    bam_counts,
    by = "sample_id"
  )
main_tsv_bam_table$inverse_size <-
  mean(main_tsv_bam_table$mapped_dedup_records) /
  main_tsv_bam_table$mapped_dedup_records

read_tsv_gz_file <-
  function(tsv_file) {
    df_raw <-
      readr::read_tsv(
        tsv_file,
        col_types = cols(
          .default = col_character(),
          contig = col_character(),
          position = col_integer(),
          variantID = col_character(),
          refAllele = col_character(),
          altAllele = col_character(),
          refCount = col_integer(),
          altCount = col_integer(),
          totalCount = col_integer(),
          lowMAPQDepth = col_integer(),
          lowBaseQDepth = col_integer(),
          rawDepth = col_integer(),
          otherBases = col_logical(),
          improperPairs = col_logical()
        )
      )
    df_raw <-
      df_raw |>
      dplyr::filter(
        variantID != ".",
        !(contig %in% c("chrM", "chrX", "chrY")),
        !otherBases,
        !improperPairs,
        refCount > 0,
        altCount > 0
      ) |>
      dplyr::select(1:6)
    return(as.data.frame(df_raw))
  }

# read in all tsv files and assemble a list of data frames, with names as sample IDs
# (reuses the cluster created above)
tsv_df_list <-
  foreach(
    tsv_file = main_tsv_bam_table$tsv_file,
    .packages = "readr"
  ) %dopar%
  {
    read_tsv_gz_file(tsv_file)
  }

stopCluster(cl)
names(tsv_df_list) <-
  main_tsv_bam_table$sample_id

qs_save(
  tsv_df_list,
  file = file.path(data_dir, "tsv_df_list.qs2"),
  nthreads = 8
)
qs_save(
  main_tsv_bam_table,
  file = file.path(data_dir, "main_tsv_bam_table.qs2"),
  nthreads = 8
)

normalised_df_tsv_list <-
  lapply(
    names(tsv_df_list),
    function(sample_id) {
      df <- tsv_df_list[[sample_id]]
      df$inverse_size <- main_tsv_bam_table$inverse_size[main_tsv_bam_table$sample_id == sample_id]
      df$norm_refCount <- round(df$refCount * df$inverse_size)
      df$norm_altCount <- round(df$altCount * df$inverse_size)
      # df$norm_totalCount <- 
      #   df$norm_refCount + 
      #   df$norm_altCount
      df$one_line_index <- paste0(
        df$contig, 
        ":", 
        df$position, 
        ":", 
        df$variantID, 
        ":",
        df$refAllele, 
        ":", 
        df$altAllele
      )
      df_2_return <- 
      df[, 
      c(
        "one_line_index", 
        "norm_refCount", 
        "norm_altCount"
        )]
      return(df_2_return)
    }
  )
names(normalised_df_tsv_list) <-
  names(tsv_df_list)
qs_save(
  normalised_df_tsv_list,
  file = file.path(data_dir, "normalised_df_tsv_list.qs2"),
  nthreads = 8
)