#! /usr/bin/env Rscript

# one-line conclusion: all 3 methods are biased.
# either CFRSB or CF would work (CFRSB wins a bit if only considering FDR-significant ones),
# CR does not work and produced very high bias after FDR correction

# For count pattern _CFRSB: ASoC percRef > 0.5: 14901 of 44599 variants. of which are also FDR significant: 506 of 1280 variants.
# For count pattern _CR: ASoC percRef > 0.5: 16688 of 49897 variants. of which are also FDR significant: 987 of 2687 variants.
# For count pattern _CF: ASoC percRef > 0.5: 14902 of 44605 variants. of which are also FDR significant: 496 of 1259 variants.

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

  library(GenomicRanges)

  library(dplyr)
  library(stringr)
  library(magrittr)
  library(qs2)

  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(paletteer)
  library(EnhancedVolcano)

  library(Gviz)

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
working_dir = "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
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

# count bam sizes (only need to do once)

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

print(paste0(
  "Counting BAM records in parallel using ",
  cluster_size,
  " workers..."
))
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

# define all count pattern to be tested
all_count_patterns <-
  c(
    "_CFRSB",
    "_CR",
    "_CF"
  )

# file to collect the per-pattern summaries
summary_output_file <-
  file.path(data_dir, "benchmark_count_summary.txt")
# truncate / create the file before the loop so each run starts fresh
cat("", file = summary_output_file, append = FALSE)

# loop through all count patterns and save the count results
for (count_pattern in all_count_patterns) {
  print(paste0("Processing count pattern: ", count_pattern))

  df_tsv_list <-
    list.files(
      path = ASoC_output_dir,
      pattern = paste0(count_pattern, "\\.tsv\\.gz$"),
      recursive = T,
      full.names = T
    )
  df_tsv_list <-
    sort(df_tsv_list)

  tsv_list_4_match <-
    data.frame(
      tsv_file = df_tsv_list,
      sample_id = str_split_i(
        df_tsv_list,
        pattern = "/",
        i = -1
      ) %>%
        str_remove(pattern = paste0(count_pattern, "\\.tsv\\.gz$")) %>%
        str_remove(pattern = "_ASE.*$")
    )

  if (!all(tsv_list_4_match$sample_id %in% bam_counts$sample_id)) {
    stop(
      "Not all sample IDs from TSV list are present in BAM counts. Check the sample_id extraction logic."
    )
  }

  main_tsv_bam_table <-
    merge(
      tsv_list_4_match,
      bam_counts,
      by = "sample_id"
    )
  main_tsv_bam_table$inverse_size <-
    mean(main_tsv_bam_table$mapped_dedup_records) /
    main_tsv_bam_table$mapped_dedup_records

  print(
    "Finished counting BAM records and matching with TSV files. Now reading TSV files and normalizing counts..."
  )
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
        dplyr::select(1:7)
      return(as.data.frame(df_raw))
    }

  # read in all tsv files and assemble a list of data frames, with names as sample IDs
  # (reuses the cluster created above)
  print("Reading TSV files in parallel...")
  tsv_df_list <-
    foreach(
      tsv_file = main_tsv_bam_table$tsv_file,
      .packages = "readr"
    ) %dopar%
    {
      read_tsv_gz_file(tsv_file)
    }

  # stopCluster(cl)
  names(tsv_df_list) <-
    main_tsv_bam_table$sample_id

  normalised_df_tsv_list <-
    lapply(
      names(tsv_df_list),
      function(sample_id) {
        df <- tsv_df_list[[sample_id]]
        df$inverse_size <- main_tsv_bam_table$inverse_size[
          main_tsv_bam_table$sample_id == sample_id
        ]
        df$norm_refCount <- df$refCount * df$inverse_size
        df$norm_altCount <- df$altCount * df$inverse_size
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
            )
          ]
        df_2_return <-
          df_2_return[
            (df_2_return$norm_refCount > 0) & (df_2_return$norm_altCount > 0),
          ]
        return(df_2_return)
      }
    )
  names(normalised_df_tsv_list) <-
    names(tsv_df_list)

  # full outer merge across all samples on one_line_index ####
  # row-bind every data frame, then collapse on one_line_index:
  #   - indices unique to a single source df are carried over unchanged
  #   - indices shared across source dfs have their norm_refCount / norm_altCount
  #     summed (na.rm = TRUE)
  # this is equivalent to Reduce(merge(all.x = T, all.y = T)) + summing matched
  # counts, but scales to any number of data frames in one pass.
  print("Merging normalized counts across samples...")
  merged_normalised_df <-
    dplyr::bind_rows(normalised_df_tsv_list) %>%
    dplyr::group_by(one_line_index) %>%
    dplyr::summarise(
      norm_refCount = sum(norm_refCount, na.rm = TRUE),
      norm_altCount = sum(norm_altCount, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  merged_normalised_df$norm_refCount <-
    round(merged_normalised_df$norm_refCount)
  merged_normalised_df$norm_altCount <-
    round(merged_normalised_df$norm_altCount)
  merged_normalised_df$norm_sumCount <-
    merged_normalised_df$norm_refCount + merged_normalised_df$norm_altCount

  # sum(
  #   merged_normalised_df$norm_refCount >= 2 &
  #     merged_normalised_df$norm_altCount >= 2
  # ) # 75479
  # mean(merged_normalised_df$norm_sumCount) # 13.84079
  # median(merged_normalised_df$norm_sumCount) # 3

  merged_normalised_df <-
    merged_normalised_df[
      merged_normalised_df$norm_refCount >= 2 &
        merged_normalised_df$norm_altCount >= 2,
    ]
  # mean(merged_normalised_df$norm_sumCount) # 34.46173
  # median(merged_normalised_df$norm_sumCount) # 12
  # hist(merged_normalised_df$norm_sumCount, breaks = 500, xlim = c(0, 100)) # most have low counts, but some have very high counts (up to 400+)

  merged_normalised_df <-
    merged_normalised_df[
      merged_normalised_df$norm_sumCount >= 10,
    ]
  # mean(merged_normalised_df$norm_sumCount) # 54.13875
  # median(merged_normalised_df$norm_sumCount) # 24

  print("calculating p-values and percentages...")
  merged_normalised_df$pVal <-
    mapply(
      function(refCount, altCount) {
        totalCount <- refCount + altCount
        if (totalCount == 0) {
          return(NA)
        }
        # Perform binomial test against null hypothesis of 0.5
        test_result <- binom.test(altCount, totalCount, p = 0.5)
        return(test_result$p.value)
      },
      merged_normalised_df$norm_refCount,
      merged_normalised_df$norm_altCount
    )
  merged_normalised_df$percRef <-
    merged_normalised_df$norm_refCount / merged_normalised_df$norm_sumCount
  # sum(merged_normalised_df$percRef > 0.5) # 14901
  merged_normalised_df$percAlt <-
    merged_normalised_df$norm_altCount / merged_normalised_df$norm_sumCount
  # sum(merged_normalised_df$percAlt > 0.5) # 26083

  merged_normalised_df$adjPVal <-
    p.adjust(
      merged_normalised_df$pVal,
      method = "fdr"
    )
  print("Finished calculating p-values and percentages.")

  merged_normalised_df <-
    merged_normalised_df[
      !is.na(merged_normalised_df$adjPVal),
    ]

  summary_message <-
    paste0(
      "For count pattern ",
      count_pattern,
      ": ",
      "ASoC percRef > 0.5: ",
      sum(merged_normalised_df$percRef > 0.5),
      " of ",
      nrow(merged_normalised_df),
      " variants.\n",
      "of which are also FDR significant: ",
      sum(
        merged_normalised_df$percRef > 0.5 & merged_normalised_df$adjPVal < 0.05
      ),
      " of ",
      sum(merged_normalised_df$adjPVal < 0.05),
      " variants.\n",
      "#################################\n"
    )

  print(summary_message)

  # append this iteration's summary to the shared output file
  cat(
    summary_message,
    "\n",
    file = summary_output_file,
    append = TRUE
  )
}

stopCluster(cl)
