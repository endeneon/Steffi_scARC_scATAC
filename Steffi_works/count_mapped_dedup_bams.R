#! /usr/bin/env Rscript

# init ####

{
  library(Rsamtools)
  library(doParallel)
  library(foreach)
}

# Count mapped, deduplicated reads for each BAM in df_bam_list ####
# - isUnmappedQuery = FALSE keeps only mapped reads
# - isDuplicate = FALSE drops reads flagged as PCR/optical duplicates
# Assumes df_bam_list and workers_2_use already exist in the session.

count_param <-
  ScanBamParam(
    flag = scanBamFlag(
      isUnmappedQuery = FALSE,
      isDuplicate = FALSE
    )
  )

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

cl <-
  makeCluster(cluster_size)
registerDoParallel(cl)

bam_counts <-
  foreach(
    bam_file = df_bam_list,
    .combine = rbind,
    .packages = "Rsamtools"
  ) %dopar%
  {
    res <-
      countBam(
        bam_file,
        param = count_param
      )

    data.frame(
      bam_file = bam_file,
      mapped_dedup_records = res$records,
      stringsAsFactors = FALSE
    )
  }

stopCluster(cl)

rownames(bam_counts) <-
  NULL

print(bam_counts)
