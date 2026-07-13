#! /usr/bin/env Rscript

library(stringr)

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/R_ASoC_analysis"
)

# find all csv barcode files in the current working directory and its subdirectories
csv_barcode_files <-
  list.files(
    path = multiome_bams_dir,
    pattern = "per_barcode_metrics.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
names(csv_barcode_files) <-
  paste0(
    "X__",
    str_replace(
      str_replace(
        csv_barcode_files,
        pattern = paste0(multiome_bams_dir, "/"),
        replacement = ""
      ),
      pattern = "/outs/per_barcode_metrics.csv",
      replacement = ""
    )
  )

barcode_wo_suffix_dir <- "barcode_wo_suffix"
if (!dir.exists(barcode_wo_suffix_dir)) {
  dir.create(barcode_wo_suffix_dir)
}

for (each_csv_barcode_file in csv_barcode_files) {
  barcode_df <- read.csv(each_csv_barcode_file, header = FALSE)
  barcode_df <-
  str_split(barcode_df, pattern = "-", simplify = TRUE)[, 1]
  write.table(
    barcode_df,
    file = file.path(
        barcode_wo_suffix_dir, 
        paste0(
          str_replace(
            str_replace(
              each_csv_barcode_file,
              pattern = paste0(multiome_bams_dir, "/"),
              replacement = ""
            ),
            pattern = "/outs/per_barcode_metrics.csv",
            replacement = "_barcodes_wo_suffix.txt"
          )
        )
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}