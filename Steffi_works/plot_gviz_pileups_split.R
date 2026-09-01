#! /usr/bin/env Rscript

# Split the hepatocyte ASoC SNP table into N roughly-even parts for the
# scatter/gather MPI pipeline (see plot_gviz_pileups_mpi_worker.R).
#
# Each part (except possibly the last) holds a MULTIPLE OF 4 SNPs so that every
# chunk fills whole 2x2 pages when the panels are later combined 4-per-page.
# The final part takes the remainder.
#
# Outputs, under <writeout_dir>/parts/:
#   part_01.tsv ... part_NN.tsv   one SNP sub-table per job-array element
#   manifest.txt                  single integer: the number of parts written
#
# Run standalone (the guardian script calls it for you):
#   Rscript plot_gviz_pileups_split.R [n_parts]

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)

# ---- parameters ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
n_parts <- if (length(args) >= 1L) as.integer(args[[1]]) else 8L
stopifnot(is.finite(n_parts), n_parts >= 1L)

input_tsv <- "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv"
writeout_dir <- "gviz_hepatocyte_SNP_pileups"
parts_dir <- file.path(writeout_dir, "parts")
panels_per_page <- 4L # keep chunk sizes a multiple of this

# ---- read ------------------------------------------------------------------
# quote = "" / comment.char = "": fields in the motif / TF columns contain
# apostrophes and quotes; with R's default quoting those rows merge and only
# ~148 of the 315 records parse. Disabling quote handling reads all 315.
df_sig_snp_list <- read.table(
  input_tsv,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
n_snps <- nrow(df_sig_snp_list)
if (!n_snps) {
  stop("No SNPs found in ", input_tsv)
}

# ---- chunk size: ceil(n / n_parts) rounded UP to a multiple of 4 -----------
# Rounding the chunk (not the count) up to a multiple of 4 keeps every part but
# the last aligned to whole 2x2 pages. The actual number of parts written is
# then ceil(n_snps / chunk), which may be <= n_parts.
base_chunk <- ceiling(n_snps / n_parts)
chunk <- ceiling(base_chunk / panels_per_page) * panels_per_page
chunk <- max(chunk, panels_per_page)

starts <- seq(1L, n_snps, by = chunk)
actual_parts <- length(starts)

# ---- (re)create the parts directory ----------------------------------------
if (dir.exists(parts_dir)) {
  old <- list.files(
    parts_dir,
    pattern = "^(part_\\d+\\.tsv|chunk_\\d+\\.pdf|manifest\\.txt)$",
    full.names = TRUE
  )
  if (length(old)) {
    file.remove(old)
  }
} else {
  dir.create(parts_dir, recursive = TRUE)
}

# ---- write the parts -------------------------------------------------------
for (p in seq_len(actual_parts)) {
  i0 <- starts[[p]]
  i1 <- min(i0 + chunk - 1L, n_snps)
  part_df <- df_sig_snp_list[i0:i1, , drop = FALSE]
  out_tsv <- file.path(parts_dir, sprintf("part_%02d.tsv", p))
  write.table(
    part_df,
    file = out_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  message(sprintf(
    "part %02d: %d SNPs (rows %d-%d) -> %s",
    p,
    i1 - i0 + 1L,
    i0,
    i1,
    out_tsv
  ))
}

# ---- manifest (number of parts) so the array / guardian can size itself ----
writeLines(as.character(actual_parts), file.path(parts_dir, "manifest.txt"))

message(sprintf(
  "Wrote %d part(s) of up to %d SNPs each (%d SNPs total) into %s",
  actual_parts,
  chunk,
  n_snps,
  parts_dir
))
