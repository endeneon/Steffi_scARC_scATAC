#! /usr/bin/env Rscript

# Build a GRanges object from the ENCODE hg38 blacklist BED file.
# BED coordinates are 0-based, half-open [start, end); GRanges are 1-based,
# closed [start, end], so we add 1 to the start when converting.

library(GenomicRanges)
library(GenomeInfoDb)

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/R_ASoC_analysis"
)

bed_file <- "hg38-blacklist.v2.bed"

# The 4th column ("Low Mappability" / "High Signal Region") contains spaces, so
# read as strictly tab-delimited to keep it as a single field.
blacklist_df <-
  read.table(
    bed_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    col.names = c("chrom", "start", "end", "reason")
  )

blacklist_gr <-
  GenomicRanges::GRanges(
    seqnames = blacklist_df$chrom,
    ranges = IRanges::IRanges(
      start = blacklist_df$start + 1L, # BED 0-based -> GRanges 1-based
      end = blacklist_df$end
    ),
    strand = "*",
    reason = blacklist_df$reason
  )

# Tag the genome build for downstream compatibility (e.g. ArchR / overlap ops).
GenomeInfoDb::genome(blacklist_gr) <- "hg38"

# Sort by chromosome then coordinate for tidy, reproducible output.
blacklist_gr <- GenomeInfoDb::sortSeqlevels(blacklist_gr)
blacklist_gr <- sort(blacklist_gr)

print(blacklist_gr)
print(paste0("Created GRanges with ", length(blacklist_gr), " blacklist regions."))

# Optionally persist the object for reuse.
saveRDS(blacklist_gr, file = "hg38_blacklist_v2_granges.rds")
