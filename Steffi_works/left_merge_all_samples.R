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

options(ucscChromosomeNames = FALSE)
# constants ####
ASoC_output_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/R_ASoC_analysis"
data_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/data_dir"

tsv_df_list <-
  qs_read(
    file.path(
      data_dir,
      "tsv_df_list.qs2"
    ),
    nthreads = 4
  )

main_tsv_bam_table <-
  qs_read(
    file.path(
      data_dir,
      "main_tsv_bam_table.qs2"
    ),
    nthreads = 4
  )

# only take round after added up all reads
# hist(main_tsv_bam_table$mapped_dedup_records, breaks = 5)
# mean(main_tsv_bam_table$mapped_dedup_records)
# median(main_tsv_bam_table$mapped_dedup_records)

main_tsv_bam_table$median_inverse_size <-
  median(main_tsv_bam_table$mapped_dedup_records) /
  main_tsv_bam_table$mapped_dedup_records
# mean(main_tsv_bam_table$median_inverse_size) # 1.039025
# mean(main_tsv_bam_table$inverse_size) # 1.409165, use this one

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
sum(merged_normalised_df$percRef > 0.5) # 14901
merged_normalised_df$percAlt <-
  merged_normalised_df$norm_altCount / merged_normalised_df$norm_sumCount
sum(merged_normalised_df$percAlt > 0.5) # 26083


merged_normalised_df$adjPVal <-
  p.adjust(
    merged_normalised_df$pVal,
    method = "fdr"
  )

merged_normalised_df <-
  merged_normalised_df[
    !is.na(merged_normalised_df$adjPVal),
  ]

merged_normalised_df <-
  merged_normalised_df[
    order(merged_normalised_df$pVal),
  ]

ASoC_df <-
  as.data.frame(
    cbind(
      str_split(merged_normalised_df$one_line_index, ":", simplify = TRUE),
      merged_normalised_df$norm_refCount,
      merged_normalised_df$norm_altCount,
      merged_normalised_df$norm_sumCount,
      merged_normalised_df$percRef,
      merged_normalised_df$percAlt,
      merged_normalised_df$pVal,
      merged_normalised_df$adjPVal
    )
  )
colnames(ASoC_df) <-
  c(
    "chromosome",
    "position",
    "variantID",
    "refAllele",
    "altAllele",
    "norm_refCount",
    "norm_altCount",
    "norm_sumCount",
    "percRef",
    "percAlt",
    "pVal",
    "adjPVal"
  )

qs_save(
  ASoC_df,
  file.path(
    data_dir,
    "ASoC_df.qs2"
  ),
  nthreads = 4
)

# ASoC_df <-
# qs_read(
#   file.path(
#     data_dir,
#     "ASoC_df.qs2"
#   ),
#   nthreads = 4
# )
# write.table(
#   ASoC_df,
#   file.path(
#     data_dir,
#     "ASoC_df.tsv"
#   ),
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

# coerce ASoC_df into a GRangesList, one element per chromosome ####
# coerce metadata columns back to their proper types (cbind made all character)
ASoC_df$position <- as.integer(ASoC_df$position)
for (numeric_col in c(
  "norm_refCount",
  "norm_altCount",
  "norm_sumCount",
  "percRef",
  "pVal",
  "adjPVal"
)) {
  ASoC_df[[numeric_col]] <- as.numeric(ASoC_df[[numeric_col]])
}

# build a single GRanges, then split it by chromosome into a GRangesList
ASoC_Granges <-
  GenomicRanges::GRanges(
    seqnames = ASoC_df$chromosome,
    ranges = IRanges::IRanges(
      start = ASoC_df$position,
      end = ASoC_df$position
    ),
    strand = "*",
    ASoC_df[,
      c(
        "variantID",
        "refAllele",
        "altAllele",
        "norm_refCount",
        "norm_altCount",
        "norm_sumCount",
        "percRef",
        "pVal",
        "adjPVal"
      )
    ]
  )

# order seqlevels as chr1, chr2, ... chr21, chr22 so the split list follows suit
chr_order <- paste0("chr", 1:22)
chr_order <- chr_order[chr_order %in% GenomeInfoDb::seqlevels(ASoC_Granges)]
GenomeInfoDb::seqlevels(ASoC_Granges) <- chr_order

ASoC_GrangesList <-
  GenomicRanges::split(
    ASoC_Granges,
    GenomicRanges::seqnames(ASoC_Granges)
  )

qs_save(
  ASoC_GrangesList,
  file.path(
    data_dir,
    "ASoC_GrangesList.qs2"
  ),
  nthreads = 4
)

# genome-wide Manhattan plot with Gviz ####
# Gviz draws along a single coordinate axis, so the chromosomes are laid
# end-to-end on one synthetic chromosome and coloured per chromosome.

# flatten the GRangesList back into a single GRanges (chr1..chr22 order kept)
manhattan_gr <- unlist(ASoC_GrangesList)
names(manhattan_gr) <- NULL

# y axis value: 0 - log10(adjPVal)
manhattan_gr$negLogAdjPVal <- 0 - log10(manhattan_gr$adjPVal)

# chromosome factor in chr1, chr2, ... chr22 order
chr_levels <- paste0("chr", 1:22)
chr_levels <-
  chr_levels[
    chr_levels %in% as.character(GenomicRanges::seqnames(manhattan_gr))
  ]
chrom_factor <-
  factor(
    as.character(GenomicRanges::seqnames(manhattan_gr)),
    levels = chr_levels
  )

# cumulative per-chromosome offset so positions do not overlap across chrs
chrom_max <-
  tapply(
    GenomicRanges::end(manhattan_gr),
    chrom_factor,
    max
  )
chrom_offset <-
  cumsum(as.numeric(c(0, chrom_max[-length(chrom_max)])))
names(chrom_offset) <- chr_levels
cumulative_pos <-
  GenomicRanges::start(manhattan_gr) +
  chrom_offset[as.character(chrom_factor)]

# IRanges coordinates must stay below 2^31; the whole genome (~3e9 bp) exceeds
# that, so scale the synthetic coordinates down to fit (only relative spacing
# matters for the plot)
coord_scale_factor <-
  ceiling(max(cumulative_pos) / (2^31 - 2))
cumulative_pos <-
  as.integer(round(cumulative_pos / coord_scale_factor))

# place every variant on one synthetic chromosome with cumulative coordinates
manhattan_plot_gr <-
  GenomicRanges::GRanges(
    seqnames = "genome",
    ranges = IRanges::IRanges(
      start = cumulative_pos,
      width = 1
    ),
    strand = "*"
  )

# qualitative discrete palette with enough distinct colours for every
# chromosome (ggsci::default_igv has 51), no colour reused across chromosomes
chrom_colours <-
  paletteer::paletteer_d(
    "ggsci::default_igv",
    n = length(chr_levels)
  )

# Gviz's DataTrack `groups` argument groups data *rows* (samples), not points,
# so build one DataTrack per chromosome and overlay them; each track carries its
# own colour, giving the per-chromosome colouring with no colour reused.
manhattan_track_list <-
  lapply(
    seq_along(chr_levels),
    function(i) {
      keep <- chrom_factor == chr_levels[i]
      Gviz::DataTrack(
        range = manhattan_plot_gr[keep],
        data = manhattan_gr$negLogAdjPVal[keep],
        name = "-log10(adjPVal)",
        col = chrom_colours[i],
        ,
        cex = 0.5,
        type = "p"
      )
    }
  )

manhattan_track <-
  Gviz::OverlayTrack(
    trackList = manhattan_track_list,
    name = "-log10(adjPVal)"
  )

# FDR = 0.05 threshold on the y = 0 - log10(adjPVal) scale
fdr_threshold <- 0 - log10(0.05)

# per-chromosome label track: one feature spanning each chromosome's segment on
# the synthetic genome, labelled with the chromosome name (rotated 315 degrees
# so labels on small chromosomes do not overlap)
chrom_segment_starts <-
  as.integer(round(chrom_offset[chr_levels] / coord_scale_factor)) + 1L
chrom_segment_ends <-
  as.integer(round(
    (chrom_offset[chr_levels] + chrom_max[chr_levels]) / coord_scale_factor
  ))
chrom_label_track <-
  Gviz::AnnotationTrack(
    range = GenomicRanges::GRanges(
      seqnames = "genome",
      ranges = IRanges::IRanges(
        start = chrom_segment_starts,
        end = chrom_segment_ends
      ),
      strand = "*"
    ),
    id = chr_levels,
    name = "chr",
    featureAnnotation = "id",
    fontcolor.feature = "black",
    cex.feature = 0.7,
    rotation.item = 315,
    # rotation.title = 0,
    fill = "transparent",
    col = "transparent"
  )

# dev.off(which = which(names(dev.list()) == "httpgd"))
# graphics.off()

Gviz::plotTracks(
  list(
    manhattan_track,
    chrom_label_track
  ),
  chromosome = "genome",
  from = min(cumulative_pos),
  to = max(cumulative_pos) + 0.03 * diff(range(cumulative_pos)), # right pad so chr22 label is not cropped
  ylim = c(0, max(manhattan_gr$negLogAdjPVal, na.rm = TRUE)),
  baseline = fdr_threshold,
  col.baseline = "darkred",
  lty.baseline = 2,
  lwd.baseline = 1,
  sizes = c(5, 0.5),
  title.width = 1, # shrink the left title strip
  background.title = "white", # remove the grey background
  col.title = "black", # black title text
  fontcolor.title = "black", # black track-name font
  col.axis = "black", # black axis lines and tick marks
  fontcolor.axis = "black", # black axis tick labels
  # rotation.title = 0,       # horizontal track-title labels
  cex.axis = 1.2, # font size of y-axis tick labels
  cex.title = 1.2 # font size of y-axis track title
)

# volcano-style plot with EnhancedVolcano ####
# x axis: percRef, y axis: 0 - log10(pVal)
# EnhancedVolcano transforms y as -log10(y), so feed it 10^(log10(pVal)) to make
# the plotted height equal 0 - log10(pVal), as requested.
EnhancedVolcano_df <- ASoC_df

# colour points inside the [0.25, 0.75] percRef band grey; colour points
# outside the band (and significant) red, the rest black. colCustom is a named
# vector (one colour per point); the names become the legend labels.
volcano_colours <-
  ifelse(
    EnhancedVolcano_df$percRef >= 0.4 & EnhancedVolcano_df$percRef <= 0.6,
    "grey",
    ifelse(
      EnhancedVolcano_df$pVal < 0.05,
      "darkred",
      "grey20"
    )
  )
names(volcano_colours) <-
  ifelse(
    volcano_colours == "grey20",
    "0.4 <= percRef <= 0.6",
    ifelse(
      volcano_colours == "darkred",
      "outside band & p < 0.05",
      "outside band & NS"
    )
  )

EnhancedVolcano::EnhancedVolcano(
  EnhancedVolcano_df,
  lab = EnhancedVolcano_df$variantID,
  x = "percRef",
  y = "pVal",
  xlab = "percRef",
  ylab = "-log10(pVal)",
  xlim = c(0, 1),
  pointSize = ifelse(
    EnhancedVolcano_df$pVal < 0.05 &
      (EnhancedVolcano_df$percRef < 0.4 |
        EnhancedVolcano_df$percRef > 0.6),
    0.5,
    0.2
  ), # larger points for significant variants
  pCutoff = 0.05,
  FCcutoff = Inf, # disable the symmetric (around 0) FC cutoff lines
  vline = c(0.4, 0.5, 0.6), # custom vertical x-axis cutoff lines
  vlineCol = "black",
  vlineType = c("dashed", "solid", "dashed"),
  colCustom = volcano_colours,
  title = NULL,
  subtitle = NULL,
  caption = paste0("total = ", nrow(EnhancedVolcano_df), " SNPs"),
  legendPosition = "top",
  # legendJustification = "left",  # left-align the stacked legend under the title
  legendLabSize = 10,
  legendIconSize = 3
) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(ncol = 1) # keep items vertically stacked
  )
