#! /usr/bin/env Rscript

# Per-SNP ASoC peak x genotype x sample-category composition in hepatocytes.
#
# A variant of run_ASoC_peak_genotype_lm.R / run_ASoC_peak_category_lm.R. For
# every SNP in hepatocyte_ASoC_genotypes_GT_summary_min3.tsv:
#   * find the PeakMatrix peak(s) it falls inside (SNPs with no peak are dropped)
#   * assign each SAMPLE a genotype (0/0, 0/1, 1/1) from the per-sample calls in
#     hepatocyte_ASoC_genotypes_GT_only.tsv and a category
#     (projHepatocytes@cellColData$category: Primary or Resistant)
#   * draw a dodged column graph: x = genotype, and within each genotype the
#     number of samples that are Primary vs Resistant (Dark2 fill)
#   * fit a binomial additive model of category on allele dosage
#     (mgcv::gam(category ~ dosage, family = binomial)) and print the effect
#     size (log-odds per alt allele) and its Wald p-value in the panel title
#
# There is no peak accessibility (RPGC) in the model here: the outcome is the
# sample's Primary/Resistant category, the predictor is genotype dosage. No mean
# is computed (the bars are integer sample counts).
#
# init
{
  library(ArchR)
  library(future)
  library(stringr)

  library(parallel)
  library(foreach)
  library(doParallel)

  library(Matrix)

  library(qs2)

  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)

  library(mgcv)

  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
}

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)

# LSF does not expose its core allocation to future::availableCores(), so read
# LSB_DJOB_NUMPROC directly and fall back to parallelly's detection.
lsf_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = NA))
available_cores <-
  if (!is.na(lsf_cores) && lsf_cores >= 1) {
    lsf_cores
  } else {
    parallelly::availableCores()
  }
workers_2_use <- max(1L, min(as.integer(available_cores) - 1L, 30L))

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = FALSE)
  options(useUCSCChromosomeNames = FALSE)
  set.seed(42)
}
addArchRThreads(threads = 1)
addArchRGenome("hg38")

print(paste0(
  "Detected ",
  available_cores,
  " cores; using ",
  workers_2_use,
  " parallel workers."
))

# ---- parameters -------------------------------------------------------------
# snp_summary_file <- "hepatocyte_ASoC_genotypes_GT_summary_min3.tsv"
# snp_genotype_file <- "hepatocyte_ASoC_genotypes_GT_only.tsv"
# snp_annotation_file <-
#   "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv"
# archr_project_dir <- "ArchR_hepato"

snp_summary_file <- "test_ASoC_w_WASP/ASoC_Macrophage_genotyping_output/Macrophage_ASoC_genotypes_GT_summary_min3.tsv"
snp_genotype_file <- "test_ASoC_w_WASP/ASoC_Macrophage_genotyping_output/Macrophage_ASoC_genotypes_GT_only.tsv"
snp_annotation_file <-
  "sig_ASoC_by_celltype/sig_ASoC_in_Macrophage_annotated.tsv"
archr_project_dir <- "ArchR_macrophages"

# Fail-fast: every required input must exist before we start. Report ALL
# missing paths at once (files checked as files, the ArchR project as a dir),
# resolved to absolute paths so it's obvious WHERE we looked, rather than dying
# on the first one.
local({
  required_files <- c(snp_summary_file, snp_genotype_file, snp_annotation_file)
  missing_files <- required_files[!file.exists(required_files)]
  missing_dirs <- archr_project_dir[!dir.exists(archr_project_dir)]
  if (length(missing_files) || length(missing_dirs)) {
    abs <- function(p) {
      # normalizePath() leaves nonexistent relative paths relative, so anchor
      # them to the working dir ourselves before normalising.
      p <- ifelse(
        grepl("^(/|[A-Za-z]:)", p),
        p,
        file.path(getwd(), p)
      )
      normalizePath(p, mustWork = FALSE)
    }
    stop(
      "Missing required input(s) (working dir: ",
      getwd(),
      "):\n",
      paste0("  [file] ", abs(missing_files), collapse = "\n"),
      if (length(missing_files) && length(missing_dirs)) "\n",
      paste0("  [dir]  ", abs(missing_dirs), collapse = "\n"),
      call. = FALSE
    )
  }
})
writeout_dir <- "macrophage_ASoC_peak_genotype_category_gam"
if (!dir.exists(writeout_dir)) {
  dir.create(writeout_dir, recursive = TRUE)
}
plot_dir <- file.path(writeout_dir, "plots")

gt_levels <- c("0/0", "0/1", "1/1")
# "Resistant" is the reference (0) level; "Primary" is the modelled event (1).
category_levels <- c("Primary", "Resistant")
min_samples_per_fit <- 6L

panels_per_row <- 6L
panels_per_col <- 6L
panels_per_page <- panels_per_row * panels_per_col
page_width_in <- 11
page_height_in <- 8.5

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x
# also catches NA lookups from named-vector indexing
`%NA%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# ---- helpers ----------------------------------------------------------------

# Unphased calls: "1/0" and "0/1" are the same level. Anything that is not a
# biallelic 0/1 call (missing "./.", multi-allelic "0/2", "1/2", ...) becomes NA
# and that sample is dropped from the SNP's plot and model.
normalise_gt <-
  function(gt) {
    vapply(
      strsplit(as.character(gt), "[/|]"),
      function(alleles) {
        if (length(alleles) != 2L || !all(alleles %in% c("0", "1"))) {
          return(NA_character_)
        }
        paste(sort(as.integer(alleles)), collapse = "/")
      },
      character(1)
    )
  }

# ---- 1. SNP list + per-sample genotypes -------------------------------------
df_sig_snp_list <-
  read.table(
    snp_summary_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
print(paste0("Read ", nrow(df_sig_snp_list), " SNPs from ", snp_summary_file))

df_genotypes <-
  read.table(
    snp_genotype_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )

gt_cols <- grep("_GT$", colnames(df_genotypes), value = TRUE)
# genotype header "multiome_11t_GT" / "atac_pt11_pre_GT" -> ArchR sample "X__11t"
gt_sample_key <- sub("^(atac|multiome)_", "", sub("_GT$", "", gt_cols))
names(gt_cols) <- gt_sample_key

gt_key <- paste(df_genotypes$CHROM, df_genotypes$POS, sep = ":")
snp_key <- paste(df_sig_snp_list$CHROM, df_sig_snp_list$POS, sep = ":")
missing_gt <- setdiff(snp_key, gt_key)
if (length(missing_gt)) {
  stop(
    length(missing_gt),
    " SNP(s) from the summary table have no per-sample genotypes, e.g. ",
    paste(head(missing_gt, 3), collapse = ", ")
  )
}

# rsID / gene / peak-annotation lookup. A locus can appear more than once in the
# annotated table (e.g. one row per alt allele), so distinct values are collapsed
# with ";" rather than silently taking the first.
df_snp_annot <-
  read.table(
    snp_annotation_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
annot_key <- paste(df_snp_annot$seqnames, df_snp_annot$start, sep = ":")

collapse_annot <-
  function(column) {
    vals <- split(as.character(df_snp_annot[[column]]), annot_key)
    vapply(
      vals,
      function(x) {
        x <- unique(x[!is.na(x) & nzchar(x)])
        if (!length(x)) NA_character_ else paste(x, collapse = ";")
      },
      character(1)
    )
  }

annot_lookup <- list(
  variantID = collapse_annot("variantID"),
  SYMBOL = collapse_annot("SYMBOL"),
  annotation = collapse_annot("annotation")
)

# ---- 2. project, peaks and SNP -> peak overlaps ------------------------------
# The peak accessibility matrix is not needed here (category, not RPGC, is the
# outcome), but the SNP -> peak overlap is still used so the panel layout and
# per-pair looping match the sibling scripts.
projHepatocytes <- ArchR::loadArchRProject(path = archr_project_dir)

cell_meta <- as.data.frame(
  projHepatocytes@cellColData[, c("Sample", "nFrags", "category")]
)
cell_meta$cellName <- rownames(cell_meta)
cell_meta$Sample <- as.character(cell_meta$Sample)
cell_meta$category <- as.character(cell_meta$category)

arrow_files <- ArchR::getArrowFiles(projHepatocytes)
if (!length(arrow_files)) {
  stop("No Arrow files found in ", archr_project_dir)
}
# arrow file basename "X__11t.arrow" is the ArchR sample name
names(arrow_files) <- sub("\\.arrow$", "", basename(arrow_files))

archr_samples <- sort(unique(cell_meta$Sample))
sample_key <- sub("^X__", "", archr_samples)
names(archr_samples) <- sample_key

shared_key <- intersect(sample_key, gt_sample_key)
if (!length(shared_key)) {
  stop("No ArchR sample name could be matched to a genotype column.")
}
if (length(shared_key) < length(archr_samples)) {
  message(
    length(archr_samples) - length(shared_key),
    " ArchR sample(s) have no genotype column and were dropped: ",
    paste(setdiff(sample_key, shared_key), collapse = ", ")
  )
}
use_samples <- unname(archr_samples[shared_key])
use_gt_cols <- unname(gt_cols[shared_key])
names(use_gt_cols) <- use_samples

sample_category <- vapply(
  use_samples,
  function(s) unique(cell_meta$category[cell_meta$Sample == s])[[1]],
  character(1)
)
if (!all(sample_category %in% category_levels)) {
  stop(
    "Unexpected category value(s): ",
    paste(setdiff(unique(sample_category), category_levels), collapse = ", ")
  )
}
print(table(sample_category))

# The PeakMatrix FeatureDF is the authoritative row index of the peak matrix,
# so overlaps are computed against it rather than against getPeakSet().
feat_all <- ArchR:::.getFeatureDF(unname(arrow_files), "PeakMatrix")
feat_all <- as.data.frame(feat_all)
feat_all$seqnames <- as.character(feat_all$seqnames)

peak_gr <- GenomicRanges::GRanges(
  seqnames = feat_all$seqnames,
  ranges = IRanges::IRanges(start = feat_all$start, end = feat_all$end)
)

snp_gr <- GenomicRanges::GRanges(
  seqnames = df_sig_snp_list$CHROM,
  ranges = IRanges::IRanges(
    start = df_sig_snp_list$POS,
    width = 1L
  )
)

hits <- GenomicRanges::findOverlaps(snp_gr, peak_gr, ignore.strand = TRUE)
if (!length(hits)) {
  stop("None of the SNPs overlap a called peak.")
}

df_pairs <- data.frame(
  snp_row = S4Vectors::queryHits(hits),
  peak_row = S4Vectors::subjectHits(hits),
  stringsAsFactors = FALSE
)
df_pairs <- cbind(
  df_sig_snp_list[df_pairs$snp_row, c("CHROM", "POS", "REF", "ALT")],
  df_pairs[, c("snp_row", "peak_row")]
)
pair_key <- paste(df_pairs$CHROM, df_pairs$POS, sep = ":")
df_pairs$variantID <- unname(annot_lookup$variantID[pair_key])
df_pairs$SYMBOL <- unname(annot_lookup$SYMBOL[pair_key])
df_pairs$annotation <- unname(annot_lookup$annotation[pair_key])
if (anyNA(df_pairs$variantID)) {
  message(
    sum(is.na(df_pairs$variantID)),
    " SNP-peak pair(s) had no match in ",
    snp_annotation_file,
    "; variantID/SYMBOL/annotation left as NA."
  )
}
df_pairs$peak_chr <- feat_all$seqnames[df_pairs$peak_row]
df_pairs$peak_start <- feat_all$start[df_pairs$peak_row]
df_pairs$peak_end <- feat_all$end[df_pairs$peak_row]
df_pairs$peak_id <- paste0(
  df_pairs$peak_chr,
  ":",
  df_pairs$peak_start,
  "-",
  df_pairs$peak_end
)
rownames(df_pairs) <- NULL

print(paste0(
  nrow(df_pairs),
  " SNP-peak pair(s) from ",
  length(unique(df_pairs$snp_row)),
  " SNP(s); ",
  nrow(df_sig_snp_list) - length(unique(df_pairs$snp_row)),
  " SNP(s) had no overlapping peak and were dropped."
))

# ---- 3. per-pair genotype x category composition + binomial GAM -------------
# For each SNP-peak pair: build the per-sample (genotype, category) table, fit
# gam(category ~ dosage, family = binomial). Dosage takes only 0/1/2, so a
# smooth s(dosage) is unidentifiable; a linear dosage term gives an
# interpretable log-odds-per-allele effect size with a Wald p-value.
prep_one <-
  function(i) {
    gt_row <- match(
      paste(df_pairs$CHROM[[i]], df_pairs$POS[[i]], sep = ":"),
      gt_key
    )
    gt_raw <- unlist(df_genotypes[gt_row, use_gt_cols[use_samples]])
    gt <- normalise_gt(gt_raw)

    keep <- !is.na(gt)
    d <- data.frame(
      sample = use_samples[keep],
      genotype = factor(gt[keep], levels = gt_levels),
      dosage = match(gt[keep], gt_levels) - 1L,
      category = factor(
        unname(sample_category[use_samples[keep]]),
        levels = category_levels
      ),
      stringsAsFactors = FALSE
    )
    d <- d[!is.na(d$genotype) & !is.na(d$category), , drop = FALSE]

    res <- data.frame(
      CHROM = df_pairs$CHROM[[i]],
      POS = df_pairs$POS[[i]],
      REF = df_pairs$REF[[i]],
      ALT = df_pairs$ALT[[i]],
      variantID = df_pairs$variantID[[i]],
      SYMBOL = df_pairs$SYMBOL[[i]],
      annotation = df_pairs$annotation[[i]],
      peak_chr = df_pairs$peak_chr[[i]],
      peak_start = df_pairs$peak_start[[i]],
      peak_end = df_pairs$peak_end[[i]],
      peak_id = df_pairs$peak_id[[i]],
      n_samples = nrow(d),
      n_00 = sum(d$genotype == "0/0"),
      n_01 = sum(d$genotype == "0/1"),
      n_11 = sum(d$genotype == "1/1"),
      n_primary = sum(d$category == "Primary"),
      n_resistant = sum(d$category == "Resistant"),
      effect_logodds = NA_real_,
      effect_se = NA_real_,
      p_value_gam = NA_real_,
      stringsAsFactors = FALSE
    )

    # Need both categories present and >=2 dosage levels for the effect to be
    # estimable, plus a minimum sample size.
    usable <- nrow(d) >= min_samples_per_fit &&
      length(unique(d$dosage)) >= 2L &&
      length(unique(d$category)) == 2L
    if (!usable) {
      return(list(stats = res, data = d))
    }

    # category as 0/1: Resistant (reference) = 0, Primary (event) = 1.
    d_fit <- d
    d_fit$y <- as.integer(d_fit$category == "Primary")

    fit <- tryCatch(
      mgcv::gam(y ~ dosage, family = stats::binomial(), data = d_fit),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      sm <- summary(fit)
      pt <- sm$p.table
      if ("dosage" %in% rownames(pt)) {
        res$effect_logodds <- unname(pt["dosage", "Estimate"])
        res$effect_se <- unname(pt["dosage", "Std. Error"])
        res$p_value_gam <- unname(pt["dosage", ncol(pt)])
      }
    }

    list(stats = res, data = d)
  }

fits <- lapply(seq_len(nrow(df_pairs)), prep_one)
df_results <- do.call(rbind, lapply(fits, `[[`, "stats"))
df_results$fdr_gam <- stats::p.adjust(df_results$p_value_gam, method = "BH")
df_results <- df_results[
  order(df_results$p_value_gam, na.last = TRUE),
  ,
  drop = FALSE
]
rownames(df_results) <- NULL

write.table(
  df_results,
  file.path(writeout_dir, "ASoC_peak_genotype_category_gam_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
qs2::qs_save(
  df_results,
  file.path(writeout_dir, "ASoC_peak_genotype_category_gam_results.qs2"),
  nthreads = min(4L, workers_2_use)
)
print(paste0(
  "Wrote ",
  nrow(df_results),
  " SNP-peak model result(s) to ",
  writeout_dir
))

# ---- 4. per-peak genotype x category count column panels ---------------------
plot_data <- lapply(fits, `[[`, "data")
plot_stats <- lapply(fits, `[[`, "stats")
keep_panel <- vapply(plot_data, function(d) nrow(d) > 0L, logical(1))
plot_data <- plot_data[keep_panel]
plot_stats <- plot_stats[keep_panel]

# "Intron (ENST00000537821.2/2212, intron 3 of 4)" overruns a 1/6-page panel;
# the transcript detail is dropped for the label only, not for the table.
short_annotation <-
  function(x) sub(" \\(ENST[^)]*\\)$", "", x)

# GAM effect size (log-odds per alt allele) and Wald p in place of R2/p.
panel_titles <- vapply(
  plot_stats,
  function(st) {
    sprintf(
      "%s (%s)\n%s:%s %s>%s | %s\npeak %s\nbeta=%s p=%s",
      st$variantID %NA% "NA",
      st$SYMBOL %NA% "NA",
      st$CHROM,
      st$POS,
      st$REF,
      st$ALT,
      short_annotation(st$annotation %NA% "NA"),
      st$peak_id,
      ifelse(
        is.na(st$effect_logodds),
        "NA",
        formatC(st$effect_logodds, format = "f", digits = 3)
      ),
      ifelse(
        is.na(st$p_value_gam),
        "NA",
        formatC(st$p_value_gam, format = "g", digits = 3)
      )
    )
  },
  character(1)
)

n_pages <- max(1L, ceiling(length(plot_data) / panels_per_page))
page_index <- lapply(seq_len(n_pages), function(pg) {
  seq(
    (pg - 1L) * panels_per_page + 1L,
    min(pg * panels_per_page, length(plot_data))
  )
})

n_par <- max(1L, min(workers_2_use, n_pages))
par_cl <- parallel::makePSOCKcluster(n_par)
doParallel::registerDoParallel(par_cl)
on.exit(try(parallel::stopCluster(par_cl), silent = TRUE), add = TRUE)

pg <- NULL # foreach iterator
page_files <-
  foreach(
    pg = seq_len(n_pages),
    .combine = c,
    .errorhandling = "pass",
    .packages = c("ggplot2", "cowplot", "RColorBrewer")
  ) %dopar%
  {
    idx <- page_index[[pg]]
    panels <- lapply(idx, function(k) {
      d <- plot_data[[k]]
      d$genotype <- droplevels(d$genotype)
      # sample counts per genotype x category (all levels kept so empty
      # category bars still occupy a slot within each genotype)
      counts <- as.data.frame(
        table(genotype = d$genotype, category = d$category)
      )
      counts$category <- factor(counts$category, levels = category_levels)

      ggplot2::ggplot(
        counts,
        ggplot2::aes(x = genotype, y = Freq, fill = category)
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(preserve = "single"),
          width = 0.7,
          linewidth = 0.2,
          colour = "grey20"
        ) +
        ggplot2::scale_fill_brewer(palette = "Dark2", drop = FALSE) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.1))
        ) +
        ggplot2::labs(
          title = panel_titles[[k]],
          x = NULL,
          y = "sample count"
        ) +
        ggplot2::theme_bw(base_size = 5) +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(size = 3.4, lineheight = 1.05),
          axis.text = ggplot2::element_text(size = 4),
          axis.title.y = ggplot2::element_text(size = 4),
          panel.grid.minor = ggplot2::element_blank()
        )
    })

    grid_pg <- cowplot::plot_grid(
      plotlist = panels,
      ncol = panels_per_row,
      nrow = panels_per_col
    )
    out_file <- file.path(
      plot_dir,
      sprintf("ASoC_peak_genotype_category_columns_page_%02d.pdf", pg)
    )
    ggplot2::ggsave(
      out_file,
      plot = grid_pg,
      width = page_width_in,
      height = page_height_in,
      units = "in",
      device = "pdf"
    )
    out_file
  }

parallel::stopCluster(par_cl)

print(paste0(
  "Wrote ",
  length(plot_data),
  " panel(s) across ",
  n_pages,
  " PDF page(s) in ",
  plot_dir
))
print("Done.")
