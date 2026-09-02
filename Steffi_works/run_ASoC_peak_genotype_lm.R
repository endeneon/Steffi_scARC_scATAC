#! /usr/bin/env Rscript

# Per-SNP ASoC peak x genotype association in hepatocytes.
#
# For every SNP in hepatocyte_ASoC_genotypes_GT_summary_min3.tsv:
#   * find the PeakMatrix peak(s) it falls inside (SNPs with no peak are dropped)
#   * sum the peak's fragment counts per SAMPLE over the hepatocyte cells and
#     convert to deepTools-style RPGC
#   * regress the per-sample RPGC on the sample's genotype (additive dosage and
#     3-level factor) using the per-sample calls in
#     hepatocyte_ASoC_genotypes_GT_only.tsv
#   * write the R^2 / p-value table and one box-and-whisker panel per peak
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
snp_summary_file <- "hepatocyte_ASoC_genotypes_GT_summary_min3.tsv"
snp_genotype_file <- "hepatocyte_ASoC_genotypes_GT_only.tsv"
snp_annotation_file <-
  "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv"
archr_project_dir <- "ArchR_hepato"
writeout_dir <- "hepatocyte_ASoC_peak_genotype_lm"
plot_dir <- file.path(writeout_dir, "plots")

gt_levels <- c("0/0", "0/1", "1/1")
frag_len_cells <- 300L # cells sampled per arrow when estimating fragment length
frag_len_chrs <- c("chr1", "chr2")
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
# and that sample is dropped from the SNP's model and plot.
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

# Effective genome size for RPGC: ArchR chromSizes minus blacklisted bp.
.effective_genome_size <-
  function(ref_ArchR_obj) {
    ga <- ArchR::getGenomeAnnotation(ref_ArchR_obj)
    total <- sum(as.numeric(GenomicRanges::width(ga$chromSizes)))
    bl <- ga$blacklist
    if (!is.null(bl) && length(bl)) {
      total <- total -
        sum(as.numeric(GenomicRanges::width(GenomicRanges::reduce(bl))))
    }
    total
  }

# Mean fragment width for one Arrow file, from a random sample of its cells.
# Reading every fragment genome-wide is not tractable, but the mean stabilises
# quickly, so a few hundred cells on a couple of chromosomes is enough.
.estimate_frag_length_arrow <-
  function(arrow_file, cells, chrs, n_sample = 300L) {
    avail <- ArchR:::.availableCells(arrow_file, "PeakMatrix")
    cl <- intersect(cells, avail)
    if (!length(cl)) {
      return(NA_real_)
    }
    if (length(cl) > n_sample) {
      cl <- sample(cl, n_sample)
    }
    widths <- unlist(lapply(chrs, function(cc) {
      fr <- tryCatch(
        ArchR:::.getFragsFromArrow(
          arrow_file,
          chr = cc,
          out = "GRanges",
          cellNames = cl
        ),
        error = function(e) NULL
      )
      if (is.null(fr) || !length(fr)) {
        return(numeric(0))
      }
      as.numeric(GenomicRanges::width(fr))
    }))
    if (!length(widths)) {
      return(NA_real_)
    }
    mean(widths)
  }

# Overall model p-value from the F statistic; NULL for an intercept-only fit.
.model_p <-
  function(fit) {
    fs <- summary(fit)$fstatistic
    if (is.null(fs)) {
      return(NA_real_)
    }
    unname(stats::pf(fs[[1]], fs[[2]], fs[[3]], lower.tail = FALSE))
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

# unique peaks to pull out of the Arrow files
peak_rows <- sort(unique(df_pairs$peak_row))
feat_sel <- feat_all[peak_rows, c("seqnames", "idx", "start", "end")]
feat_sel_ids <- paste0(
  feat_sel$seqnames,
  ":",
  feat_sel$start,
  "-",
  feat_sel$end
)

# ---- 3. per-sample peak counts -> RPGC ---------------------------------------
egs <- .effective_genome_size(projHepatocytes)

cells_by_sample <- split(cell_meta$cellName, cell_meta$Sample)
nfrags_by_sample <- vapply(
  use_samples,
  function(s) sum(as.numeric(cell_meta$nFrags[cell_meta$Sample == s])),
  numeric(1)
)

n_par <- max(1L, min(workers_2_use, length(use_samples)))
par_cl <- parallel::makePSOCKcluster(n_par)
doParallel::registerDoParallel(par_cl)
on.exit(try(parallel::stopCluster(par_cl), silent = TRUE), add = TRUE)

# PSOCK workers are fresh R processes; foreach auto-exports the objects the loop
# body references but not the helper function, so name that one explicitly.
s <- NULL # foreach iterator; silences R CMD check
per_sample <-
  foreach(
    s = use_samples,
    .errorhandling = "pass",
    .export = ".estimate_frag_length_arrow",
    .packages = c("ArchR", "Matrix", "GenomicRanges")
  ) %dopar%
  {
    ArchR::addArchRThreads(threads = 1, force = TRUE)
    foreach::registerDoSEQ()
    set.seed(5813L)

    af <- arrow_files[[s]]
    cells_s <- intersect(
      cells_by_sample[[s]],
      ArchR:::.availableCells(af, "PeakMatrix")
    )
    if (!length(cells_s)) {
      return(list(sample = s, counts = NULL, frag_len = NA_real_))
    }

    m <- ArchR:::.getMatFromArrow(
      ArrowFile = af,
      featureDF = feat_sel,
      binarize = FALSE,
      useMatrix = "PeakMatrix",
      cellNames = cells_s
    )
    # .getMatFromArrow reorders rows back to the featureDF passed in and then
    # drops the rownames, so the result is positionally aligned with feat_sel.
    cnt <- Matrix::rowSums(m)
    stopifnot(length(cnt) == nrow(feat_sel))

    list(
      sample = s,
      counts = as.numeric(cnt),
      n_cells = length(cells_s),
      frag_len = .estimate_frag_length_arrow(
        af,
        cells_s,
        chrs = frag_len_chrs,
        n_sample = frag_len_cells
      )
    )
  }

bad <- vapply(per_sample, function(x) inherits(x, "error"), logical(1))
if (any(bad)) {
  stop(
    "Peak extraction failed for ",
    sum(bad),
    " sample(s): ",
    conditionMessage(per_sample[[which(bad)[1]]])
  )
}
names(per_sample) <- vapply(per_sample, `[[`, character(1), "sample")

mat_counts <- matrix(
  0,
  nrow = nrow(feat_sel),
  ncol = length(use_samples),
  dimnames = list(feat_sel_ids, use_samples)
)
frag_len_by_sample <- setNames(
  rep(NA_real_, length(use_samples)),
  use_samples
)
for (s in use_samples) {
  res <- per_sample[[s]]
  if (is.null(res$counts)) {
    next
  }
  mat_counts[, s] <- res$counts
  frag_len_by_sample[[s]] <- res$frag_len
}

# fall back to the cohort mean where a sample yielded no fragments to measure
frag_len_by_sample[!is.finite(frag_len_by_sample)] <-
  mean(frag_len_by_sample[is.finite(frag_len_by_sample)]) %||% 100

# deepTools RPGC ("1x genomic coverage"):
#   scale = effective_genome_size / (n_fragments * mean_fragment_length)
scale_by_sample <- egs / (nfrags_by_sample[use_samples] * frag_len_by_sample)
message(
  "RPGC scaling: effective genome size ",
  format(egs, big.mark = ","),
  " bp; mean fragment length ",
  round(mean(frag_len_by_sample), 1),
  " bp; scale factors ",
  paste(
    paste0(use_samples, "=", signif(scale_by_sample, 3)),
    collapse = ", "
  )
)

mat_rpgc <- sweep(mat_counts, 2, scale_by_sample, `*`)

df_rpgc_out <- data.frame(
  peak_id = rownames(mat_rpgc),
  mat_rpgc,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(
  df_rpgc_out,
  file.path(writeout_dir, "peak_RPGC_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---- 4. per-peak linear models ----------------------------------------------
peak_row_lookup <- setNames(seq_along(peak_rows), as.character(peak_rows))

fit_one <-
  function(i) {
    snp_i <- df_pairs$snp_row[[i]]
    peak_i <- peak_row_lookup[[as.character(df_pairs$peak_row[[i]])]]

    gt_row <- match(
      paste(df_pairs$CHROM[[i]], df_pairs$POS[[i]], sep = ":"),
      gt_key
    )
    gt_raw <- unlist(df_genotypes[gt_row, use_gt_cols[use_samples]])
    gt <- normalise_gt(gt_raw)

    keep <- !is.na(gt)
    d <- data.frame(
      sample = use_samples[keep],
      category = unname(sample_category[use_samples[keep]]),
      genotype = factor(gt[keep], levels = gt_levels),
      dosage = match(gt[keep], gt_levels) - 1L,
      rpgc = as.numeric(mat_rpgc[peak_i, use_samples[keep]]),
      stringsAsFactors = FALSE
    )
    d <- d[!is.na(d$genotype) & is.finite(d$rpgc), , drop = FALSE]

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
      beta_additive = NA_real_,
      r_squared_additive = NA_real_,
      p_value_additive = NA_real_,
      r_squared_factor = NA_real_,
      p_value_factor = NA_real_,
      stringsAsFactors = FALSE
    )

    usable <- nrow(d) >= min_samples_per_fit &&
      length(unique(d$dosage)) >= 2L &&
      stats::sd(d$rpgc) > 0
    if (!usable) {
      return(list(stats = res, data = d))
    }

    fit_add <- stats::lm(rpgc ~ dosage, data = d)
    sm_add <- summary(fit_add)
    res$beta_additive <- unname(stats::coef(fit_add)[["dosage"]])
    res$r_squared_additive <- sm_add$r.squared
    res$p_value_additive <- .model_p(fit_add)

    d_fac <- d
    d_fac$genotype <- droplevels(d_fac$genotype)
    fit_fac <- stats::lm(rpgc ~ genotype, data = d_fac)
    res$r_squared_factor <- summary(fit_fac)$r.squared
    res$p_value_factor <- .model_p(fit_fac)

    list(stats = res, data = d)
  }

fits <- lapply(seq_len(nrow(df_pairs)), fit_one)
df_results <- do.call(rbind, lapply(fits, `[[`, "stats"))
df_results$fdr_additive <- stats::p.adjust(
  df_results$p_value_additive,
  method = "BH"
)
df_results$fdr_factor <- stats::p.adjust(
  df_results$p_value_factor,
  method = "BH"
)
df_results <- df_results[
  order(df_results$p_value_additive, na.last = TRUE),
  ,
  drop = FALSE
]
rownames(df_results) <- NULL

write.table(
  df_results,
  file.path(writeout_dir, "ASoC_peak_genotype_lm_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
qs2::qs_save(
  df_results,
  file.path(writeout_dir, "ASoC_peak_genotype_lm_results.qs2"),
  nthreads = min(4L, workers_2_use)
)
print(paste0(
  "Wrote ",
  nrow(df_results),
  " SNP-peak model result(s) to ",
  writeout_dir
))

# ---- 5. per-peak box-and-whisker panels --------------------------------------
plot_data <- lapply(fits, `[[`, "data")
plot_stats <- lapply(fits, `[[`, "stats")
keep_panel <- vapply(plot_data, function(d) nrow(d) > 0L, logical(1))
plot_data <- plot_data[keep_panel]
plot_stats <- plot_stats[keep_panel]

# "Intron (ENST00000537821.2/2212, intron 3 of 4)" overruns a 1/6-page panel;
# the transcript detail is dropped for the label only, not for the table.
short_annotation <-
  function(x) sub(" \\(ENST[^)]*\\)$", "", x)

panel_titles <- vapply(
  plot_stats,
  function(st) {
    sprintf(
      "%s (%s)\n%s:%s %s>%s | %s\npeak %s\nR2=%s p=%s",
      st$variantID %NA% "NA",
      st$SYMBOL %NA% "NA",
      st$CHROM,
      st$POS,
      st$REF,
      st$ALT,
      short_annotation(st$annotation %NA% "NA"),
      st$peak_id,
      ifelse(
        is.na(st$r_squared_additive),
        "NA",
        formatC(st$r_squared_additive, format = "f", digits = 3)
      ),
      ifelse(
        is.na(st$p_value_additive),
        "NA",
        formatC(st$p_value_additive, format = "g", digits = 3)
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
      mu <- stats::aggregate(rpgc ~ genotype, data = d, FUN = mean)

      ggplot2::ggplot(d, ggplot2::aes(x = genotype, y = rpgc)) +
        ggplot2::geom_boxplot(
          ggplot2::aes(fill = genotype),
          outlier.shape = NA,
          width = 0.65,
          linewidth = 0.25
        ) +
        ggplot2::geom_point(
          ggplot2::aes(shape = category),
          position = ggplot2::position_jitter(width = 0.15, height = 0),
          size = 0.7,
          colour = "grey20"
        ) +
        ggplot2::geom_point(
          data = mu,
          ggplot2::aes(x = genotype, y = rpgc),
          shape = 4,
          size = 1.8,
          stroke = 0.7,
          colour = "darkred",
          inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_brewer(palette = "Set2", drop = FALSE) +
        ggplot2::scale_shape_manual(
          values = c(Primary = 16, Resistant = 17),
          drop = FALSE
        ) +
        ggplot2::labs(
          title = panel_titles[[k]],
          x = NULL,
          y = "RPGC"
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
      sprintf("ASoC_peak_genotype_boxplots_page_%02d.pdf", pg)
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
