# Chat Session Summary

## Session Metadata

- **Date**: 2026-09-02 (work), summarized 2026-09-04
- **Workspace**: `/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main`
- **Working directory for all scripts**: `Multiome_main/Steffi_works`
- **Repository**: `endeneon/Steffi_scARC_scATAC`, branch `main`
- **Compute environment**: St. Jude HPC, LSF scheduler, conda env
  `/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR`
  (R 4.3.0, ggplot2 3.5.2, cowplot 1.2.0)
- **Goal**: Build a new R pipeline that tests whether chromatin accessibility at
  ASoC-SNP-overlapping ATAC peaks in hepatocytes depends on donor genotype, then
  visualize each peak and submit the job to LSF.

## Tasks Completed

1. **Clarified ambiguous requirements before writing code.**
   The user's original request asked for a per-sample value of
   `Primary / (Primary + Resistant)`. Inspection of the ArchR project showed
   every sample belongs to exactly one `category` (14 Primary, 16 Resistant), so
   that ratio is undefined per sample. Six clarifying questions were asked; the
   answers fixed the design (see *Key Decisions & Rationale*).
   **Outcome**: unambiguous spec before implementation.

2. **Explored the ArchR project and input files to ground the design.**
   - `ArchR_hepato/ArrowFiles/` contains 30 Arrow files named `X__11t.arrow`,
     `X__pt11_pre.arrow`, etc.
   - `h5ls` confirmed each Arrow file stores a `PeakMatrix` group, and
     `PeakMatrix/Info/FeatureDF` is an HDF5 compound dataset with fields
     `seqnames`, `idx`, `start`, `end` (271,608 peaks).
   - `cellColData` has 291,150 cells and already carries the `category` column
     (`Primary` for `*t` samples, `Resistant` for `*pre` samples).
   - Sample-name mapping verified 1:1: ArchR `X__11t` ↔ genotype column
     `multiome_11t_GT`; ArchR `X__pt11_pre` ↔ `atac_pt11_pre_GT`.
   - All required R packages (`ggplot2`, `cowplot`, `qs2`, `ArchR`, `foreach`,
     `doParallel`, `GenomicRanges`, `RColorBrewer`) were already installed, so
     **no package installation was needed**.
   **Outcome**: no guesswork in the implementation; no `mamba install` required.

3. **Created the analysis script `Steffi_works/run_ASoC_peak_genotype_lm.R`.**
   Implements the five requested steps plus parallelism. Details in
   *Code Changes*.
   **Outcome**: runs end-to-end; 96 of 97 SNPs overlap a peak (1 dropped),
   producing 96 SNP–peak models and 3 PDF pages.

4. **Created the LSF wrapper `Steffi_works/bsub_run_ASoC_peak_genotype_lm.sh`.**
   Modeled on the existing `bsub_run_generate_atac.sh`.
   **Outcome**: submittable job script; `bash -n` syntax check passes.

5. **Diagnosed and fixed an all-zero peak-count bug.**
   The first run produced a results table where every `r_squared`/`p_value` was
   `NA` and every RPGC value was `0`. Root cause: `ArchR:::.getMatFromArrow`
   assigns per-chromosome rownames `f1..fN` internally, but its final two lines
   are `mat <- mat[rownames(featureDF), , drop = FALSE]` followed by
   `rownames(mat) <- NULL`. The initial implementation indexed the row sums by
   name (`Matrix::rowSums(m)[paste0("f", seq_len(nrow(feat_sel)))]`), which
   returned all `NA`, and the `cnt[is.na(cnt)] <- 0` guard silently converted
   them to zeros.
   **Fix**: index positionally, because the function already guarantees the row
   order matches the supplied `featureDF`; added a `stopifnot` length guard.
   **Outcome**: 0 of 96 rows `NA`; top hit `chr19:39194052` (rs11668945,
   `NCCRP1`) with `R² = 0.375`, `p = 3.25e-4`.

6. **Silenced redundant `foreach` export warnings.**
   `foreach` auto-detects and exports globals referenced in the loop body, so
   listing them again in `.export` produced
   `already exporting variable(s): ...` warnings. Trimmed `.export` to the one
   symbol *not* auto-detected (`.estimate_frag_length_arrow`, a function) and
   removed `.export` entirely from the plotting loop.
   **Outcome**: clean run log.

7. **Confirmed the `mkdir -p main_log` placement in the LSF wrapper.**
   The user stated their LSF accepts shell commands before `#BSUB` directives,
   so `mkdir -p main_log` was moved back above the directive block to match the
   house style of `bsub_run_generate_atac.sh` (and to guarantee the `-o`/`-e`
   log directory exists at submission time).
   **Outcome**: wrapper matches existing project convention.

8. **Added `variantID`, `SYMBOL` and `annotation` columns.**
   Joined from `sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv` on
   `seqnames:start` → `CHROM:POS`. All 96 SNP–peak pairs matched. The three
   columns appear immediately after `ALT` in the output table and on every plot
   panel title.
   **Outcome**: results are human-interpretable without a separate lookup.

9. **Shortened long annotation strings in plot labels only.**
   Strings such as `Intron (ENST00000537821.2/2212, intron 3 of 4)` overran a
   1/6-page panel and bled into the neighbouring panel. The plot label strips
   the trailing `(ENST...)` detail; the table keeps the full string.
   **Outcome**: verified by rendering page 1 with `pdftoppm` and viewing it.

## Key Decisions & Rationale

| Question                                         | Decision                                                                                           | Rationale                                                                                                                                         |
| ------------------------------------------------ | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| Regression response variable                     | **Per-sample peak RPGC**; `category` kept only as an annotation (plot point shape)                 | Each sample is either Primary or Resistant, so a per-sample `Primary/(Primary+Resistant)` ratio is mathematically undefined. User chose option A. |
| Genotype coding                                  | **Both** additive dosage (0/1/2) *and* 3-level factor                                              | User requested both; additive gives a directional β, factor allows non-additive effects.                                                          |
| Non-biallelic / missing GT (`./.`, `0/2`, `1/2`) | **Drop that sample** from the SNP's model and plot                                                 | Avoids inventing an allele mapping. Affects e.g. `chr6:29927225` (ALT `C,A`).                                                                     |
| Unphased genotypes                               | **`1/0` normalised to `0/1`** — alleles sorted before pasting                                      | Explicit user follow-up: the data are unphased, so heterozygotes are one level.                                                                   |
| RPGC formula                                     | `counts × EGS / (nFrags_sample × mean_frag_len_sample)`; **no** division by peak width             | User chose to match the existing `plot_gviz_pileups_by_category.R` convention rather than convert to per-bp coverage.                             |
| Effective genome size                            | ArchR `chromSizes` total minus reduced blacklist = **2,860,760,506 bp**                            | Same helper as the existing Gviz script.                                                                                                          |
| Mean fragment length                             | Estimated **per sample** from ≤300 random cells on `chr1`+`chr2`                                   | Genome-wide fragment reads are intractable; the mean stabilises quickly. Observed ≈150.1 bp.                                                      |
| Plot geometry                                    | **`geom_boxplot(outlier.shape = NA)`** + jittered `geom_point` + `shape = 4` `darkred` mean marker | User confirmed "bar+whisker" meant box-and-whisker with outliers hidden.                                                                          |
| Genotype-level columns                           | Columns **5–7** (`n_0/0`, `n_0/1`, `n_1/1`)                                                        | The user wrote "column 4–7"; column 4 is `ALT` and column 8 is `n_other`, so 5–7 are the three genotype levels. Confirmed.                        |
| Multiple peaks per SNP                           | **Keep all** overlapping peaks, one result row each                                                | User choice; no arbitrary tie-breaking.                                                                                                           |
| Peak coordinate source                           | `ArchR:::.getFeatureDF(arrowFiles, "PeakMatrix")`, **not** `getPeakSet()`                          | The FeatureDF *is* the row index of the PeakMatrix, so overlaps computed against it cannot desynchronise from the extracted matrix rows.          |
| Parallel backend                                 | **PSOCK** + `doParallel`                                                                           | Arrow/HDF5 handles cannot be shared across forks; PSOCK workers open their own. Matches the existing project pattern.                             |
| PDF device                                       | Base `"pdf"` rather than `cairo_pdf`                                                               | Avoids a hard cairo dependency; one page per file so multi-page support is irrelevant.                                                            |
| BLAS/OMP threads in wrapper                      | Pinned to 1                                                                                        | `foreach` owns the parallelism; prevents oversubscription of the 32 LSF slots.                                                                    |

## Code Changes

### 1. `Steffi_works/run_ASoC_peak_genotype_lm.R` (created)

New ~700-line analysis script. Structure:

- **init** — package loads, `setwd()`, LSF core detection via `LSB_DJOB_NUMPROC`,
  `addArchRThreads(1)`, `addArchRGenome("hg38")`.
- **parameters** — input paths, `gt_levels`, panel layout, output directories.
- **helpers** — `normalise_gt`, `.effective_genome_size`,
  `.estimate_frag_length_arrow`, `.model_p`, `%||%`, `%NA%`.
- **step 1** — read the SNP summary, per-sample genotypes and the annotation table.
- **step 2** — load the ArchR project, map sample names, build peak/SNP `GRanges`, `findOverlaps`.
- **step 3** — parallel per-sample PeakMatrix extraction → RPGC matrix.
- **step 4** — per-peak `lm` (additive + factor), BH FDR, TSV + qs2 output.
- **step 5** — parallel per-page ggplot rendering, one PDF per page.

Key excerpts:

```r
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
```

```r
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
```

```r
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
```

Annotation join (rsID / gene symbol / genomic annotation):

```r
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
```

Parallel per-sample PeakMatrix extraction (contains the bug fix):

```r
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
```

RPGC scaling:

```r
# fall back to the cohort mean where a sample yielded no fragments to measure
frag_len_by_sample[!is.finite(frag_len_by_sample)] <-
  mean(frag_len_by_sample[is.finite(frag_len_by_sample)]) %||% 100

# deepTools RPGC ("1x genomic coverage"):
#   scale = effective_genome_size / (n_fragments * mean_fragment_length)
scale_by_sample <- egs / (nfrags_by_sample[use_samples] * frag_len_by_sample)

mat_rpgc <- sweep(mat_counts, 2, scale_by_sample, `*`)
```

Per-peak linear models:

```r
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
```

Panel title construction:

```r
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
```

Per-panel plot specification:

```r
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
```

### 2. `Steffi_works/bsub_run_ASoC_peak_genotype_lm.sh` (created)

```bash
#! /bin/bash

mkdir -p main_log
#BSUB -n 32
#BSUB -R "rusage[mem=20G]"
# Keep all slots on ONE host so the PSOCK workers and their ArchR/HDF5 reads
# stay local; without this LSF may spread -n 32 across nodes.
#BSUB -R "span[hosts=1]"
#BSUB -q "standard"
#BSUB -J "run_ASoC_peak_genotype_lm"
#BSUB -o main_log/ASoC_peak_genotype_lm_%J.out
#BSUB -e main_log/ASoC_peak_genotype_lm_%J.err

# Regress per-sample peak RPGC on ASoC SNP genotype and draw the per-peak
# box-and-whisker panels. See run_ASoC_peak_genotype_lm.R.

# set error trap
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# Report the failing command on non-zero exit only.
_on_exit() {
	local ec=$?
	[[ ${ec} -ne 0 ]] && echo "\"${last_command}\" command failed with exit code ${ec}." >&2
}
trap '_on_exit' EXIT
shopt -s nullglob
shopt -s extglob

module load conda3/202402
conda activate /research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/Anaconda/miniconda3/envs/jwen_scRNA_singCellaR

# 0. set up base dir
{
	base_dir="/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
	if [[ ! -d "${base_dir}" ]]; then
		echo "ERROR: base directory not found: ${base_dir}" >&2
		exit 1
	fi
	cd "${base_dir}"
	mkdir -p "${base_dir}/main_log"
}

# 1. set up the environment
# HDF5/BLAS are called from inside each PSOCK worker, so keep them single
# threaded and let foreach own the parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript run_ASoC_peak_genotype_lm.R

set +e
```

## Inputs and Outputs

### Inputs (all relative to `Steffi_works/`)

| Path                                                        | Role                                                                                                                            |
| ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `hepatocyte_ASoC_genotypes_GT_summary_min3.tsv`             | 97 SNPs; source of `df_sig_snp_list`. Columns: `CHROM POS REF ALT n_0/0 n_0/1 n_1/1 n_other`.                                   |
| `hepatocyte_ASoC_genotypes_GT_only.tsv`                     | Per-sample genotype calls; 30 `*_GT` columns.                                                                                   |
| `sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv` | Source of `variantID`, `SYMBOL`, `annotation`; joined on `seqnames:start`. Must be read with `quote = ""`, `comment.char = ""`. |
| `ArchR_hepato/`                                             | ArchR project (30 Arrow files, 291,150 hepatocyte cells, `PeakMatrix` with 271,608 peaks).                                      |

### Outputs — directory `hepatocyte_ASoC_peak_genotype_lm/`

| File                                                | Contents                                                      |
| --------------------------------------------------- | ------------------------------------------------------------- |
| `ASoC_peak_genotype_lm_results.tsv`                 | 96 rows, 22 columns (below). Sorted by `p_value_additive`.    |
| `ASoC_peak_genotype_lm_results.qs2`                 | Same data frame via `qs2::qs_save`.                           |
| `peak_RPGC_by_sample.tsv`                           | Peak × sample RPGC matrix (`peak_id` + 30 sample columns).    |
| `plots/ASoC_peak_genotype_boxplots_page_01..03.pdf` | 6×6 panels per page, 11 × 8.5 in landscape, one PDF per page. |

Result table column order:

```
CHROM, POS, REF, ALT, variantID, SYMBOL, annotation,
peak_chr, peak_start, peak_end, peak_id,
n_samples, n_00, n_01, n_11,
beta_additive, r_squared_additive, p_value_additive,
r_squared_factor, p_value_factor,
fdr_additive, fdr_factor
```

### Validation run results

- 30 samples matched (14 Primary, 16 Resistant).
- 97 SNPs read → **96 SNP–peak pairs from 96 SNPs; 1 SNP had no overlapping peak**.
- RPGC: effective genome size 2,860,760,506 bp; mean fragment length 150.1 bp.
- 0 of 96 result rows are `NA`.
- Top associations by `p_value_additive`:

  | variantID  | SYMBOL    | annotation        | locus              | β       | R²    | p       |
  | ---------- | --------- | ----------------- | ------------------ | ------- | ----- | ------- |
  | rs11668945 | NCCRP1    | Promoter (2-3kb)  | chr19:39194052 T>G | 73.54   | 0.375 | 3.25e-4 |
  | rs907868   | SDC1      | Distal Intergenic | chr2:20171406 C>T  | 78.83   | 0.290 | 2.14e-3 |
  | rs79328924 | LINC01000 | Promoter (<=1kb)  | chr7:128531226 C>G | -173.29 | 0.257 | 4.21e-3 |
  | rs2257153  | SMIM8     | Promoter (<=1kb)  | chr6:87322684 C>A  | -120.62 | 0.234 | 6.71e-3 |
  | rs7775082  | CLIC5     | Promoter (<=1kb)  | chr6:46129892 T>C  | 202.12  | 0.175 | 2.14e-2 |

## Outstanding Issues / Next Steps

- **The job has not been submitted to LSF.** All validation runs were executed
  interactively on the login node with `LSB_DJOB_NUMPROC=9`. Submit with
  `bsub < bsub_run_ASoC_peak_genotype_lm.sh` from `Steffi_works/` to use the
  full 32-slot allocation.
- **No multiple-testing survivors were checked.** `fdr_additive` /
  `fdr_factor` are computed but the run log did not report how many pass a
  given threshold. Worth summarising.
- **Category is not modelled.** `category` (Primary/Resistant) is currently only
  a plot aesthetic. If a genotype effect adjusted for treatment status is
  wanted, add `lm(rpgc ~ dosage + category)` or a `dosage * category`
  interaction term.
- **Peak width is not normalised out.** Per the user's decision, RPGC is a
  scaled fragment count, not per-bp coverage, so values are not comparable
  across peaks of different widths. Only within-peak comparisons across
  genotypes are valid — which is all the current models do.
- **Fragment length is estimated from `chr1` + `chr2` only.** Adequate for a
  scale factor but could be widened if a more precise estimate is needed.
- **Sample size is small (n = 30).** With genotype classes as small as 3–5
  samples, individual `lm` fits are underpowered; `min_samples_per_fit = 6L`
  only guards against degenerate fits, not low power.

## Context for LLM Handoff

This session built a complete, working R pipeline at
`Multiome_main/Steffi_works/run_ASoC_peak_genotype_lm.R` (plus the LSF wrapper
`bsub_run_ASoC_peak_genotype_lm.sh`) that tests whether hepatocyte ATAC peak
accessibility depends on donor genotype at ASoC SNPs. The pipeline reads 97 SNPs
from `hepatocyte_ASoC_genotypes_GT_summary_min3.tsv`, intersects them with the
`PeakMatrix` peaks of the `ArchR_hepato` ArchR project (96 SNPs hit a peak),
sums each peak's fragment counts per sample across that sample's hepatocyte
cells, converts to deepTools-style RPGC
(`counts × effective_genome_size / (nFrags_sample × mean_frag_len_sample)`,
no peak-width normalisation by explicit user choice), then fits two linear
models per peak — additive genotype dosage and a 3-level genotype factor — using
the per-sample calls in `hepatocyte_ASoC_genotypes_GT_only.tsv`. Genotypes are
unphased so `1/0` is normalised to `0/1`; non-biallelic or missing calls drop
that sample from the SNP. Results (β, R², p, BH FDR) plus `variantID`, `SYMBOL`
and `annotation` joined from `sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv`
are written to `hepatocyte_ASoC_peak_genotype_lm/ASoC_peak_genotype_lm_results.tsv`
and `.qs2`, and each peak gets a box-and-whisker panel (outliers hidden, `Set2`
fill, jittered per-sample points shaped by Primary/Resistant, a `shape = 4`
`darkred` mean marker) laid out 6×6 per landscape 11 × 8.5 in page, one PDF per
page. Parallelism uses PSOCK + `doParallel` over samples for Arrow reads and
over pages for PDF writing. The one non-obvious pitfall encountered and fixed:
`ArchR:::.getMatFromArrow` reorders its result rows to match the supplied
`featureDF` and then sets `rownames(mat) <- NULL`, so row sums must be indexed
positionally — an earlier name-based lookup silently yielded all-zero counts and
an all-`NA` results table. The pipeline has been validated end-to-end
interactively (exit 0, 96 results, 3 PDF pages visually inspected) but has not
yet been submitted to LSF.
