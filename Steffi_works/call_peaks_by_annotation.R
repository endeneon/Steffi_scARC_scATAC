#! /usr/bin/env Rscript
# use the conda env jwen_scRNA_singCellaR
# init
{
  library(Seurat)
  library(gplots)
  library(ArchR)
  library(future)
  library(stringr)
  # library(pheatmap)

  library(BiocParallel)
  # library(BiocParallel.FutureParam)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(doFuture)
  library(snow)

  library(Matrix)
  library(matrixStats)

  library(qs2)
  library(fs)

  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg38)

  library(rtracklayer)
  library(ChIPseeker)
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome)
  library(motifbreakR)
  library(MotifDb)
  library(dplyr)

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


setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works"
)
# determine if R is running in RSTUDIO/VSCode/Positron
if (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode")) {
  print("Running under RStudio/VSCode/Positron IDE, use plan(multisession)")
  session_plan <- "multisession"
} else {
  print("Running under Rscript, use plan(multicore)")
  session_plan <- "multicore"
}

# preload functions ####
get_available_workers <-
  function(x) {
    future::plan(session_plan) # check here!
    return(future::nbrOfFreeWorkers())
  }

# LSF does not expose its core allocation to future::availableCores() by
# default, so it falls back to reporting 1. Read LSB_DJOB_NUMPROC directly,
# falling back to parallelly's detection when not running under LSF.
lsf_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = NA))
available_cores <-
  if (!is.na(lsf_cores) && lsf_cores >= 1) {
    lsf_cores
  } else {
    parallelly::availableCores()
  }
print(paste0(
  "Detected ",
  available_cores,
  " available cores (LSF_DJOB_NUMPROC = ",
  lsf_cores,
  ")."
))

# OpenMP-backed code (e.g. Rtsne) is single-process / shared-memory, so it can
# only use cores that live on ONE host. Keep all LSF slots on a single node
# (`#BSUB -R "span[hosts=1]"`), otherwise the -n 40 slots get spread across hosts
# and only the master host's share is actually usable.
# Rtsne() also IGNORES OMP_NUM_THREADS unless you pass num_threads = 0; its
# num_threads argument (default 1) otherwise wins -- which is why the log showed
# "OpenMP is working. 1 threads.". We therefore set an explicit thread count and
# forward it to each RunTSNE(num_threads = omp_threads) call below. Cap it (e.g.
# `min(8L, ...)`) if you prefer fewer threads than the full allocation.
omp_threads <-
  min(
    8L,
    max(
      1L,
      as.integer(available_cores / 2)
    )
  )
Sys.setenv(OMP_NUM_THREADS = omp_threads)
print(paste0(
  "Detected ",
  omp_threads,
  " OpenMP threads (OMP_NUM_THREADS = ",
  Sys.getenv("OMP_NUM_THREADS"),
  ")."
))

if (session_plan == "multisession") {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        16
      )
    )
  doFuture::registerDoFuture()
} else {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        32
      )
    )
  # multicore_workers <- MulticoreParam(
  #   workers = workers_2_use - 1,                # Number of allocated CPU cores
  #   progressbar = TRUE,         # Show visual progress bars
  #   stop.on.error = TRUE        # Halt execution if a core errors out
  # )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  options(useUCSCChromosomeNames = FALSE)
  set.seed(42)
  # Seurat's parallelized steps (e.g. IntegrateLayers / FindIntegrationAnchors)
  # call future_lapply() without future.seed = TRUE, so future warns about
  # "unreliable" RNG. We cannot pass future.seed through Seurat, and our
  # stochastic steps are already explicitly seeded (set.seed(42) above, plus
  # seed.use = 42 / random.seed = 42 on each Seurat call), so silence the check.
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = workers_2_use * 20 * 1024^3) # 20 G per thread
  future::plan(
    session_plan, # Do NOT use "multisession" here if use LSF, use "multicore" instead
    workers = workers_2_use
  )
}
print(paste0(
  "R session plan set to ",
  session_plan,
  " with ",
  workers_2_use,
  " workers."
))

# ArchR settings ####
addArchRThreads(threads = workers_2_use)
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

# getwd()
projMultiome <-
  loadArchRProject(
    path = "ArchR_multiome_obj",
    showLogo = FALSE
  )

head(getCellNames(projMultiome))

seurat_annotated <-
  qs_read("annotated_seurat.qs2")
head(colnames(seurat_annotated))
table(seurat_annotated$celltype)

colnames(seurat_annotated) <-
  sub(
    "([^_])_(?!_)",
    "\\1#",
    colnames(seurat_annotated),
    perl = TRUE
  )
rownames(seurat_annotated@meta.data)
length(getCellNames(projMultiome)) # 205540
ncol(seurat_annotated) # 165115
sum(colnames(seurat_annotated) %in% getCellNames(projMultiome))

# Shared cells, ordered by Seurat's current column order
shared_cells <-
  colnames(seurat_annotated)[
    colnames(seurat_annotated) %in% getCellNames(projMultiome)
  ]

seurat_shared <- seurat_annotated[, shared_cells]
projMultiome_shared <- projMultiome[shared_cells, ]

cells <- getCellNames(projMultiome_shared)

for (col in c("celltype", "celltype_fine", "treatment_group")) {
  projMultiome_shared <- addCellColData(
    ArchRProj = projMultiome_shared,
    data = as.character(seurat_shared@meta.data[cells, col]),
    name = col,
    cells = cells,
    force = TRUE
  )
}

getGenomeAnnotation(projMultiome_shared)
#  B/Plasma cell        Cycling    Endothelial Fibroblast/HSC     Hepatocyte     Macrophage       Monocyte         T cell      T/NK cell
#           3446           9484          24749           4881          54229          23701           3815          15712          11838
# Prepare for peak calling by adding a pseudo-bulk replicate for each cluster, using the cluster_name_4_peak_calling determined above.
projMultiome_shared <-
  ArchR::addGroupCoverages(
    ArchRProj = projMultiome_shared,
    groupBy = "celltype",
    maxCells = 5000,
    excludeChr = c("chrM", "chrX", "chrY"),
    force = TRUE
  )

pathToMacs2 <- findMacs2()
projMultiome_annotated_celltype <-
  ArchR::addReproduciblePeakSet(
    ArchRProj = projMultiome_shared,
    groupBy = "celltype",
    pathToMacs2 = pathToMacs2,
    excludeChr = c("chrM", "chrX", "chrY"),
    force = TRUE
  )
projMultiome_annotated_celltype <-
  saveArchRProject(
    ArchRProj = projMultiome_annotated_celltype,
    outputDirectory = "ArchR_multiome_annotated_celltype_obj",
    load = TRUE,
    overwrite = TRUE
  )
print("Reproducible peak set added and ArchRProject saved successfully.")

#####
projMultiome_annotated_celltype <-
  loadArchRProject(
    path = "ArchR_multiome_annotated_celltype_obj",
    showLogo = FALSE
  )

# extract the new peak sets from the ArchR object as GRanges
peak_sets <- getPeakSet(projMultiome_annotated_celltype)
peak_sets$Group <-
  sub(
    "\\._\\..*$",
    "",
    peak_sets$GroupReplicate
  )
peaks_by_group <- split(peak_sets, peak_sets$Group)
qs_save(
  peaks_by_group,
  "peaks_by_group.qs2",
  nthreads = 8
)
# Peaks per group
lengths(peaks_by_group)

# unique(peak_sets$GroupReplicate)
# peak_sets@metadata$PeakCallSummary$Group

out_dir <- "peak_bed_by_group"
dir.create(out_dir, showWarnings = FALSE)

for (grp in names(peaks_by_group)) {
  # Sanitize group name for a safe filename (e.g. "T/NK cell" -> "T_NK_cell")
  fname <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9]+", "_", grp), ".bed"))
  export(peaks_by_group[[grp]], fname, format = "BED")
}

ASoC_df <-
  qs_read("data_dir/ASoC_df.qs2")

ASoC_gr <- GRanges(
  seqnames = ASoC_df$chromosome,
  ranges = IRanges(start = as.integer(ASoC_df$position), width = 1),
  variantID = ASoC_df$variantID,
  refAllele = ASoC_df$refAllele,
  altAllele = ASoC_df$altAllele,
  norm_refCount = as.numeric(ASoC_df$norm_refCount),
  norm_altCount = as.numeric(ASoC_df$norm_altCount),
  norm_sumCount = as.numeric(ASoC_df$norm_sumCount),
  percRef = as.numeric(ASoC_df$percRef),
  pVal = as.numeric(ASoC_df$pVal),
  adjPVal = as.numeric(ASoC_df$adjPVal)
)

names(peaks_by_group) <- gsub("[^A-Za-z0-9]+", "_", names(peaks_by_group))

library(BSgenome.Hsapiens.UCSC.hg38)

hg38_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
hg38_seqinfo <- keepSeqlevels(hg38_seqinfo, seqlevels(peak_sets))

seqlevels(ASoC_gr) <- seqlevels(hg38_seqinfo)
seqinfo(ASoC_gr) <- hg38_seqinfo

seqinfo(peak_sets) <- hg38_seqinfo
seqinfo(peaks_by_group) <- hg38_seqinfo

# Shared annotation resources (computed once, reused for every cell type) ####
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# JASPAR2020 human CORE PWMs for motifmatchr
pwm_set <- getMatrixSet(
  JASPAR2020,
  opts = list(species = 9606, collection = "CORE")
)
length(pwm_set)

# Human HOCOMOCO motifs from MotifDb for motifbreakR
motifs_human <- MotifDb::query(MotifDb, andStrings = c("hsapiens", "HOCOMOCO"))
length(motifs_human)

# motifbreakR parallel backend: SerialParam() when run interactively in an IDE,
# MulticoreParam() when run as a batch Rscript from the command line.
mb_bpparam <-
  if (
    Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode")
  ) {
    BiocParallel::SerialParam()
  } else {
    BiocParallel::MulticoreParam()
  }

# Output directory for per-cell-type results
out_dir_asoc <- "sig_ASoC_by_celltype"
dir.create(out_dir_asoc, showWarnings = FALSE, recursive = TRUE)

# Accumulate one summary row per cell type for a combined cross-cell-type table
summary_rows <- list()

# Iterate over every cell type, annotate its significant ASoC variants, and
# save one qs2 + one TSV per cell type ####
for (celltype in peaks_by_group@partitioning@NAMES) {
  message("Processing cell type: ", celltype)

  # Overlap ASoC variants with this cell type's peak set
  hits <- findOverlaps(ASoC_gr, peaks_by_group[[celltype]])

  ASoC_in_peaks <- ASoC_gr[queryHits(hits)]
  ASoC_in_peaks$peak_start <- start(peaks_by_group[[celltype]])[subjectHits(
    hits
  )]
  ASoC_in_peaks$peak_end <- end(peaks_by_group[[celltype]])[subjectHits(hits)]

  # Keep only significant ASoC variants
  sig_ASoC_in_peaks <- ASoC_in_peaks[ASoC_in_peaks$adjPVal < 0.05]

  if (length(sig_ASoC_in_peaks) == 0) {
    message("  No significant ASoC variants in ", celltype, " peaks; skipping.")
    summary_rows[[celltype]] <- data.frame(
      celltype = celltype,
      n_peaks = length(peaks_by_group[[celltype]]),
      n_ASoC_in_peaks = length(ASoC_in_peaks),
      n_sig_ASoC = 0L,
      n_with_motif_match = 0L,
      n_with_any_disruption = 0L,
      n_with_strong_disruption = 0L,
      pct_strong_disruption = NA_real_
    )
    next
  }

  # ChIPseeker genomic annotation
  sig_ASoC_anno <- annotatePeak(
    sig_ASoC_in_peaks,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
  sig_ASoC_anno_df <- as.data.frame(sig_ASoC_anno)
  mcols(sig_ASoC_in_peaks) <- cbind(
    as.data.frame(mcols(sig_ASoC_in_peaks)),
    sig_ASoC_anno_df[, c(
      "annotation",
      "geneId",
      "transcriptId",
      "distanceToTSS",
      "ENSEMBL",
      "SYMBOL",
      "GENENAME"
    )]
  )

  # motifmatchr: which JASPAR motifs are present in the surrounding peak
  peak_ranges <- GRanges(
    seqnames = seqnames(sig_ASoC_in_peaks),
    ranges = IRanges(
      start = sig_ASoC_in_peaks$peak_start,
      end = sig_ASoC_in_peaks$peak_end
    ),
    seqinfo = seqinfo(sig_ASoC_in_peaks)
  )
  motif_matches <- matchMotifs(
    pwm_set,
    peak_ranges,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    out = "matches"
  )
  match_mat <- motifMatches(motif_matches)
  motif_names <- colData(motif_matches)$name
  motif_list <- apply(match_mat, 1, function(row) motif_names[row])
  sig_ASoC_in_peaks$n_motifs_matched <- lengths(motif_list)
  sig_ASoC_in_peaks$motif_names <-
    vapply(motif_list, function(x) paste(x, collapse = ";"), character(1))

  # motifbreakR: allele-specific motif disruption.
  # motifbreakR expects REF/ALT DNAStringSet columns, hg38 seqlengths, and an
  # attributes(snpList)$genome.package pointing to the BSgenome package;
  # without the last one motifbreakR() errors with "subscript out of bounds".
  snp_gr_min <- GRanges(
    seqnames = seqnames(sig_ASoC_in_peaks),
    ranges = ranges(sig_ASoC_in_peaks),
    strand = "+"
  )
  names(snp_gr_min) <- sig_ASoC_in_peaks$variantID
  snp_gr_min$REF <- Biostrings::DNAStringSet(sig_ASoC_in_peaks$refAllele)
  snp_gr_min$ALT <- Biostrings::DNAStringSet(sig_ASoC_in_peaks$altAllele)
  snp_gr_min$SNP_id <- sig_ASoC_in_peaks$variantID
  genome(snp_gr_min) <- "hg38"
  seqinfo(snp_gr_min) <- hg38_seqinfo
  attributes(snp_gr_min)$genome.package <- "BSgenome.Hsapiens.UCSC.hg38"

  set.seed(4827)
  mb_results <- motifbreakR(
    snpList = snp_gr_min,
    pwmList = motifs_human,
    filterp = TRUE,
    threshold = 1e-4,
    method = "ic",
    bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    BPPARAM = mb_bpparam
  )

  # Per-SNP summary of motif disruption
  mb_df <- as.data.frame(mb_results, row.names = NULL) |>
    as_tibble() |>
    mutate(SNP_id = names(mb_results))
  mb_summary <- mb_df |>
    group_by(SNP_id) |>
    summarise(
      n_motifs_disrupted = n(),
      n_strong = sum(effect == "strong"),
      disrupted_TFs = paste(
        sort(unique(geneSymbol[effect == "strong"])),
        collapse = ";"
      ),
      max_abs_alleleDiff = max(abs(alleleDiff)),
      .groups = "drop"
    )

  # Integrate motifbreakR summary back into the GRanges by variantID
  m <- match(sig_ASoC_in_peaks$variantID, mb_summary$SNP_id)
  sig_ASoC_in_peaks$mb_n_motifs_disrupted <- mb_summary$n_motifs_disrupted[m]
  sig_ASoC_in_peaks$mb_n_strong <- mb_summary$n_strong[m]
  sig_ASoC_in_peaks$mb_disrupted_TFs <- mb_summary$disrupted_TFs[m]
  sig_ASoC_in_peaks$mb_max_abs_alleleDiff <- mb_summary$max_abs_alleleDiff[m]

  # SNPs overlapping no motif get 0 counts (score-diff columns remain NA)
  sig_ASoC_in_peaks$mb_n_motifs_disrupted[is.na(
    sig_ASoC_in_peaks$mb_n_motifs_disrupted
  )] <- 0L
  sig_ASoC_in_peaks$mb_n_strong[is.na(sig_ASoC_in_peaks$mb_n_strong)] <- 0L

  # Save one qs2 + one TSV per cell type (sanitized names already applied above)
  out_base <- file.path(
    out_dir_asoc,
    paste0("sig_ASoC_in_", celltype, "_annotated")
  )
  qs_save(
    sig_ASoC_in_peaks,
    paste0(out_base, ".qs2"),
    nthreads = 8
  )
  write.table(
    as.data.frame(sig_ASoC_in_peaks, row.names = NULL),
    file = paste0(out_base, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Record a summary row for the combined cross-cell-type table
  summary_rows[[celltype]] <- data.frame(
    celltype = celltype,
    n_peaks = length(peaks_by_group[[celltype]]),
    n_ASoC_in_peaks = length(ASoC_in_peaks),
    n_sig_ASoC = length(sig_ASoC_in_peaks),
    n_with_motif_match = sum(sig_ASoC_in_peaks$n_motifs_matched > 0),
    n_with_any_disruption = sum(sig_ASoC_in_peaks$mb_n_motifs_disrupted > 0),
    n_with_strong_disruption = sum(sig_ASoC_in_peaks$mb_n_strong > 0),
    pct_strong_disruption = round(
      100 * mean(sig_ASoC_in_peaks$mb_n_strong > 0),
      1
    )
  )

  message(
    "  ",
    celltype,
    ": ",
    length(sig_ASoC_in_peaks),
    " significant variants; ",
    sum(sig_ASoC_in_peaks$mb_n_strong > 0),
    " with a strong motif disruption. Saved to ",
    out_base,
    ".{qs2,tsv}"
  )
}

# Combined cross-cell-type summary table ####
celltype_summary <- do.call(rbind, summary_rows)
rownames(celltype_summary) <- NULL
celltype_summary <- celltype_summary[
  order(
    -celltype_summary$n_sig_ASoC
  ),
]

write.table(
  celltype_summary,
  file = file.path(out_dir_asoc, "sig_ASoC_celltype_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
celltype_summary
