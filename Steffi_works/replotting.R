#! /usr/bin/env Rscript

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
  library(Rsamtools)

  library(ggplot2)
  library(Gviz)
  # library(Gviz)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(org.Hs.eg.db)

  if (
    interactive() &&
      (Sys.getenv("TERM_PROGRAM") == "vscode") &&
      (Sys.getenv("POSITRON") != "1")
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
  future::plan(future::multisession, workers = workers_2_use)
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
  future::plan(future::multicore, workers = workers_2_use)
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
}
print(paste0(
  "R session plan set to ",
  session_plan,
  " with ",
  workers_2_use,
  " workers."
))

# ArchR settings ####
if (
  interactive() &&
    (Sys.getenv("TERM_PROGRAM") == "vscode")
) {
  print("Running under IDE, use 1 ArchR Thread")
  addArchRThreads(threads = 1)
} else {
  print("Running under Rscript, use all usable ArchR Threads")
  addArchRThreads(threads = workers_2_use)
}
addArchRGenome("hg38")
print(paste0(
  "ArchR threads set to ",
  getArchRThreads(),
  " and genome set to ",
  getArchRGenome()
))
print("All settings initialized successfully.")

# projMacrophages <-
#   ArchR::loadArchRProject(path = "ArchR_macrophages")

# seurat_multiome_obj <-
#   qs2::qs_read(
#     "annotation_script_package/marker_based/annotated_seurat_marker_based.qs2",
#     nthreads = 8
#   )
# SeuratObject::Idents(seurat_multiome_obj) <- "celltype_broad"
# seurat_macrophage_obj <-
#   subset(
#     seurat_multiome_obj,
#     cells = SeuratObject::WhichCells(
#       seurat_multiome_obj,
#       idents = "Macrophage"
#     )
#   )
# SeuratObject::Idents(seurat_macrophage_obj) <- "preparation"
# seurat_macrophage_obj <-
#   subset(
#     seurat_macrophage_obj,
#     cells = SeuratObject::WhichCells(
#       seurat_macrophage_obj,
#       idents = "multiome"
#     )
#   )

# projMultiome <-
#   ArchR::loadArchRProject(
#     path = "ArchR_merged_ATAC_multiome_obj_multiple_cell_types",
#     showLogo = FALSE
#   )
# ArchR_macrophage_obj <-
#   ArchR::subsetArchRProject(
#     ArchRProj = projMultiome,
#     cells = projMultiome$cellNames[
#       projMultiome$projected_barcodes_multiple_cell_types == "Macrophage"
#     ],
#     outputDirectory = "ArchR_macrophages_projected",
#     dropCells = TRUE,
#     logFile = createLogFile("ArchR_macrophages_projected", logDir = "log")
#   )
ArchR_macrophage_obj <-
  ArchR::loadArchRProject(
    path = "ArchR_merged_ATAC_multiome_obj_multiple_cell_types",
    showLogo = FALSE
  )

df_genotype_raw <-
  read.table(
    "scQTLtools/INPP4A_ASoC_genotypes_extended_GT_only.tsv",
    header = TRUE,
    sep = "\t"
  )

df_metadata <-
  as.data.frame(ArchR_macrophage_obj@cellColData)

df_metadata <-
  df_metadata[!duplicated(df_metadata$Sample), , drop = FALSE]
df_metadata <-
  df_metadata[order(df_metadata$Sample), , drop = FALSE]

colnames(df_genotype_raw) <-
  stringr::str_replace(
    colnames(df_genotype_raw),
    "multiome_",
    ""
  )
colnames(df_genotype_raw) <-
  stringr::str_replace(
    colnames(df_genotype_raw),
    "atac_",
    ""
  )
colnames(df_genotype_raw) <-
  stringr::str_replace(
    colnames(df_genotype_raw),
    "_GT",
    ""
  )

df_genotype_multiome <- df_genotype_raw
df_genotype_multiome <-
  as.data.frame(t(df_genotype_multiome[, -c(1:4)])) # transpose to have samples as columns
rownames(df_genotype_multiome) <- colnames(df_genotype_raw)[-c(1:4)]
df_genotype_multiome <-
  df_genotype_multiome[, order(colnames(df_genotype_multiome)), drop = FALSE]
df_genotype_multiome$Sample <-
  rownames(df_genotype_multiome) # add Sample column for merging with metadata
colnames(df_genotype_multiome)[1:2] <-
  as.vector(stringr::str_c(
    df_genotype_raw[, 1],
    df_genotype_raw[, 2],
    df_genotype_raw[, 3],
    df_genotype_raw[, 4],
    sep = ":"
  )) # rename first two columns to match df_genotype_raw

df_metadata$Sample <-
  stringr::str_replace(
    df_metadata$Sample,
    "X__",
    ""
  )
all(df_metadata$Sample %in% df_genotype_multiome$Sample) # TRUE

df_2_plot <-
  merge(
    df_metadata,
    df_genotype_multiome,
    by = "Sample",
    all.x = TRUE,
    sort = FALSE
  )

df_2_plot$category <- "Primary"
df_2_plot$category[str_detect(
  as.character(df_2_plot$Sample),
  pattern = "pre$",
  negate = FALSE
)] <- "Resistant"

colnames(df_2_plot)[(ncol(df_2_plot) - 2):(ncol(df_2_plot) - 1)] <-
  c("PICALM", "INPP4A")

df_2_plot |>
  dplyr::count(INPP4A, category) |>
  # keep empty combinations (e.g. 1/1 Resistant) as 0 so every x has two bars
  tidyr::complete(INPP4A, category, fill = list(n = 0)) |>
  ggplot(aes(x = INPP4A, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.3
  ) +
  scale_fill_manual(values = brewer.pal(n = 8, name = "Dark2")[c(1, 2)]) +
  theme_bw()


rownames(df_genotype_multiome) <-
  colnames(df_genotype_multiome) <-
    stringr::str_replace(
      colnames(df_genotype_multiome),
      "multiome_",
      ""
    )
colnames(df_genotype_multiome) <-
  stringr::str_replace(
    colnames(df_genotype_multiome),
    "atac_",
    ""
  )
colnames(df_genotype_multiome) <-
  stringr::str_replace(
    colnames(df_genotype_multiome),
    "_GT",
    ""
  )

# Recode genotype strings to integer codes; anything unmatched (incl. NA) -> 0.
# 1 = homozygous reference (0/0), 2 = homozygous alternate (1/1),
# 3 = heterozygous (0/1 or 1/0), 0 = missing / unknown
df_genotype_multiome[] <- lapply(
  df_genotype_multiome,
  function(col) {
    dplyr::case_when(
      col == "0/0" ~ 1L,
      col %in% c("0/1", "1/0") ~ 3L,
      col == "1/1" ~ 2L,
      .default = 0L
    )
  }
)

all_barcode_list <-
  rownames(ArchR_macrophage_obj@cellColData)
df_full_genotype <-
  as.data.frame(matrix(
    0,
    nrow = nrow(df_genotype_multiome),
    ncol = length(all_barcode_list)
  ))
rownames(df_full_genotype) <- rownames(df_genotype_multiome)
colnames(df_full_genotype) <- all_barcode_list

# Fill each cell column from its sample's genotype column. The sample prefix is
# everything before the final "_<barcode>" segment of the cell name. Cells whose
# sample has no genotype column (e.g. X__20t) are dropped from df_full_genotype.
cell_sample <- sub("#[^_]+$", "", colnames(df_full_genotype))
cell_sample <-
  stringr::str_replace_all(
    cell_sample,
    pattern = "X__",
    replacement = ""
  ) # remove multiome_ prefix
# cell_sample <- colnames(df_full_genotype)
keep <- cell_sample %in% colnames(df_genotype_multiome)
df_full_genotype <- df_full_genotype[, keep, drop = FALSE]
df_full_genotype[] <- df_genotype_multiome[, cell_sample[keep], drop = FALSE]
all(colnames(ArchR_macrophage_obj) == colnames(df_full_genotype)) # TRUE

colnames(df_full_genotype) <- stringr::str_replace(
  colnames(df_full_genotype),
  "X__",
  ""
)
# Sample ID = prefix before the first "_" in each column name
sample_ids <- str_split_i(colnames(df_full_genotype), "\\#", 1)
# Named list of data frames, one per sample
list_genotype_by_sample <- lapply(
  split(seq_along(sample_ids), sample_ids),
  function(idx) df_full_genotype[, idx, drop = FALSE]
)
list_genotype_by_sample <- lapply(
  list_genotype_by_sample,
  function(df) {
    colnames(df) <- str_split_i(colnames(df), "#", 2)
    colnames(df) <- str_split_i(colnames(df), "-", 1)
    return(df)
  }
)

list_genotype_by_sample <-
  list_genotype_by_sample[order(names(list_genotype_by_sample))] # sort by sample name

# qs2::qs_save(
#   list_genotype_by_sample,
#   file = "scQTLtools/list_genotype_by_sample_macrophages_all.qs2",
#   nthreads = 8
# )

WASPed_bam_list <-
  list.files(
    path = "test_ASoC_w_WASP/WASPed_bams",
    pattern = "WASPed.bam$",
    full.names = TRUE,
    recursive = TRUE
  )
WASPed_bam_list <-
  WASPed_bam_list[order(basename(WASPed_bam_list))] # sort by sample name
names(WASPed_bam_list) <- stringr::str_split_i(
  basename(WASPed_bam_list),
  "_bwa",
  1
)
# intersect(names(WASPed_bam_list), names(list_genotype_by_sample))
WASPed_bam_list <-
  WASPed_bam_list[names(WASPed_bam_list) %in% names(list_genotype_by_sample)] # keep only samples with genotype data
all(names(list_genotype_by_sample) == names(WASPed_bam_list)) # TRUE

# Allele-specific RPGC pileups at heterozygous SNPs ####
library(Rsamtools)
library(GenomicAlignments)
library(data.table)
library(future.apply)

bam_split_dir <- "/lustre_scratch/user_scratch/szhang37/bam_split"
pileup_out_dir <- "allelic_pileups"
dir.create(bam_split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pileup_out_dir, recursive = TRUE, showWarnings = FALSE)

effective_genome_size <- 2913022398 # deepTools hg38 EGS
# Window fetched around each SNP so the non-covering mate of a pair is also read
fetch_flank <- 1000L
plot_flank <- 50L

read_flag <- scanBamFlag(
  isPaired = TRUE,
  isProperPair = TRUE,
  isUnmappedQuery = FALSE,
  hasUnmappedMate = FALSE,
  isSecondaryAlignment = FALSE,
  isNotPassingQualityControls = FALSE,
  isDuplicate = FALSE,
  isSupplementaryAlignment = FALSE
)

stopifnot(all(names(list_genotype_by_sample) %in% names(WASPed_bam_list)))

# Step 2: keep only reads from the sample's selected cells, restricted to the
# het-SNP windows (+/- fetch_flank) since nothing outside them is read later.
# BAM CB tags carry a "-<n>" GEM-well suffix that the genotype colnames do not.
subset_bam_by_cells <- function(sample_name) {
  geno <- list_genotype_by_sample[[sample_name]]
  het_snps <- rownames(geno)[rowMeans(geno) == 3]
  if (!length(het_snps)) {
    return(NULL)
  }
  parts <- stringr::str_split_fixed(het_snps, ":", 4)
  pos <- as.integer(parts[, 2])
  # reduce() merges overlapping windows so no read is written twice
  windows <- reduce(GRanges(
    parts[, 1],
    IRanges(pmax(1L, pos - fetch_flank), pos + fetch_flank)
  ))

  barcodes <- colnames(geno)
  dest <- file.path(bam_split_dir, paste0(sample_name, "_cellsubset.bam"))
  unlink(c(dest, paste0(dest, ".bai")))
  filterBam(
    WASPed_bam_list[[sample_name]],
    destination = dest,
    filter = FilterRules(list(
      cell = function(x) sub("-\\d+$", "", x$CB) %in% barcodes
    )),
    param = ScanBamParam(tag = "CB", which = windows),
    indexDestination = TRUE
  )
}

# RPGC scale factor: EGS / (filtered mapped reads x mean aligned read length),
# using the original (full) WASPed BAM as the parent.
rpgc_scale_factor <- function(full_bam) {
  n_mapped <- countBam(full_bam, param = ScanBamParam(flag = read_flag))$records
  bf <- BamFile(full_bam, yieldSize = 1e5)
  open(bf)
  on.exit(close(bf))
  read_len <- mean(width(readGAlignments(
    bf,
    param = ScanBamParam(flag = read_flag)
  )))
  effective_genome_size / (n_mapped * read_len)
}

# Base on the reference coordinate `pos` for each alignment (NA if not covered
# or deleted). sequenceLayer() projects the read onto reference space.
base_at_pos <- function(ga, pos) {
  out <- rep(NA_character_, length(ga))
  hit <- start(ga) <= pos & end(ga) >= pos
  if (any(hit)) {
    layer <- sequenceLayer(mcols(ga)$seq[hit], cigar(ga)[hit])
    out[hit] <- as.character(Biostrings::subseq(
      layer,
      start = pos - start(ga)[hit] + 1L,
      width = 1L
    ))
  }
  out[out %in% c("-", ".")] <- NA_character_
  out
}

# Per-bp read coverage of a set of pairs; overlapping mates are merged so the
# shared stretch is counted once per fragment.
pair_coverage_bins <- function(gal, chr) {
  mate_ranges <- unlist(reduce(grglist(gal), ignore.strand = TRUE))
  cov_gr <- as(coverage(mate_ranges)[chr], "GRanges")
  cov_gr <- cov_gr[cov_gr$score > 0]
  w <- width(cov_gr)
  data.table(
    chrom = chr,
    pos = unlist(mapply(seq, start(cov_gr), end(cov_gr), SIMPLIFY = FALSE)),
    count = rep(cov_gr$score, w)
  )
}

# Step 3: split each heterozygous SNP's reads by allele, write allele BAMs and
# RPGC-normalised 1-bp BED6 files. Returns a manifest of files written.
split_sample_by_allele <- function(sample_name, subset_bam) {
  geno <- list_genotype_by_sample[[sample_name]]
  het_snps <- rownames(geno)[rowMeans(geno) == 3]
  parts <- stringr::str_split_fixed(het_snps, ":", 4)
  is_snv <- nchar(parts[, 3]) == 1 & nchar(parts[, 4]) == 1
  if (any(!is_snv)) {
    message(sample_name, ": skipping ", sum(!is_snv), " non-SNV het site(s)")
  }
  if (!any(is_snv)) {
    return(NULL)
  }

  scale_factor <- rpgc_scale_factor(WASPed_bam_list[[sample_name]])

  rbindlist(lapply(which(is_snv), function(i) {
    snp <- het_snps[i]
    chr <- parts[i, 1]
    pos <- as.integer(parts[i, 2])
    ref <- parts[i, 3]
    alt <- parts[i, 4]
    fetch_gr <- GRanges(
      chr,
      IRanges(max(1L, pos - fetch_flank), pos + fetch_flank)
    )

    gal <- readGAlignmentPairs(
      subset_bam,
      use.names = TRUE,
      param = ScanBamParam(which = fetch_gr, flag = read_flag, what = "seq")
    )
    if (!length(gal)) {
      return(NULL)
    }

    call_base <- function(b) {
      ifelse(b == ref, "ref", ifelse(b == alt, "alt", "other"))
    }
    c1 <- call_base(base_at_pos(first(gal), pos))
    c2 <- call_base(base_at_pos(second(gal), pos))
    # Pair allele = covering mate(s); discard discordant mates and non-ref/alt
    pair_call <- dplyr::coalesce(c1, c2)
    pair_call[!is.na(c1) & !is.na(c2) & c1 != c2] <- NA
    pair_call[pair_call %in% "other"] <- NA

    file_stub <- paste(sample_name, gsub(":", "_", snp), sep = "__")

    rbindlist(lapply(c("ref", "alt"), function(allele) {
      keep_pairs <- which(pair_call %in% allele)
      if (!length(keep_pairs)) {
        return(NULL)
      }
      qnames <- unique(names(gal)[keep_pairs])

      allele_bam <- file.path(
        bam_split_dir,
        paste0(file_stub, "__", allele, ".bam")
      )
      unlink(c(allele_bam, paste0(allele_bam, ".bai")))
      filterBam(
        subset_bam,
        destination = allele_bam,
        filter = FilterRules(list(qn = function(x) x$qname %in% qnames)),
        param = ScanBamParam(
          which = fetch_gr,
          flag = read_flag,
          what = "qname"
        ),
        indexDestination = TRUE
      )

      bins <- pair_coverage_bins(gal[keep_pairs], chr)
      bed_file <- file.path(
        bam_split_dir,
        paste0(file_stub, "__", allele, "_RPGC.bed")
      )
      fwrite(
        data.table(
          chrom = bins$chrom,
          start = bins$pos - 1L,
          end = bins$pos,
          name = paste(sample_name, snp, allele, sep = "|"),
          score = signif(bins$count * scale_factor, 6),
          strand = "."
        ),
        bed_file,
        sep = "\t",
        col.names = FALSE
      )

      data.table(
        sample = sample_name,
        snp = snp,
        allele = allele,
        n_pairs = length(keep_pairs),
        bam = allele_bam,
        bed = bed_file
      )
    }))
  }))
}

allelic_manifest <- rbindlist(future_lapply(
  names(list_genotype_by_sample),
  function(sample_name) {
    subset_bam <- subset_bam_by_cells(sample_name)
    if (is.null(subset_bam)) {
      return(NULL)
    }
    split_sample_by_allele(sample_name, subset_bam)
  },
  future.seed = TRUE,
  future.packages = c("Rsamtools", "GenomicAlignments", "data.table")
))
fwrite(
  allelic_manifest,
  file.path(pileup_out_dir, "allelic_split_manifest_all.tsv"),
  sep = "\t"
)

# Step 4: sum RPGC bins across samples per SNP x allele. Samples that are not
# heterozygous (or have no reads for an allele) have no BED and are skipped.
master_beds <- allelic_manifest[,
  {
    bed_cols <- c("chrom", "start", "end", "name", "score", "strand")
    summed <- rbindlist(lapply(bed, fread, col.names = bed_cols))[,
      .(score = sum(score)),
      by = .(chrom, start, end)
    ]
    setorder(summed, start)
    out_file <- file.path(
      pileup_out_dir,
      paste0(gsub(":", "_", snp), "__", allele, "_merged_RPGC_all.bed")
    )
    fwrite(
      summed[, .(
        chrom,
        start,
        end,
        name = paste(snp, allele, sep = "|"),
        score,
        strand = "."
      )],
      out_file,
      sep = "\t",
      col.names = FALSE
    )
    .(master_bed = out_file, n_samples = .N)
  },
  by = .(snp, allele)
]

# Annotate each SNP with dbSNP155 rsID and nearest gene symbol (hg38 UCSC genes)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

snp_parts <- stringr::str_split_fixed(unique(master_beds$snp), ":", 4)
snp_annot <- data.table(
  snp = unique(master_beds$snp),
  chrom = snp_parts[, 1],
  pos = as.integer(snp_parts[, 2]),
  ref = snp_parts[, 3],
  alt = snp_parts[, 4]
)
snp_gr <- GRanges(snp_annot$chrom, IRanges(snp_annot$pos, width = 1))

# rsID: dbSNP uses NCBI-style seqnames ("1" rather than "chr1"). Multiple rsIDs
# at one position are collapsed with ",".
snp_gr_ncbi <- snp_gr
seqlevelsStyle(snp_gr_ncbi) <- "NCBI"
dbsnp_hits <- snpsByOverlaps(
  SNPlocs.Hsapiens.dbSNP155.GRCh38,
  snp_gr_ncbi
)
hit_map <- findOverlaps(snp_gr_ncbi, dbsnp_hits)
rsid_by_query <- tapply(
  dbsnp_hits$RefSNP_id[subjectHits(hit_map)],
  queryHits(hit_map),
  function(x) paste(unique(x), collapse = ",")
)
snp_annot[, rsid := "NA"]
snp_annot[as.integer(names(rsid_by_query)), rsid := unname(rsid_by_query)]

# Nearest gene body (distance 0 if the SNP lies within a gene)
txdb_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
nearest_idx <- nearest(snp_gr, txdb_genes, ignore.strand = TRUE)
nearest_entrez <- ifelse(
  is.na(nearest_idx),
  NA_character_,
  names(txdb_genes)[nearest_idx]
)
nearest_symbol <- suppressMessages(mapIds(
  org.Hs.eg.db,
  keys = unique(na.omit(nearest_entrez)),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
))
snp_annot[, gene_symbol := unname(nearest_symbol[nearest_entrez])]
snp_annot[is.na(gene_symbol), gene_symbol := "NA"]
# Distance in bp to the gene body: 0 if inside, otherwise gap to nearest end
nearest_gene_gr <- txdb_genes[nearest_idx]
snp_annot[,
  gene_distance_bp := as.integer(pmax(
    0L,
    start(nearest_gene_gr) - pos,
    pos - end(nearest_gene_gr)
  ))
]
snp_annot[,
  panel_title := paste(
    chrom,
    pos,
    rsid,
    gene_symbol,
    ref,
    alt,
    ifelse(is.na(gene_distance_bp), "NA", paste0(gene_distance_bp, "bp")),
    sep = ":"
  )
]
fwrite(
  snp_annot,
  file.path(pileup_out_dir, "snp_annotation_all.tsv"),
  sep = "\t"
)

# Step 5: overlay ref / alt coverage per SNP
cytoband_file <- "cytoBandIdeo_hg38.txt.gz"
cytobands <- if (file.exists(cytoband_file)) {
  read.table(
    cytoband_file,
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
    stringsAsFactors = FALSE
  )
}
plot_allelic_pileup <- function(snp_id) {
  parts <- strsplit(snp_id, ":", fixed = TRUE)[[1]]
  chr <- parts[1]
  pos <- as.integer(parts[2])
  files <- master_beds[snp == snp_id]

  read_master <- function(allele_id) {
    f <- files[allele == allele_id, master_bed]
    if (!length(f)) {
      return(GRanges())
    }
    d <- fread(
      f,
      col.names = c("chrom", "start", "end", "name", "score", "strand")
    )
    gr <- GRanges(d$chrom, IRanges(d$start + 1L, d$end), score = d$score)
    # Collapse runs of equal-value adjacent 1-bp bins into one bar. At ~1 px
    # per bp, Gviz's 1-bp rectangles overlap and semi-transparent fills stack.
    cov <- coverage(gr, weight = "score")[[unique(d$chrom)]]
    runs <- GRanges(unique(d$chrom), ranges(cov), score = runValue(cov))
    runs[runs$score > 0]
  }
  gr_ref <- read_master("ref")
  gr_alt <- read_master("alt")
  all_bins <- c(gr_ref, gr_alt)
  if (!length(all_bins)) {
    return(invisible(NULL))
  }

  from <- min(start(all_bins)) - plot_flank
  to <- max(end(all_bins)) + plot_flank
  y_max <- max(all_bins$score)

  make_track <- function(gr, colour) {
    DataTrack(
      gr,
      genome = "hg38",
      chromosome = chr,
      name = "RPGC",
      type = "histogram",
      # Transparency on the fill only; track-level alpha would also fade the
      # baseline
      fill.histogram = grDevices::adjustcolor(colour, alpha.f = 0.5),
      col.histogram = "transparent",
      alpha = 1,
      baseline = 0,
      col.baseline = "black",
      lwd.baseline = 1,
      ylim = c(0, y_max)
    )
  }
  data_tracks <- list()
  if (length(gr_ref)) {
    data_tracks <- c(data_tracks, make_track(gr_ref, "darkblue"))
  }
  if (length(gr_alt)) {
    data_tracks <- c(data_tracks, make_track(gr_alt, "darkred"))
  }

  snp_line <- HighlightTrack(
    trackList = list(OverlayTrack(data_tracks, name = "RPGC")),
    start = pos,
    end = pos,
    chromosome = chr,
    col = "black",
    fill = "transparent",
    lty = "dashed",
    inBackground = FALSE
  )

  # Draws into the current grid viewport (one grid cell): title strip on top,
  # tracks below.
  grid::pushViewport(grid::viewport(
    y = 1,
    height = grid::unit(1, "lines"),
    just = "top"
  ))
  grid::grid.text(
    snp_annot[snp == snp_id, panel_title],
    gp = grid::gpar(cex = 0.5, fontface = "bold")
  )
  grid::popViewport()
  grid::pushViewport(grid::viewport(
    y = 0,
    height = grid::unit(1, "npc") - grid::unit(1, "lines"),
    just = "bottom"
  ))
  # Ideogram from the local cytoband table (no UCSC query); skipped if missing
  ideo_track <- NULL
  if (file.exists(cytoband_file)) {
    bands_chr <- cytobands[cytobands$chrom == chr, , drop = FALSE]
    ideo_track <- tryCatch(
      IdeogramTrack(genome = "hg38", chromosome = chr, bands = bands_chr),
      error = function(e) NULL
    )
  }
  # labelPos = "below" stops Gviz alternating tick labels above/below the axis
  # (the "doubled" axis). Windows are only ~100-1000 bp, so pin ~3 ticks and
  # print plain bp positions; Gviz's auto ticks overlap at this width.
  axis_ticks <- pretty(c(from, to), n = 3)
  axis_ticks <- axis_ticks[axis_ticks >= from & axis_ticks <= to]
  axis_track <- GenomeAxisTrack(
    labelPos = "below",
    exponent = 0,
    ticksAt = axis_ticks,
    add53 = FALSE,
    add35 = FALSE,
    littleTicks = FALSE
  )
  track_list <- c(
    if (!is.null(ideo_track)) list(ideo_track),
    list(axis_track, snp_line)
  )
  # Relative heights: ideogram at half its previous (default 1) share
  track_sizes <- c(if (!is.null(ideo_track)) 0.5, 1.5, 4)

  plotTracks(
    track_list,
    sizes = track_sizes,
    from = from,
    to = to,
    chromosome = chr,
    add = TRUE,
    showTitle = FALSE,
    cex.axis = 0.5,
    col.axis = "black",
    fontcolor = "black"
  )
  grid::popViewport()
  invisible(snp_id)
}

# One multi-page Letter-landscape PDF, 4 rows x 3 columns of plots per page.
grid_nrow <- 4L
grid_ncol <- 3L
per_page <- grid_nrow * grid_ncol
snps_to_plot <- unique(master_beds$snp)
pileup_pdf <- file.path(pileup_out_dir, "allelic_pileups_all_snps.pdf")

pdf(pileup_pdf, width = 11, height = 8.5)
for (i in seq_along(snps_to_plot)) {
  slot <- (i - 1L) %% per_page
  if (slot == 0L) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = 0.97,
      height = 0.97,
      layout = grid::grid.layout(grid_nrow, grid_ncol)
    ))
  }
  grid::pushViewport(grid::viewport(
    layout.pos.row = slot %/% grid_ncol + 1L,
    layout.pos.col = slot %% grid_ncol + 1L
  ))
  tryCatch(
    plot_allelic_pileup(snps_to_plot[i]),
    error = function(e) {
      message("Plot failed for ", snps_to_plot[i], ": ", conditionMessage(e))
    }
  )
  grid::popViewport()
  if (slot == per_page - 1L || i == length(snps_to_plot)) {
    grid::popViewport()
  }
}
dev.off()
