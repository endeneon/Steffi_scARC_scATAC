#Merged scATAC-seq pipeline: Cell Filtering (step 01) + scCNV (step 08)
#Executes 01_Filter_Cells_v2.R logic first, then 08_Run_scCNV_v2.R logic.
#The filtered `fragments` GRanges is passed directly in memory to scCNA(),
#skipping the lengthy saveRDS/readRDS round-trip between the two steps.
#Cite Satpathy*, Granja*, et al.
#Massively parallel single-cell chromatin landscapes of human immune
#cell development and intratumoral T cell exhaustion (2019)
#Original scripts created by Jeffrey Granja
{
  library(Matrix)
  library(SummarizedExperiment)
  library(matrixStats)
  library(readr)
  library(GenomicRanges)
  library(magrittr)
  library(edgeR)
  library(Seurat)
  library(rtracklayer) # import.bed()
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
  library(Rcpp)
  library(viridis)
  library(future)
  library(future.apply)

  library(qs2)
}
set.seed(1)

"%ni%" <- Negate("%in%")

#--------------------------------------------
# Parallel backend
#--------------------------------------------
# Per-sample job allocation: 20 cores + 300 GB RAM (one LSF job per sample; see
# submit_scCNV_lsf.sh). Prefer forking (multicore) on Linux/LSF so the large
# shared objects (`fragments`, `windows`, `counts`, `countSummary`, the genome)
# are inherited copy-on-write by workers instead of being serialized to each one.
# Fall back to multisession only in interactive sessions where forking is unsafe.
n_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = NA))
if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- 20L
}
n_cores <- min(n_cores, 20L)

session_plan <- if (interactive()) "multisession" else "multicore"
future::plan(session_plan, workers = n_cores)

# future.apply size-checks the globals it captures on EVERY backend (not just
# socket clusters), so the cap must accommodate the largest per-future globals
# even when forking. Here `fragments` (~29 GiB) is captured three times by the
# same future expression -- directly, inside the FUN closure, and inside the
# enclosing environment of insertionProfileSingles_helper -- which the size
# check sums to ~86 GiB. Under multicore (forking) these all point at the same
# object and are inherited copy-on-write, so nothing is actually duplicated;
# raise the cap generously (well under the 300 GB LSF allocation) so the
# spurious sum no longer aborts the run.
options(future.globals.maxSize = 200 * 1024^3)
options(future.rng.onMisuse = "ignore")
message(sprintf("Parallel backend: %s with %d workers", session_plan, n_cores))

#============================================================================
# Step 01 functions: Cell Filtering (TSS enrichment)
#============================================================================

sourceCpp(
  code = '
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
    if(x1.size() != y1.size()){
      stop("width must equal size!");
    }
    IntegerVector x = clone(x1);
    IntegerVector y = clone(y1);
    int n = x.size();
    IntegerVector rx = seq(xmin,xmax);
    IntegerVector ry = seq(ymin,ymax);
    IntegerMatrix mat( ry.size() , rx.size() );
    int xi,yi;
    for(int i = 0; i < n; i++){
      xi = (x[i] - xmin);
      yi = (y[i] - ymin);
      if(yi >= 0 && yi < ry.size()){
        if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
        }
      }
    }
    return mat;
  }'
)

insertionProfileSingles <- function(
  feature,
  fragments,
  by = "RG",
  getInsertions = TRUE,
  fix = "center",
  flank = 2000,
  norm = 100,
  smooth = 51,
  range = 100,
  batchSize = 100
) {
  insertionProfileSingles_helper <- function(
    feature,
    fragments,
    by = "RG",
    getInsertions = TRUE,
    fix = "center",
    flank = 2000,
    norm = 100,
    smooth = 51,
    range = 100,
    batchSize = 100
  ) {
    #Convert To Insertion Sites
    if (getInsertions) {
      insertions <- c(
        GRanges(
          seqnames = seqnames(fragments),
          ranges = IRanges(start(fragments), start(fragments)),
          RG = mcols(fragments)[, by]
        ),
        GRanges(
          seqnames = seqnames(fragments),
          ranges = IRanges(end(fragments), end(fragments)),
          RG = mcols(fragments)[, by]
        )
      )
      by <- "RG"
    } else {
      insertions <- fragments
    }
    remove(fragments)
    gc()

    #center the feature
    center <- unique(resize(
      feature,
      width = 1,
      fix = fix,
      ignore.strand = FALSE
    ))

    #get overlaps between the feature and insertions only up to flank bp
    overlap <- DataFrame(findOverlaps(
      query = center,
      subject = insertions,
      maxgap = flank,
      ignore.strand = TRUE
    ))
    overlap$strand <- strand(center)[overlap[, 1]]
    overlap$name <- mcols(insertions)[overlap[, 2], by]
    overlap <- transform(overlap, id = match(name, unique(name)))
    ids <- length(unique(overlap$name))

    #distance
    overlap$dist <- NA
    minus <- which(overlap$strand == "-")
    other <- which(overlap$strand != "-")
    overlap$dist[minus] <- start(center[overlap[minus, 1]]) -
      start(insertions[overlap[minus, 2]])
    overlap$dist[other] <- start(insertions[overlap[other, 2]]) -
      start(center[overlap[other, 1]])

    #Insertion Mat
    profile_mat <- tabulate2dCpp(
      x1 = overlap$id,
      y1 = overlap$dist,
      xmin = 1,
      xmax = ids,
      ymin = -flank,
      ymax = flank
    )
    colnames(profile_mat) <- unique(overlap$name)
    profile <- rowSums(profile_mat)

    #normalize
    profile_mat_norm <- apply(profile_mat, 2, function(x) {
      x / max(mean(x[c(1:norm, (flank * 2 - norm + 1):(flank * 2 + 1))]), 0.5)
    }) #Handles low depth cells
    profile_norm <- profile /
      mean(profile[c(1:norm, (flank * 2 - norm + 1):(flank * 2 + 1))])

    #smooth
    profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) {
      zoo::rollmean(x, smooth, fill = 1)
    })
    profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)

    #enrichment
    max_finite <- function(x) {
      suppressWarnings(max(x[is.finite(x)], na.rm = TRUE))
    }
    e_mat <- apply(profile_mat_norm_smooth, 2, function(x) {
      max_finite(x[(flank - range):(flank + range)])
    })
    names(e_mat) <- colnames(profile_mat_norm_smooth)
    e <- max_finite(profile_norm_smooth[(flank - range):(flank + range)])

    #Summary
    df_mat <- data.frame(
      enrichment = e_mat,
      insertions = as.vector(table(mcols(insertions)[, by])[names(e_mat)]),
      insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
    )
    df_sum <- data.frame(
      bp = (-flank):flank,
      profile = profile,
      norm_profile = profile_norm,
      smooth_norm_profile = profile_norm_smooth,
      enrichment = e
    )
    rownames(df_sum) <- NULL

    return(list(
      df = df_sum,
      dfall = df_mat,
      profileMat = profile_mat_norm,
      profileMatSmooth = profile_mat_norm_smooth
    ))
  }

  uniqueTags <- as.character(unique(mcols(fragments)[, by]))
  splitTags <- split(uniqueTags, ceiling(seq_along(uniqueTags) / batchSize))

  # Each batch is fully independent (its own findOverlaps + tabulate2dCpp), so
  # distribute batches across workers. Forked workers inherit `fragments`,
  # `feature`, the helper, and the compiled tabulate2dCpp() copy-on-write.
  # (Progress bar dropped: it cannot update meaningfully across workers.)
  batchTSS <- future.apply::future_lapply(
    seq_along(splitTags),
    function(x) {
      insertionProfileSingles_helper(
        feature = feature,
        fragments = fragments[which(
          mcols(fragments)[, by] %in% splitTags[[x]]
        )],
        by = by,
        getInsertions = getInsertions,
        fix = fix,
        flank = flank,
        norm = norm,
        smooth = smooth,
        range = range
      )
    },
    future.seed = TRUE
  )
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind", .)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind", .)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>%
    Reduce("cbind", .)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>%
    Reduce("cbind", .)
  return(list(
    df = df,
    dfall = dfall,
    profileMat = profileMat,
    profileMatSmooth = profileMatSmooth
  ))
}

#============================================================================
# Step 08 functions: scCNV estimation
#============================================================================

countInsertions <- function(query, fragments, by = "RG") {
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(
      seqnames = seqnames(fragments),
      ranges = IRanges(start(fragments), start(fragments)),
      RG = mcols(fragments)[, by]
    ),
    GRanges(
      seqnames = seqnames(fragments),
      ranges = IRanges(end(fragments), end(fragments)),
      RG = mcols(fragments)[, by]
    )
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(
    query,
    inserts,
    ignore.strand = TRUE,
    maxgap = -1L,
    minoverlap = 0L,
    type = "any"
  ))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1],
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)),
    dims = c(length(query), length(unique(overlapDF$name)))
  )
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

makeWindows <- function(
  genome,
  blacklist,
  windowSize = 10e6,
  slidingSize = 2e6
) {
  chromSizes <- GRanges(
    names(seqlengths(genome)),
    IRanges(1, seqlengths(genome))
  )
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(
    chromSizes,
    pruning.mode = "coarse"
  )
  windows <- slidingWindows(
    x = chromSizes,
    width = windowSize,
    step = slidingSize
  ) %>%
    unlist %>%
    .[which(width(.) == windowSize), ]
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- start(windows)
  mcols(windows)$wEnd <- end(windows)
  message("Subtracting Blacklist...")
  # Independent per-window setdiff -> distribute across workers (fork = shared
  # `windows`/`blacklist` copy-on-write). Progress messages dropped in parallel.
  windowsBL <- future.apply::future_lapply(
    seq_along(windows),
    function(x) {
      gr <- GenomicRanges::setdiff(windows[x, ], blacklist)
      mcols(gr) <- mcols(windows[x, ])
      return(gr)
    },
    future.seed = TRUE
  )
  names(windowsBL) <- paste0("w", seq_along(windowsBL))
  windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- future.apply::future_lapply(
    seq_along(windowSplit),
    function(x) {
      chrSeq <- Biostrings::getSeq(
        genome,
        chromSizes[which(seqnames(chromSizes) == names(windowSplit)[x])]
      )
      grx <- windowSplit[[x]]
      aFreq <- alphabetFrequency(Biostrings::Views(
        chrSeq[[1]],
        ranges(grx)
      ))
      mcols(grx)$GC <- rowSums(aFreq[, c("G", "C")]) / rowSums(aFreq)
      mcols(grx)$AT <- rowSums(aFreq[, c("A", "T")]) / rowSums(aFreq)
      return(grx)
    },
    future.seed = TRUE
  ) %>%
    GenomicRangesList %>%
    unlist %>%
    sortSeqlevels %>%
    sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  windowNuc
}

scCNA <- function(
  windows,
  fragments,
  neighbors = 100,
  LFC = 1.5,
  FDR = 0.1,
  force = FALSE,
  remove = c("chrM", "chrX", "chrY")
) {
  #Keep only regions in filtered chromosomes
  windows <- GenomeInfoDb::keepStandardChromosomes(
    windows,
    pruning.mode = "coarse"
  )
  fragments <- GenomeInfoDb::keepStandardChromosomes(
    fragments,
    pruning.mode = "coarse"
  )
  windows <- windows[seqnames(windows) %ni% remove]
  fragments <- fragments[seqnames(fragments) %ni% remove]

  #Count Insertions in windows
  message("Getting Counts...")
  counts <- countInsertions(windows, fragments, by = "RG")[[1]]
  message("Summarizing...")
  # Each unique window name is summarized independently; parallelize and then
  # reassemble in order (future_lapply preserves ordering).
  uniqueNames <- unique(mcols(windows)$name)
  summaryList <- future.apply::future_lapply(
    seq_along(uniqueNames),
    function(x) {
      idx <- which(mcols(windows)$name == uniqueNames[x])
      wx <- windows[idx, ]
      wo <- GRanges(
        mcols(wx)$wSeq,
        ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd)
      )[1, ]
      mcols(wo)$name <- mcols(wx)$name[1]
      mcols(wo)$effectiveLength <- sum(width(wx))
      mcols(wo)$percentEffectiveLength <- 100 * sum(width(wx)) / width(wo)
      mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx)) / width(wo)
      mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx)) / width(wo)
      mcols(wo)$N <- sum(mcols(wx)$N * width(wx)) / width(wo)
      list(wo = wo, counts = Matrix::colSums(counts[idx, , drop = FALSE]))
    },
    future.seed = TRUE
  )
  windowSummary <- unlist(GenomicRangesList(lapply(summaryList, `[[`, "wo")))
  countSummary <- do.call(rbind, lapply(summaryList, `[[`, "counts"))

  #Keep only regions with less than 0.1% N
  keep <- which(windowSummary$N < 0.001)
  windowSummary <- windowSummary[keep, ]
  countSummary <- countSummary[keep, ]

  #Now determine the nearest neighbors by GC content
  message("Computing Background...")
  # Each window's background is computed independently from its GC nearest
  # neighbors -> embarrassingly parallel over rows. Compute per-row vectors in
  # parallel, then assemble the result matrices (order preserved).
  gcVec <- windowSummary$GC
  bgList <- future.apply::future_lapply(
    seq_len(nrow(countSummary)),
    function(x) {
      #Get Nearest Indices
      idxNN <- head(
        order(abs(gcVec[x] - gcVec)),
        neighbors + 1
      )
      idxNN <- idxNN[idxNN %ni% x]
      #Background
      bdgMeanx <- colMeans(countSummary[idxNN, ])
      if (any(bdgMeanx == 0)) {
        if (force) {
          message(
            "Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins"
          )
        } else {
          stop("Background Mean = 0!")
        }
      }
      bdgSdx <- matrixStats::colSds(countSummary[idxNN, ])
      log2FCx <- log2((countSummary[x, ] + 1e-5) / (bdgMeanx + 1e-5))
      zx <- (countSummary[x, ] - bdgMeanx) / bdgSdx
      pvalx <- 2 * pnorm(-abs(zx))
      list(
        bdgMean = bdgMeanx,
        bdgSd = bdgSdx,
        log2FC = log2FCx,
        z = zx,
        pval = pvalx
      )
    },
    future.seed = TRUE
  )
  bdgMean <- do.call(rbind, lapply(bgList, `[[`, "bdgMean"))
  bdgSd <- do.call(rbind, lapply(bgList, `[[`, "bdgSd"))
  log2FC <- do.call(rbind, lapply(bgList, `[[`, "log2FC"))
  z <- do.call(rbind, lapply(bgList, `[[`, "z"))
  pval <- do.call(rbind, lapply(bgList, `[[`, "pval"))
  padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
  CNA <- matrix(0, nrow = nrow(countSummary), ncol = ncol(countSummary))
  CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

  se <- SummarizedExperiment(
    assays = SimpleList(
      CNA = CNA,
      counts = countSummary,
      log2FC = log2FC,
      padj = padj,
      pval = pval,
      z = z,
      bdgMean = bdgMean,
      bdgSd = bdgSd
    ),
    rowRanges = windowSummary
  )
  colnames(se) <- colnames(counts)

  return(se)
}

#============================================================================
# Command-line arguments (ONE sample per invocation)
#   1: fragment_file  - path to a CellRanger fragments .tsv.gz
#   2: sample_label   - unique label, e.g. "multiome_1t" or "atac_24t"
#                       (used to tag barcodes and name all output files)
#   3: output_dir     - directory for results (optional; default "results")
#   4: blacklist_bed  - blacklist BED (optional; default hg38 v2 blacklist)
# Submitted one-per-sample by submit_scCNV_lsf.sh.
#============================================================================
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) < 2L) {
  stop(paste(
    "Usage: Rscript Run_scCNV_pipeline.R",
    "<fragment_file> <sample_label> [output_dir] [blacklist_bed]"
  ))
}
file_fragments <- cli_args[[1]]
name <- cli_args[[2]]
output_dir <- if (length(cli_args) >= 3L && nzchar(cli_args[[3]])) {
  cli_args[[3]]
} else {
  "results"
}
blacklist_bed <- if (length(cli_args) >= 4L && nzchar(cli_args[[4]])) {
  cli_args[[4]]
} else {
  paste0(
    "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/",
    "szhang37/Databases/Genome/Blacklist/lists/hg38-blacklist.v2.bed"
  )
}

if (!file.exists(file_fragments)) {
  stop(sprintf("Fragment file not found: %s", file_fragments))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
message(sprintf("Processing sample '%s' from %s", name, file_fragments))

#============================================================================
# Step 01 driver: Cell Filtering
#============================================================================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
minFrags <- 100
filterFrags <- 1000
filterTSS <- 8

#-----------------
# Reading Fragment Files
#-----------------
# CellRanger / cellranger-arc fragment files carry a '#' comment header block;
# skip it with comment = "#". Columns: chr, start (0-based), end, barcode, count.
message("Reading in fragment files...")
fragments <- data.frame(readr::read_tsv(
  file_fragments,
  col_names = FALSE,
  comment = "#"
))

fragments <- GRanges(
  seqnames = fragments[, 1],
  IRanges(fragments[, 2] + 1, fragments[, 3]),
  RG = fragments[, 4],
  N = fragments[, 5]
)

message("Filtering Lowly Represented Cells...")
tabRG <- table(fragments$RG)
keep <- names(tabRG)[which(tabRG >= minFrags)]
fragments <- fragments[fragments$RG %in% keep, ]
fragments <- sort(sortSeqlevels(fragments))

#-----------------
# TSS Profile
#-----------------
feature <- txdb %>%
  transcripts(.) %>%
  resize(., width = 1, fix = "start") %>%
  unique
# Smaller batches -> more independent tasks to keep all 20 workers busy and to
# shrink each worker's per-batch fragment subset (lower peak memory).
tssProfile <- insertionProfileSingles(
  feature = feature,
  fragments = fragments,
  getInsertions = TRUE,
  batchSize = 100
)
tssSingles <- tssProfile$dfall
tssSingles$uniqueFrags <- 0
tssSingles[names(tabRG), "uniqueFrags"] <- tabRG
tssSingles$cellCall <- 0
tssSingles$cellCall[
  tssSingles$uniqueFrags >= filterFrags & tssSingles$enrichment >= filterTSS
] <- 1

#-----------------
# Plot Stats (named by sample label)
#-----------------
tssSingles <- tssSingles[complete.cases(tssSingles), ]
nPass <- sum(tssSingles$cellCall == 1)
nTotal <- sum(tssSingles$uniqueFrags >= filterFrags)

qc_pdf <- file.path(output_dir, paste0(name, "_Filter-Cells.pdf"))
qc_txt <- file.path(output_dir, paste0(name, "_Filter-Cells.txt"))

pdf(qc_pdf)
print(
  ggplot(
    tssSingles[tssSingles$uniqueFrags > 500, ],
    aes(x = log10(uniqueFrags), y = enrichment)
  ) +
    geom_hex(bins = 100) +
    theme_bw() +
    scale_fill_viridis() +
    xlab("log10 Unique Fragments") +
    ylab("TSS Enrichment") +
    geom_hline(yintercept = filterTSS, lty = "dashed") +
    geom_vline(xintercept = log10(filterFrags), lty = "dashed") +
    ggtitle(sprintf(
      "%s Pass Rate : %s of %s (%s)",
      name,
      nPass,
      nTotal,
      round(100 * nPass / nTotal, 2)
    ))
)
dev.off()

write.table(tssSingles, qc_txt)

#-----------------
# Filter fragments to called cells and tag with the sample label.
# Kept in memory and handed directly to scCNA() below (no saveRDS/readRDS).
#-----------------
fragments <- fragments[
  mcols(fragments)$RG %in% rownames(tssSingles)[tssSingles$cellCall == 1]
]
fragments$RG <- paste0(name, "#", fragments$RG)

#============================================================================
# Step 08 driver: scCNV estimation (uses the in-memory `fragments` above)
#============================================================================
blacklist <- import.bed(blacklist_bed)
windows <- makeWindows(
  genome = BSgenome.Hsapiens.UCSC.hg38,
  blacklist = blacklist
)
cnaObj <- scCNA(
  windows,
  fragments,
  neighbors = 100,
  LFC = 1.5,
  FDR = 0.1,
  force = FALSE,
  remove = c("chrM", "chrX", "chrY")
)
out_file <- file.path(output_dir, paste0(name, "_CNV_LFC_GC.qs2"))
qs_save(cnaObj, out_file)
message(sprintf("Done. Wrote %s", out_file))
