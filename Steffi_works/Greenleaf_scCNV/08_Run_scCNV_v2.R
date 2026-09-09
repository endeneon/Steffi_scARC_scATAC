#Estimating Copy Number Variation in scATAC-seq
#05/02/19
#Cite Satpathy*, Granja*, et al.
#Massively parallel single-cell chromatin landscapes of human immune
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(future.apply)

library(qs2)
set.seed(1)

"%ni%" <- Negate("%in%")

#--------------------------------------------
# Parallel backend
#--------------------------------------------
# Job allocation: 16 cores + 300 GB RAM.
# Prefer forking (multicore) on Linux/LSF so large shared objects (`windows`,
# `fragments`, `counts`, `countSummary`, the genome) are inherited copy-on-write
# by workers instead of being serialized to each one.
# Fall back to multisession only in interactive sessions where forking is unsafe.
n_cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", unset = NA))
if (is.na(n_cores) || n_cores < 1L) {
	n_cores <- 16L
}
n_cores <- min(n_cores, 16L)

session_plan <- if (interactive()) "multisession" else "multicore"
future::plan(session_plan, workers = n_cores)

# 300 GB across 16 workers -> allow generous per-worker globals (used only if a
# socket backend is in play; forked workers inherit memory copy-on-write).
options(future.globals.maxSize = 16 * 1024^3)
options(future.rng.onMisuse = "ignore")
message(sprintf("Parallel backend: %s with %d workers", session_plan, n_cores))

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

#----------------------------
# Get Inputs
#----------------------------
blacklist <- import.bed("data/hg38.blacklist.bed")
windows <- makeWindows(
	genome = BSgenome.Hsapiens.UCSC.hg38,
	blacklist = blacklist
)
cnaObj <- scCNA(
	windows,
	readRDS("data/PBMC_10x-Sub25M-fragments.gr.rds"),
	neighbors = 100,
	LFC = 1.5,
	FDR = 0.1,
	force = FALSE,
	remove = c("chrM", "chrX", "chrY")
)
saveRDS(cnaObj, "results/PBMC-CNV_LFC_GC.rds")
