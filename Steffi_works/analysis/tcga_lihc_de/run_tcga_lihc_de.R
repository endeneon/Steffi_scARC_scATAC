#!/usr/bin/env Rscript
# TCGA-LIHC tumor-vs-normal differential expression for ASoC hepatocyte genes.
# Data: UCSC Xena GDC hub STAR counts (log2(count+1)) + gencode v36 probemap.
# DE: edgeR quasi-likelihood F-test (tumor 01 vs solid-normal 11).
# Output: full DE table + ASoC-gene intersection filtered by |log2FC|>=0.2 & FDR<0.05.

suppressMessages({
  library(data.table)
  library(edgeR)
})

set.seed(1)
here    <- "/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main/Steffi_works"
outdir  <- file.path(here, "analysis", "tcga_lihc_de")
datadir <- file.path(outdir, "data")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

counts_url   <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz"
probemap_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap"
counts_gz    <- file.path(datadir, "TCGA-LIHC.star_counts.tsv.gz")
probemap_f   <- file.path(datadir, "gencode.v36.gene.probemap")

dl <- function(url, dest) if (!file.exists(dest)) {
  message("Downloading ", basename(dest))
  download.file(url, dest, mode = "wb", quiet = TRUE)
}
dl(counts_url, counts_gz)
dl(probemap_url, probemap_f)

## ---- ASoC Ensembl IDs ----
asoc <- fread(file.path(here, "sig_ASoC_by_celltype", "sig_ASoC_in_Hepatocyte_annotated.tsv"),
              sep = "\t", header = TRUE)
asoc[, ENSEMBL_base := sub("\\..*$", "", ENSEMBL)]
asoc_genes <- unique(asoc[!is.na(ENSEMBL) & ENSEMBL != "NA" & ENSEMBL != "",
                          .(ENSEMBL_base, SYMBOL)])
message("ASoC annotated genes with Ensembl ID: ", nrow(asoc_genes))

## ---- Expression matrix (log2(count+1)) ----
mat <- fread(cmd = paste("zcat", shQuote(counts_gz)), sep = "\t", header = TRUE)
setnames(mat, 1, "Ensembl")
mat[, Ensembl_base := sub("\\..*$", "", Ensembl)]

samples <- setdiff(colnames(mat), c("Ensembl", "Ensembl_base"))
stype   <- substr(samples, 14, 15)                     # 01=tumor, 11=solid normal
grp     <- ifelse(stype == "11", "normal",
            ifelse(stype %in% c("01","02","03","04","05","06","07","08","09"), "tumor", NA))
keep_s  <- samples[!is.na(grp)]
grp     <- factor(grp[!is.na(grp)], levels = c("normal", "tumor"))
message("Samples: tumor=", sum(grp=="tumor"), "  normal=", sum(grp=="normal"))

## reverse log2(x+1) -> raw counts
cnt <- as.matrix(mat[, ..keep_s])
cnt <- round(2^cnt - 1)
cnt[cnt < 0] <- 0
rownames(cnt) <- mat$Ensembl_base
cnt <- cnt[!duplicated(rownames(cnt)), ]

## ---- edgeR DE ----
y  <- DGEList(counts = cnt, group = grp)
keep <- filterByExpr(y, group = grp)
y  <- y[keep, , keep.lib.sizes = FALSE]
y  <- calcNormFactors(y)
design <- model.matrix(~ grp)
y  <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)                        # tumor vs normal
res <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "Ensembl_base")

fwrite(res, file.path(outdir, "TCGA-LIHC_tumor_vs_normal_DE_full.tsv"), sep = "\t")

## ---- Intersect with ASoC genes ----
m <- merge(asoc_genes, res, by = "Ensembl_base", all.x = TRUE)
m[, direction := fifelse(is.na(logFC), "not_tested",
                  fifelse(FDR < 0.05 & logFC >=  0.2, "Up",
                  fifelse(FDR < 0.05 & logFC <= -0.2, "Down", "ns")))]
setorder(m, -logFC)
fwrite(m, file.path(outdir, "ASoC_genes_TCGA-LIHC_DE.tsv"), sep = "\t")

sig <- m[direction %in% c("Up","Down")][order(-logFC)]
message("\n===== ASoC genes significant in TCGA-LIHC (|log2FC|>=0.2, FDR<0.05) =====")
print(sig[, .(SYMBOL, Ensembl_base, logFC = round(logFC,3),
              logCPM = round(logCPM,2), FDR = signif(FDR,3), direction)])
message("\nSummary of ASoC-gene directions:")
print(m[, .N, by = direction][order(-N)])
