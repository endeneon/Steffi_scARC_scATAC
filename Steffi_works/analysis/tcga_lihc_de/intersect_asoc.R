#!/usr/bin/env Rscript
# Intersect saved TCGA-LIHC DE table with ASoC genes (fixes merge key case).
suppressMessages(library(data.table))
here   <- "/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main/Steffi_works"
outdir <- file.path(here, "analysis", "tcga_lihc_de")

res <- fread(file.path(outdir, "TCGA-LIHC_tumor_vs_normal_DE_full.tsv"))
setnames(res, "Ensembl_base", "ENSEMBL_base")

asoc <- fread(file.path(here, "sig_ASoC_by_celltype", "sig_ASoC_in_Hepatocyte_annotated.tsv"))
asoc[, ENSEMBL_base := sub("\\..*$", "", ENSEMBL)]
asoc_genes <- unique(asoc[!is.na(ENSEMBL) & ENSEMBL != "NA" & ENSEMBL != "",
                          .(ENSEMBL_base, SYMBOL)])

m <- merge(asoc_genes, res, by = "ENSEMBL_base", all.x = TRUE)
m[, direction := fifelse(is.na(logFC), "not_tested",
                  fifelse(FDR < 0.05 & logFC >=  0.2, "Up",
                  fifelse(FDR < 0.05 & logFC <= -0.2, "Down", "ns")))]
setorder(m, -logFC)
fwrite(m, file.path(outdir, "ASoC_genes_TCGA-LIHC_DE.tsv"), sep = "\t")

cat("\n== ASoC-gene direction summary (TCGA-LIHC, |log2FC|>=0.2 & FDR<0.05) ==\n")
print(m[, .N, by = direction][order(-N)])

# Panel of genes cited in Section 3 (strong TF-disruption set)
panel <- c("FZD7","USP22","NDRG1","HNF4G","CAPN2","LAMC1","TPD52","NEK6","CHEK1","YWHAZ","CPT1A",
           "TP53INP1","KLF10","MTSS1","CES1",
           "INSIG1","DUSP10","GSDME","HES1","GPAM","MAP3K2","HDAC11","LRP5")
cat("\n== Section-3 panel genes: measured TCGA-LIHC tumor-vs-normal ==\n")
p <- m[SYMBOL %in% panel, .(SYMBOL, ENSEMBL_base, log2FC = round(logFC,3),
                            logCPM = round(logCPM,2), FDR = signif(FDR,3), direction)]
setorder(p, -log2FC)
print(p, nrow = 100)
