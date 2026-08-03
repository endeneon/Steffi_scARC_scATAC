#!/usr/bin/env Rscript
# Join the per-gene curated-HCC cross-reference (Enrichr variant) back onto the
# per-SNP hepatocyte ASoC table, adding the curated-HCC columns to every variant
# row (matched on base Ensembl gene ID). Rows whose gene is absent from the
# cross-reference (e.g. ENSEMBL == NA) get NA in the added columns.

suppressMessages(library(data.table))

here   <- Sys.getenv("STEFFI_WORKS_DIR",
                     unset = "/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main/Steffi_works")
in_tsv  <- file.path(here, "sig_ASoC_by_celltype", "sig_ASoC_in_Hepatocyte_annotated.tsv")
xref_tsv<- file.path(here, "analysis", "hcc_geneset_crossref", "ASoC_genes_HCC_curated_crossref_enrichr.tsv")
out_tsv <- file.path(here, "sig_ASoC_by_celltype", "sig_ASoC_in_Hepatocyte_annotated_HCC_crossref.tsv")

snps <- fread(in_tsv, sep = "\t", header = TRUE)
xref <- fread(xref_tsv, sep = "\t", header = TRUE)

# columns to attach (everything from the cross-ref except the join key + duplicate SYMBOL)
add_cols <- setdiff(names(xref), c("ensembl", "SYMBOL"))
xref_join <- xref[, c("ensembl", add_cols), with = FALSE]

snps[, .row_order := .I]                       # preserve original variant order
m <- merge(snps, xref_join, by.x = "ENSEMBL", by.y = "ensembl", all.x = TRUE, sort = FALSE)
setorder(m, .row_order)
m[, .row_order := NULL]
setcolorder(m, c(setdiff(names(snps), ".row_order"), add_cols))  # keep original cols first

fwrite(m, out_tsv, sep = "\t")

message("Input variant rows:        ", nrow(snps))
message("Rows matched to a gene:    ", sum(!is.na(m$curated_HCC_hit)))
message("Variant rows in a curated HCC set: ", sum(m$curated_HCC_hit == "yes", na.rm = TRUE),
        " (COSMIC CGC: ", sum(m$COSMIC_CGC_HCC == "yes", na.rm = TRUE),
        ", DisGeNET: ", sum(m$DisGeNET_HCC == "yes", na.rm = TRUE), ")")
message("Columns added: ", paste(add_cols, collapse = ", "))
message("Written: ", out_tsv)
