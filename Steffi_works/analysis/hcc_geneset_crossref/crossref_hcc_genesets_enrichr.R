#!/usr/bin/env Rscript
# Cross-reference the hepatocyte ASoC gene list against curated HCC gene sets,
# using ENRICHR for the DisGeNET half (tokenless; no DisGeNET API key / rate limit).
# Sources:
#   - COSMIC Cancer Gene Census, HCC-mapped: Open Targets Platform GraphQL
#     (`cancer_gene_census` datasource for hepatocellular carcinoma, MONDO_0007256).
#   - DisGeNET "Carcinoma, Hepatocellular" (C2239176): Enrichr enrichment API
#     (addList -> enrich, backgroundType=DisGeNET). The overlap of the ASoC list
#     with the Enrichr "Liver carcinoma" (== C2239176) term gives membership plus
#     a set-level enrichment adjusted p-value. NOTE: Enrichr's DisGeNET is a fixed,
#     text-mining-inclusive snapshot and does NOT carry the native GDA score / EI
#     (use crossref_hcc_genesets_disgenet.R for those).
# Input: analysis/tcga_lihc_de/ASoC_genes_TCGA-LIHC_DE.tsv (230 ASoC genes + DE).

suppressMessages({
  library(data.table)
  library(jsonlite)
  library(curl)
})

here    <- Sys.getenv("STEFFI_WORKS_DIR",
                      unset = "/home/szhang37/CAB_workspace/pulled_git_repos/Multiome_main/Steffi_works")
outdir  <- file.path(here, "analysis", "hcc_geneset_crossref")
datadir <- file.path(outdir, "data")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

## ---- ASoC gene list (+ measured TCGA-LIHC DE direction) ----
asoc <- fread(file.path(here, "analysis", "tcga_lihc_de", "ASoC_genes_TCGA-LIHC_DE.tsv"),
              sep = "\t", header = TRUE)
setnames(asoc, "ENSEMBL_base", "ensembl")
asoc <- unique(asoc[!is.na(ensembl) & ensembl != "", .(
  ensembl, SYMBOL,
  TCGA_LIHC_log2FC = round(logFC, 3),
  TCGA_LIHC_dir    = direction
)])
message("ASoC genes with Ensembl ID: ", nrow(asoc))

## ---- (1) COSMIC Cancer Gene Census (HCC) via Open Targets ----
ot_cache <- file.path(datadir, "opentargets_hcc_assoc.json")
if (!file.exists(ot_cache)) {
  message("Querying Open Targets Platform for HCC (MONDO_0007256)")
  gql <- sprintf(
    '{ disease(efoId:"MONDO_0007256"){ name associatedTargets(Bs:[%s], page:{index:0,size:%d}){ count rows{ target{ approvedSymbol id } score datasourceScores{ id score } } } } }',
    paste0('"', asoc$ensembl, '"', collapse = ","),
    nrow(asoc))
  body <- jsonlite::toJSON(list(query = gql), auto_unbox = TRUE)
  h <- curl::new_handle()
  curl::handle_setheaders(h, "Content-Type" = "application/json")
  curl::handle_setopt(h, copypostfields = body)
  resp <- curl::curl_fetch_memory("https://api.platform.opentargets.org/api/v4/graphql", handle = h)
  if (resp$status_code != 200)
    stop("Open Targets HTTP ", resp$status_code, ": ", rawToChar(resp$content))
  writeLines(rawToChar(resp$content), ot_cache)
}
ot <- jsonlite::fromJSON(ot_cache, simplifyVector = FALSE)
ot_rows <- ot$data$disease$associatedTargets$rows
pull_ds <- function(ds_list, id) {
  v <- vapply(ds_list, function(s) if (identical(s$id, id)) s$score else NA_real_, numeric(1))
  if (all(is.na(v))) 0 else max(v, na.rm = TRUE)
}
ot_dt <- rbindlist(lapply(ot_rows, function(r) data.table(
  ensembl        = r$target$id,
  OT_HCC_overall = round(r$score, 4),
  CGC_HCC_score  = round(pull_ds(r$datasourceScores, "cancer_gene_census"), 4),
  OT_clinvar     = round(max(pull_ds(r$datasourceScores, "eva"),
                             pull_ds(r$datasourceScores, "eva_somatic")), 4),
  OT_intogen     = round(pull_ds(r$datasourceScores, "intogen"), 4),
  OT_europepmc   = round(pull_ds(r$datasourceScores, "europepmc"), 4)
)), fill = TRUE)
message("ASoC genes with any Open Targets HCC association: ", nrow(ot_dt))

## ---- (2) DisGeNET "Carcinoma, Hepatocellular" (C2239176) via Enrichr ----
enr_cache <- file.path(datadir, "enrichr_disgenet_enrich.json")
if (!file.exists(enr_cache)) {
  syms <- unique(asoc$SYMBOL[!is.na(asoc$SYMBOL) & asoc$SYMBOL != ""])
  message("Enrichr addList: ", length(syms), " ASoC genes")
  h <- curl::new_handle()
  curl::handle_setform(h, list = paste(syms, collapse = "\n"),
                       description = "ASoC hepatocyte genes")
  r <- curl::curl_fetch_memory("https://maayanlab.cloud/Enrichr/addList", handle = h)
  if (r$status_code != 200) stop("Enrichr addList HTTP ", r$status_code, ": ", rawToChar(r$content))
  uid <- jsonlite::fromJSON(rawToChar(r$content))$userListId
  message("Enrichr userListId: ", uid, "; enriching against DisGeNET")
  eu <- paste0("https://maayanlab.cloud/Enrichr/enrich?userListId=", uid, "&backgroundType=DisGeNET")
  enr_txt <- NULL
  for (attempt in seq_len(5)) {
    r2 <- curl::curl_fetch_memory(eu)
    if (r2$status_code == 200) {
      tmp <- rawToChar(r2$content)
      if (!is.null(jsonlite::fromJSON(tmp, simplifyVector = FALSE)$DisGeNET)) { enr_txt <- tmp; break }
    }
    Sys.sleep(2)
  }
  if (is.null(enr_txt)) stop("Enrichr enrich (DisGeNET) returned no results after retries")
  writeLines(enr_txt, enr_cache)
}
enr <- jsonlite::fromJSON(enr_cache, simplifyVector = FALSE)$DisGeNET
# Enrichr row = [rank, term, p, oddsRatio, combinedScore, overlapGenes, adjP, ...]
liver <- Filter(function(z) tolower(z[[2]]) == "liver carcinoma", enr)  # C2239176
if (length(liver) == 0) {
  dg_overlap <- character(0); dg_p <- NA_real_; dg_padj <- NA_real_
  message("WARNING: no 'Liver carcinoma' term overlap in Enrichr DisGeNET")
} else {
  dg_overlap <- toupper(unlist(liver[[1]][[6]]))
  dg_p       <- as.numeric(liver[[1]][[3]])
  dg_padj    <- as.numeric(liver[[1]][[7]])
  message(sprintf("Enrichr DisGeNET 'Liver carcinoma': %d ASoC genes overlap (p=%.2g, padj=%.2g)",
                  length(dg_overlap), dg_p, dg_padj))
}

## ---- Merge + classify ----
x <- merge(asoc, ot_dt, by = "ensembl", all.x = TRUE)
for (col in c("OT_HCC_overall", "CGC_HCC_score", "OT_clinvar", "OT_intogen", "OT_europepmc"))
  x[is.na(get(col)), (col) := 0]
x[, COSMIC_CGC_HCC := fifelse(CGC_HCC_score > 0, "yes", "no")]
x[, DisGeNET_HCC   := fifelse(toupper(SYMBOL) %in% dg_overlap, "yes", "no")]
x[, DisGeNET_enrich_padj := fifelse(DisGeNET_HCC == "yes", round(dg_padj, 4), NA_real_)]
x[, n_curated_sources := (COSMIC_CGC_HCC == "yes") + (DisGeNET_HCC == "yes")]
x[, curated_HCC_hit := fifelse(n_curated_sources > 0, "yes", "no")]

setcolorder(x, c("SYMBOL", "ensembl", "TCGA_LIHC_dir", "TCGA_LIHC_log2FC",
                 "COSMIC_CGC_HCC", "CGC_HCC_score", "DisGeNET_HCC", "DisGeNET_enrich_padj",
                 "n_curated_sources", "curated_HCC_hit",
                 "OT_HCC_overall", "OT_clinvar", "OT_intogen", "OT_europepmc"))
setorder(x, -n_curated_sources, -OT_HCC_overall, -CGC_HCC_score)

out_tsv <- file.path(outdir, "ASoC_genes_HCC_curated_crossref_enrichr.tsv")
fwrite(x, out_tsv, sep = "\t")

## ---- Summary ----
message("\n===== Curated HCC gene-set cross-reference (Enrichr DisGeNET) =====")
message("ASoC genes total (Ensembl):     ", nrow(x))
message("In COSMIC CGC (HCC-mapped):      ", sum(x$COSMIC_CGC_HCC == "yes"))
message("In DisGeNET C2239176 (Enrichr):  ", sum(x$DisGeNET_HCC == "yes"))
message("In BOTH curated sets:            ", sum(x$n_curated_sources == 2))
message("In EITHER curated set:           ", sum(x$curated_HCC_hit == "yes"))
message("\nGenes in >=1 curated HCC set (by curated sources, then OT score):")
print(x[curated_HCC_hit == "yes",
        .(SYMBOL, TCGA_LIHC_dir, COSMIC_CGC_HCC, DisGeNET_HCC, OT_HCC_overall)])
message("\nWritten: ", out_tsv)
