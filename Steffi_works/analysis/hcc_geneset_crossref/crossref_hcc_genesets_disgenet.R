#!/usr/bin/env Rscript
# Cross-reference the hepatocyte ASoC gene list against curated HCC gene sets,
# using the NATIVE DisGeNET API for the DisGeNET half (gives GDA score + Evidence
# Index, but needs an API key and is rate-limited; see also the tokenless Enrichr
# variant crossref_hcc_genesets_enrichr.R).
# Sources:
#   - COSMIC Cancer Gene Census, HCC-mapped: Open Targets Platform GraphQL
#     (`cancer_gene_census` datasource for hepatocellular carcinoma, MONDO_0007256).
#   - DisGeNET "Carcinoma, Hepatocellular" (C2239176): native DisGeNET REST API
#     (api.disgenet.com/api/v1/gda/summary) -> per-gene GDA `score` and `ei`.
#     Key is read at runtime from ~/.apikeys/Disgenet.key (override via
#     $DISGENET_KEYFILE, or supply the key via $DISGENET_API_KEY); never stored.
#     TRIAL accounts allow <=10 identifiers/request and rate-limit hard (retry-after
#     can be ~22 h); each batch is cached under data/disgenet_batches/ so reruns
#     resume. On an LSF batch node set DISGENET_MAX_WAIT high to wait out the reset.
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

## ---- (2) DisGeNET "Carcinoma, Hepatocellular" (C2239176) via native API ----
get_disgenet_key <- function() {
  k <- Sys.getenv("DISGENET_API_KEY")
  if (nzchar(k)) return(trimws(k))
  kf <- Sys.getenv("DISGENET_KEYFILE", unset = "~/.apikeys/Disgenet.key")
  kf <- path.expand(kf)
  if (!file.exists(kf)) stop("DisGeNET API key not found (set $DISGENET_API_KEY or $DISGENET_KEYFILE, or create ", kf, ")")
  trimws(readLines(kf, warn = FALSE)[1])
}
dg_cache    <- file.path(datadir, "disgenet_hcc_gda.json")
dg_batchdir <- file.path(datadir, "disgenet_batches")  # per-batch cache for resumable reruns
DG_MAX_WAIT <- as.numeric(Sys.getenv("DISGENET_MAX_WAIT", unset = "120"))
if (!file.exists(dg_cache)) {
  dir.create(dg_batchdir, showWarnings = FALSE, recursive = TRUE)
  key  <- get_disgenet_key()
  syms <- unique(asoc$SYMBOL[!is.na(asoc$SYMBOL) & asoc$SYMBOL != ""])
  batches <- split(syms, ceiling(seq_along(syms) / 10L))  # TRIAL profile: max 10 identifiers/request
  nb <- length(batches)
  message("Querying DisGeNET API for C2239176 GDAs: ", length(syms),
          " genes in ", nb, " batches (<=10 genes/request, TRIAL limit)")
  dg_get <- function(url, key, label, max_try = 4) {
    for (attempt in seq_len(max_try)) {
      hh <- curl::new_handle()
      curl::handle_setheaders(hh, "Authorization" = key, "Accept" = "application/json")
      r <- curl::curl_fetch_memory(url, handle = hh)
      if (r$status_code == 200) return(jsonlite::fromJSON(rawToChar(r$content), simplifyVector = FALSE))
      if (r$status_code == 429) {
        hdr  <- tryCatch(curl::parse_headers_list(r$headers), error = function(e) list())
        wait <- suppressWarnings(as.numeric(hdr[["x-rate-limit-retry-after-seconds"]]))
        if (is.na(wait)) wait <- suppressWarnings(as.numeric(hdr[["retry-after"]]))
        if (is.na(wait)) wait <- 20
        if (wait > DG_MAX_WAIT)
          stop("DisGeNET rate limit exhausted at ", label, ": retry-after ",
               round(wait), "s (~", round(wait/3600, 1), " h). Progress is cached in ",
               dg_batchdir, "; rerun this script after the TRIAL quota resets to resume.")
        message("    ", label, ": 429, waiting ", round(wait), "s then retry (", attempt, "/", max_try, ")")
        Sys.sleep(wait + 1); next
      }
      stop("DisGeNET HTTP ", r$status_code, " at ", label, ": ", rawToChar(r$content))
    }
    stop("DisGeNET still rate-limited at ", label, " after ", max_try, " attempts")
  }
  for (i in seq_along(batches)) {
    bf <- file.path(dg_batchdir, sprintf("batch_%02d.json", i))
    if (file.exists(bf)) {
      message(sprintf("[%d/%d] batch cached, skipping (%s)", i, nb, basename(bf)))
      next
    }
    label <- sprintf("batch %d/%d", i, nb)
    message(sprintf("[%d/%d] querying %d genes: %s", i, nb, length(batches[[i]]),
                    paste(batches[[i]], collapse = ",")))
    url <- paste0("https://api.disgenet.com/api/v1/gda/summary?disease=UMLS_C2239176",
                  "&gene_symbol=", paste(batches[[i]], collapse = ","))
    pj  <- dg_get(url, key, label)
    writeLines(jsonlite::toJSON(pj$payload, auto_unbox = TRUE), bf)
    hits <- vapply(pj$payload, function(r) r$symbolOfGene, character(1))
    message(sprintf("[%d/%d] done: %d/%d genes have a C2239176 GDA%s", i, nb,
                    length(hits), length(batches[[i]]),
                    if (length(hits)) paste0(" (", paste(hits, collapse = ","), ")") else ""))
    Sys.sleep(8)  # stay under the TRIAL request-rate cap
  }
  payload <- list()
  for (i in seq_along(batches))
    payload <- c(payload, jsonlite::fromJSON(file.path(dg_batchdir, sprintf("batch_%02d.json", i)),
                                             simplifyVector = FALSE))
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE), dg_cache)
  message("DisGeNET fetch complete: ", length(payload), " GDA records assembled")
}
dg_payload <- jsonlite::fromJSON(dg_cache, simplifyVector = FALSE)
dg_dt <- if (length(dg_payload) == 0) {
  data.table(SYMBOL = character(), DisGeNET_score = numeric(), DisGeNET_EI = numeric(), DisGeNET_nPMID = integer())
} else rbindlist(lapply(dg_payload, function(r) data.table(
  SYMBOL         = r$symbolOfGene,
  DisGeNET_score = as.numeric(r$score),
  DisGeNET_EI    = if (is.null(r$ei)) NA_real_ else as.numeric(r$ei),
  DisGeNET_nPMID = if (is.null(r$numPMIDs)) NA_integer_ else as.integer(r$numPMIDs)
)), fill = TRUE)
dg_dt <- unique(dg_dt, by = "SYMBOL")
message("ASoC genes with a DisGeNET C2239176 association: ", nrow(dg_dt))

## ---- Merge + classify ----
x <- merge(asoc, ot_dt, by = "ensembl", all.x = TRUE)
x <- merge(x, dg_dt, by = "SYMBOL", all.x = TRUE)
for (col in c("OT_HCC_overall", "CGC_HCC_score", "OT_clinvar", "OT_intogen", "OT_europepmc"))
  x[is.na(get(col)), (col) := 0]
x[, COSMIC_CGC_HCC := fifelse(CGC_HCC_score > 0, "yes", "no")]
x[, DisGeNET_HCC   := fifelse(!is.na(DisGeNET_score), "yes", "no")]
x[, n_curated_sources := (COSMIC_CGC_HCC == "yes") + (DisGeNET_HCC == "yes")]
x[, curated_HCC_hit := fifelse(n_curated_sources > 0, "yes", "no")]

setcolorder(x, c("SYMBOL", "ensembl", "TCGA_LIHC_dir", "TCGA_LIHC_log2FC",
                 "COSMIC_CGC_HCC", "CGC_HCC_score", "DisGeNET_HCC",
                 "DisGeNET_score", "DisGeNET_EI", "DisGeNET_nPMID",
                 "n_curated_sources", "curated_HCC_hit",
                 "OT_HCC_overall", "OT_clinvar", "OT_intogen", "OT_europepmc"))
setorder(x, -n_curated_sources, -DisGeNET_score, -OT_HCC_overall, -CGC_HCC_score, na.last = TRUE)

out_tsv <- file.path(outdir, "ASoC_genes_HCC_curated_crossref_disgenet.tsv")
fwrite(x, out_tsv, sep = "\t")

## ---- Summary ----
message("\n===== Curated HCC gene-set cross-reference (native DisGeNET) =====")
message("ASoC genes total (Ensembl):     ", nrow(x))
message("In COSMIC CGC (HCC-mapped):      ", sum(x$COSMIC_CGC_HCC == "yes"))
message("In DisGeNET C2239176:            ", sum(x$DisGeNET_HCC == "yes"))
message("In BOTH curated sets:            ", sum(x$n_curated_sources == 2))
message("In EITHER curated set:           ", sum(x$curated_HCC_hit == "yes"))
message("\nGenes in >=1 curated HCC set (by DisGeNET score, then OT score):")
print(x[curated_HCC_hit == "yes",
        .(SYMBOL, TCGA_LIHC_dir, COSMIC_CGC_HCC, DisGeNET_HCC,
          DisGeNET_score, DisGeNET_EI, OT_HCC_overall)])
message("\nWritten: ", out_tsv)
