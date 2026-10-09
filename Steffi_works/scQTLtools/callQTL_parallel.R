# Parallel drop-in for scQTLtools::callQTL() (scQTLtools 1.2.4).
#
# Same arguments and same output as callQTL(), plus `n_cores`. The original
# loops over SNPs with lapply(); here the per-SNP work is run with
# parallel::mclapply() (forked workers), reusing scQTLtools' own internal
# per-SNP functions so the statistics are unchanged.
#
# Notes
# - Forking only works on Linux/macOS and is unsafe inside RStudio/Positron
#   GUI sessions; run under Rscript (e.g. bsub) for the real job.
# - Workers share the parent's matrices via copy-on-write, so memory grows
#   with the per-worker results, not with n_cores x matrix size.
# - The internal functions are accessed with `:::`; if a scQTLtools update
#   renames them this will fail loudly rather than silently change results.

callQTL_parallel <-
  function(
    eQTLObject,
    gene_ids = NULL,
    downstream = NULL,
    upstream = NULL,
    gene_mart = NULL,
    snp_mart = NULL,
    pAdjustMethod = "bonferroni",
    useModel = "zinb",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.1,
    n_cores = 1L
  ) {
    ns <- asNamespace("scQTLtools")
    internal <- function(fn) get(fn, envir = ns, inherits = FALSE)

    filter_data <- scQTLtools::get_filter_data(eQTLObject)
    if (length(filter_data) == 0) {
      stop("Please filter the data first.")
    }
    expressionMatrix <- filter_data[["expMat"]]
    snpMatrix <- filter_data[["snpMat"]]
    biClassify <- scQTLtools::load_biclassify_info(eQTLObject)
    species <- scQTLtools::load_species_info(eQTLObject)

    if (is.null(species) || species == "") {
      stop("The 'species' variable is NULL or empty.")
    }
    datasets <-
      switch(
        species,
        human = c("hsapiens_snp", "hsapiens_gene_ensembl", "org.Hs.eg.db"),
        mouse = c("mmusculus_snp", "mmusculus_gene_ensembl", "org.Mm.eg.db"),
        stop("Please enter 'human' or 'mouse'.")
      )
    if (!useModel %in% c("zinb", "poisson", "linear")) {
      stop("Invalid model Please choose from 'zinb','poisson',or 'linear'.")
    }
    if (
      !pAdjustMethod %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")
    ) {
      stop("Invalid p-adjusted method.")
    }

    matchID <-
      internal("match_gene_snp")(
        gene_ids,
        upstream,
        downstream,
        rownames(snpMatrix),
        rownames(expressionMatrix),
        gene_mart,
        snp_mart,
        datasets[2],
        datasets[1],
        datasets[3]
      )
    geneIDs <- matchID[["matched_gene"]]
    snpIDs <- matchID[["matched_snps"]]

    # Pseudo-counts as in the original (zinb, poisson and linear all use this)
    expressionMatrix <- round(expressionMatrix * 1000)
    # Restrict to tested genes up front: smaller matrix shared with workers
    expressionMatrix <- expressionMatrix[geneIDs, , drop = FALSE]

    group_info <- scQTLtools::load_group_info(eQTLObject)
    unique_group <- unique(group_info[["group"]])

    message(
      "Start the sc-eQTLs calling: ",
      length(snpIDs),
      " SNPs x ",
      length(geneIDs),
      " genes x ",
      length(unique_group),
      " group(s), ",
      n_cores,
      " cores."
    )

    result_all <-
      do.call(
        rbind,
        lapply(unique_group, function(group) {
          split_cells <- rownames(group_info)[group_info[["group"]] == group]
          exp_split <- expressionMatrix[, split_cells, drop = FALSE]
          snp_split <- snpMatrix[, split_cells, drop = FALSE]

          if (useModel == "zinb") {
            if (biClassify) {
              snp_split[snp_split == 3] <- 2
            }
            per_snp <- function(i) {
              internal("eQTLcalling")(
                i,
                snpIDs,
                exp_split,
                snp_split,
                geneIDs,
                biClassify
              )
            }
            snp_index <- seq_along(snpIDs)
          } else {
            snp_split <- internal("process_snp_matrix")(snp_split, biClassify)
            per_snp <- function(snpID) {
              internal("compute_gene_results_common")(
                snpID,
                geneIDs,
                exp_split,
                snp_split,
                group,
                useModel
              )
            }
            snp_index <- snpIDs
          }

          message(group, ": started at ", format(Sys.time()))
          # mc.preschedule = FALSE: SNPs differ a lot in run time, so hand
          # them out one at a time for better load balance
          res_list <-
            parallel::mclapply(
              snp_index,
              per_snp,
              mc.cores = n_cores,
              mc.preschedule = FALSE
            )
          failed <- vapply(res_list, inherits, logical(1), what = "try-error")
          if (any(failed)) {
            stop(
              sum(failed),
              " SNP(s) failed in group ",
              group,
              ". First error: ",
              as.character(res_list[[which(failed)[1]]])
            )
          }
          message(group, ": finished at ", format(Sys.time()))

          res <- do.call(rbind, res_list)
          if (useModel == "zinb") {
            res <- internal("rename_columns")(res, biClassify)
            res$group <- group
          }
          res
        })
      )

    message("Finished!")
    result_all <-
      internal("adjust_pvalues")(result_all, pAdjustMethod, pAdjustThreshold)
    if (useModel != "zinb") {
      result_all <- internal("filter_by_abs_b")(result_all, logfcThreshold)
    }

    eQTLObject <- scQTLtools::set_model_info(eQTLObject, useModel)
    eQTLObject <- scQTLtools::set_result_info(eQTLObject, result_all)
    eQTLObject
  }
