#! /usr/bin/env Rscript

# Siwei 24 Mar 2026
# init ####
{
  library(BiocParallel)
  library(glmGamPoi)
  library(future)
  library(parallel)
  library(foreach)
  library(doParallel)

  library(snow)
  library(harmony)
  # library(optparse)

  library(dplyr)
  library(data.table)

  library(edgeR)

  # library(scuttle)
  library(Matrix)
  library(matrixStats)

  library(Seurat)
  library(Signac)
  # library(loomR)
  # library(anndata)
  # library(hdf5r)
  # library(arrow)
  # library(rhdf5)

  library(stringr)

  library(qs2)
  library(fs)

  # library(SingleR)

  library(ggplot2)
  library(gplots)
  library(patchwork)
  library(scCustomize)

  # library(scuttle)

  # library(scMerge)
  # library(scater)
}

# setwd(
#   "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/
#    szhang37/projects/szhang_dev/STEREO_seq/Human_liver/R_liver"
# )
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
print(available_cores)

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

if (session_plan == "multisession") {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        16
      )
    )
} else {
  workers_2_use <-
    max(
      1,
      min(
        available_cores - 1,
        32
      )
    )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
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

# functions ####
## extract raw count matrix from Seurat object
extract_raw_count_mtx <-
  function(
    object,
    assay = "RNA"
  ) {
    object <-
      UpdateSeuratObject(object)
    DefaultAssay(object) <- assay
    dimnames_count_mtx <-
      list(
        gene_id = rownames(object),
        barcodes = colnames(object)
      )
    if (
      is.null(dimnames_count_mtx$gene_id) |
        is.null(dimnames_count_mtx$barcodes)
    ) {
      return(NA)
    } else {
      raw_count_mtx <-
        object[[assay]]@layers$counts
      rownames(raw_count_mtx) <-
        dimnames_count_mtx$gene_id
      colnames(raw_count_mtx) <-
        dimnames_count_mtx$barcodes
      return(raw_count_mtx)
    }
  }

# find PCs
# use Surbhi's function with changes
# this function will return npcs only, no figures
findPcsElbow <-
  function(obj) {
    # credit: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
    # calculate where the principal components start to elbow by taking the larger value of:
    #  1. The point where the principal components only contribute 5% of standard deviation and
    # 		the principal components cumulatively contribute 90% of the standard deviation.
    #  2. The point where the percent change in variation between the consecutive PCs is less than 0.1%
    # 		Determine percent of variation associated with each PC
    stopifnot(is(obj, "Seurat"))

    if (!"pca" %in% Reductions(obj)) {
      obj <-
        obj %>%
        NormalizeData(verbose = T) %>%
        FindVariableFeatures(verbose = T) %>%
        ScaleData(verbose = T) %>%
        RunPCA(verbose = T)
    }

    pct <-
      obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <-
      which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <-
      sort(
        which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
        decreasing = T
      )[1] +
      1
    # Minimum of the two calculation; i.e. PCs covering the majority of the variation in the data.
    pcs <-
      min(co1, co2)
    return(ndims = pcs)
  }

############
IterateIntegrateLayers <-
  function(
    object,
    normalization_method = "SCT",
    n_dims = 20
  ) {
    GEX_seurat <- object
    # need to run both SCT and log-normalization for all integration methods to work
    # NOTE: do NOT JoinLayers here. Seurat v5 IntegrateLayers() needs the RNA
    # assay to keep its per-sample layers split so it can find pairwise anchors
    # (and so SCTransform builds one SCT model per sample). Joining the layers
    # collapses the object to a single dataset, which makes FindIntegrationAnchors
    # return NULL and crashes with:
    #   "no applicable method for 'Assays' applied to an object of class NULL".
    # Join layers AFTER all IntegrateLayers() calls (e.g. before marker DE).
    ## run Harmony separately ####
    DefaultAssay(GEX_seurat) <- "RNA"

    GEX_seurat <-
      GEX_seurat %>%
      NormalizeData(
        assay = "RNA",
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = T
      ) %>%
      FindVariableFeatures(verbose = T, selection.method = "vst") %>%
      ScaleData(verbose = T) %>%
      RunPCA(verbose = T, npcs = n_dims)

    GEX_seurat <-
      RunHarmony(
        GEX_seurat,
        group.by.vars = c("orig.ident", "preparation"),
        reduction.use = "pca",
        dims.use = 1:n_dims,
        reduction.save = "integrated.harmony",
        verbose = T
      )

    ## run SCT separately ####
    GEX_seurat <-
      GEX_seurat %>%
      SCTransform(
        verbose = T,
        method = "glmGamPoi",
        vars.to.regress = c("orig.ident", "preparation")
      ) %>%
      FindVariableFeatures(verbose = T, selection.method = "vst") %>%
      ScaleData(verbose = T) %>%
      RunPCA(verbose = T, npcs = n_dims)

    GEX_seurat <-
      GEX_seurat %>%
      FindNeighbors(reduction = "pca", dims = 1:n_dims) %>%
      FindClusters(resolution = 1, random.seed = 42) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "pca",
        reduction.name = "tsne.unintegrated",
        num_threads = omp_threads,
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunUMAP(
        dims = 1:n_dims,
        reduction = "pca",
        reduction.name = "umap.unintegrated",
        verbose = T,
        seed.use = 42
      ) %>%
      IntegrateLayers(
        method = CCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.cca",
        normalization.method = normalization_method,
        verbose = T,
        dims = 1:n_dims,
        preserve.order = T
      ) %>%
      IntegrateLayers(
        method = RPCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.rpca",
        normalization.method = normalization_method,
        verbose = T,
        dims = 1:n_dims,
        preserve.order = T
      ) %>%
      IntegrateLayers(
        method = JointPCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.jointpca",
        normalization.method = normalization_method,
        verbose = T,
        dims = 1:n_dims,
        preserve.order = T
      )

    DefaultAssay(GEX_seurat) <- normalization_method
    GEX_seurat <-
      GEX_seurat %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.cca",
        reduction.name = "tsne.cca",
        num_threads = omp_threads,
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.rpca",
        reduction.name = "tsne.rpca",
        num_threads = omp_threads,
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.jointpca",
        reduction.name = "tsne.jointpca",
        num_threads = omp_threads,
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.harmony",
        reduction.name = "tsne.harmony",
        num_threads = omp_threads,
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunUMAP(
        dims = 1:n_dims,
        reduction = "integrated.cca",
        reduction.name = "umap.cca",
        verbose = T,
        seed.use = 42
      ) %>%
      RunUMAP(
        dims = 1:n_dims,
        reduction = "integrated.rpca",
        reduction.name = "umap.rpca",
        verbose = T,
        seed.use = 42
      ) %>%
      RunUMAP(
        dims = 1:n_dims,
        reduction = "integrated.jointpca",
        reduction.name = "umap.jointpca",
        verbose = T,
        seed.use = 42
      ) %>%
      RunUMAP(
        dims = 1:n_dims,
        reduction = "integrated.harmony",
        reduction.name = "umap.harmony",
        verbose = T,
        seed.use = 42
      )
    return(GEX_seurat)
  }

IntegrateKmeansClustering <-
  function(
    object,
    clusters = 12,
    integrationMethod = "cca",
    seed.use = 42
  ) {
    GEX_seurat <- object
    set.seed(seed = seed.use)
    integration_dim <-
      paste0(
        'integrated.',
        integrationMethod
      )
    print(integration_dim)
    matrix_4_kmeans <-
      GEX_seurat@reductions[[integration_dim]]@cell.embeddings
    kmeans_result <-
      stats::kmeans(
        x = matrix_4_kmeans,
        centers = clusters,
        nstart = 50
      )
    GEX_seurat$kmeans_clusters <-
      kmeans_result$cluster
    return(GEX_seurat)
  }

# load data ####

# liver_seurat_list <-
#   lapply(
#     query_liver_obj,
#     FUN = extract_raw_count_mtx,
#     assay = "RNA"
#   )
# names(liver_seurat_list) <-
#   names(query_liver_obj)

# for (i in seq_along(liver_seurat_list)) {
#   liver_seurat_list[[i]] <-
#     CreateSeuratObject(
#       counts = liver_seurat_list[[i]],
#       assay = "RNA",
#       project = names(liver_seurat_list)[i]
#     )
# }

# load multiome GEX (Gene Expression) data ####
## 1.1) locate all cellranger-arc filtered_feature_bc_matrix.h5 files
multiome_bams_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/test_ASoC_w_WASP/multiome_bams"
multiome_h5_file_vec <-
  list.files(
    path = multiome_bams_dir,
    pattern = "^filtered_feature_bc_matrix\\.h5$",
    full.names = TRUE,
    recursive = TRUE
  )
names(multiome_h5_file_vec) <-
  paste0(
    "X__",
    basename(dirname(dirname(multiome_h5_file_vec)))
  )
## 1.2) locate all FLEX filtered_feature_bc_matrix.h5 files
flex_bams_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/FLEX_samples"
flex_h5_file_vec <-
  list.files(
    path = flex_bams_dir,
    pattern = "^sample_filtered_feature_bc_matrix\\.h5$",
    full.names = TRUE,
    recursive = TRUE
  )
names(flex_h5_file_vec) <-
  basename(dirname(dirname(flex_h5_file_vec)))

## 1.3) locate all GEX filtered_feature_bc_matrix.h5 files
gex_bams_dir <-
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/GEX_samples"
gex_h5_file_vec <-
  list.files(
    path = gex_bams_dir,
    pattern = "^filtered_feature_bc_matrix\\.h5$",
    full.names = TRUE,
    recursive = TRUE
  )
names(gex_h5_file_vec) <-
  basename(dirname(dirname(gex_h5_file_vec)))

## 2) merge all h5_file_vec to a master vector
h5_file_vec <- c(multiome_h5_file_vec, flex_h5_file_vec, gex_h5_file_vec)
names(h5_file_vec) <- c(
  names(multiome_h5_file_vec),
  names(flex_h5_file_vec),
  names(gex_h5_file_vec)
)
#   paste0(
#     "X__",
#     basename(dirname(dirname(h5_file_vec)))
#   )
# test_data <-
#   Read10X_h5("GEX_samples/GEX_pt1pre/outs/filtered_feature_bc_matrix.h5")
## 3)-4) load GEX-only data, build the Seurat object, and add QC percentages
load_gex_seurat <-
  function(h5_path, sample_name) {
    data_10x <-
      Read10X_h5(filename = h5_path)
    ## extract the GEX counts only from the "Gene Expression" slot
    rna_counts <-
      if ("Gene Expression" %in% names(data_10x)) {
        data_10x[["Gene Expression"]]
      } else {
        data_10x
      }
    seurat_obj <-
      CreateSeuratObject(
        counts = rna_counts,
        assay = "RNA",
        project = sample_name
      )
    ## set orig.ident to the (prefixed) sample name
    seurat_obj$orig.ident <- sample_name
    seurat_obj$preparation <-
      dplyr::case_when(
        startsWith(sample_name, "X__") ~ "multiome",
        startsWith(sample_name, "BC") ~ "FLEX",
        TRUE ~ "GEX"
      )

    ## human mitochondrial and ribosomal gene patterns
    seurat_obj[["percent.mt"]] <-
      PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj[["percent.ribo"]] <-
      PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
    return(seurat_obj)
  }

## 5) process all samples in parallel and return a named list of Seurat objects

n_cluster_workers <-
  max(
    1,
    min(
      length(h5_file_vec),
      workers_2_use - 2
    )
  )
gex_cluster <-
  parallel::makeCluster(n_cluster_workers)
doParallel::registerDoParallel(gex_cluster)

raw_gex_seurat_list <-
  foreach(
    i = seq_along(h5_file_vec),
    .packages = c("Seurat")
  ) %dopar%
  {
    load_gex_seurat(
      h5_path = h5_file_vec[i],
      sample_name = names(h5_file_vec)[i]
    )
  }

parallel::stopCluster(gex_cluster)

names(raw_gex_seurat_list) <-
  names(h5_file_vec)

qs_save(
  raw_gex_seurat_list,
  "raw_all_samples_seurat_list.qs2",
  nthreads = 4
)

# merged_liver_obj <-
#   merge(
#     x = raw_gex_seurat_list[[1]],
#     y = raw_gex_seurat_list[2:length(raw_gex_seurat_list)],
#     add.cell.ids = names(raw_gex_seurat_list)
#   )

# merged_liver_obj[["RNA"]] <-
#   JoinLayers(merged_liver_obj[["RNA"]])
# merged_liver_obj <-
#   merged_liver_obj %>%
#   NormalizeData(
#     assay = "RNA",
#     normalization.method = "LogNormalize",
#     scale.factor = 10000,
#     verbose = T
#   )

n_dims <-
  findPcsElbow(merged_liver_obj)

merged_integrated_liver_obj <-
  IterateIntegrateLayers(
    object = merged_liver_obj,
    normalization_method = "SCT",
    n_dims = n_dims
  )

qs_save(
  merged_integrated_liver_obj,
  "merged_all_samples_integrated_seurat_obj.qs2",
  nthreads = 4
)
# merged_integrated_liver_obj <-
#   qs_read(
#     "merged_integrated_liver_obj.qs2",
#     nthreads = 4
#   )

# scCustomize::DimPlot_scCustom(
#   merged_integrated_liver_obj,
#   reduction = "umap.cca",
#   group.by = "orig.ident",
#   label = F,
#   label.size = 4,
#   pt.size = 0.2
# )

# merged_integrated_liver_obj <-
#   merged_integrated_liver_obj %>%
#   FindNeighbors(
#     reduction = "integrated.cca",
#     dims = 1:n_dims,
#     verbose = T
#   ) %>%
#   FindClusters(
#     resolution = c(0.1, 0.2, 0.5),
#     verbose = T
#   )

# merged_integrated_liver_obj <-
#   IntegrateKmeansClustering(
#     object = merged_integrated_liver_obj,
#     clusters = 6,
#     integrationMethod = "cca",
#     seed.use = 42
#   )

# sapply(merged_integrated_liver_obj@meta.data, dplyr::n_distinct)

# scCustomize::DimPlot_scCustom(
#   merged_integrated_liver_obj,
#   reduction = "umap.cca",
#   group.by = "SCT_snn_res.0.1",
#   split.by = "orig.ident",
#   colors_use = DiscretePalette_scCustomize(
#     num_colors = length(unique(merged_integrated_liver_obj$SCT_snn_res.0.1)),
#     palette = "glasbey"
#   ),
#   label = F,
#   label.size = 4,
#   pt.size = 0.2
# )

# scCustomize::DimPlot_scCustom(
#   merged_integrated_liver_obj,
#   reduction = "umap.cca",
#   group.by = "kmeans_clusters",
#   split.by = "orig.ident",
#   colors_use = DiscretePalette_scCustomize(
#     n = 6,
#     palette = "glasbey"
#   ),
#   label = F,
#   label.size = 4,
#   pt.size = 0.2
# )

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# merged_integrated_liver_obj[["RNA"]] <-
#   JoinLayers(merged_integrated_liver_obj[["RNA"]])
# merged_integrated_liver_obj <-
#   merged_integrated_liver_obj %>%
#   NormalizeData(
#     assay = "RNA",
#     normalization.method = "LogNormalize",
#     scale.factor = 10000,
#     verbose = T
#   )
# cluster_levels <-
#   levels(factor(merged_integrated_liver_obj$SCT_snn_res.0.1))

# markers_list <-
#   lapply(
#     cluster_levels,
#     function(cl) {
#       obj_sub <-
#         subset(
#           merged_integrated_liver_obj,
#           subset = SCT_snn_res.0.1 == cl
#         )
#       # Skip clusters that lack both conditions or have too few cells
#       cond_counts <- table(obj_sub$hypoxia.ident)
#       if (
#         !all(c("IH", "RA") %in% names(cond_counts)) ||
#           any(cond_counts[c("IH", "RA")] < 3)
#       ) {
#         return(NULL)
#       }
#       res <-
#         tryCatch(
#           FindMarkers(
#             obj_sub,
#             assay = "RNA",
#             slot = "data",
#             test.use = "MAST",
#             ident.1 = "IH",
#             ident.2 = "RA", # reference level, so that the logFC is positive for genes upregulated in IH
#             group.by = "hypoxia.ident",
#             min.cells.group = 3,
#             latent.vars = c("sex"),
#             logfc.threshold = 0.05,
#             min.pct = 0.05,
#             only.pos = F,
#             verbose = T
#           ),
#           error = function(e) {
#             warning(
#               "FindMarkers failed for cluster ",
#               cl,
#               ": ",
#               conditionMessage(e)
#             )
#             NULL
#           }
#         )
#       if (is.null(res) || nrow(res) == 0) {
#         return(NULL)
#       }
#       res$cluster <- cl
#       tibble::rownames_to_column(res, "gene")
#     }
#   ) |>
#   # Name each list item by its cluster and drop skipped (NULL) clusters
#   setNames(cluster_levels) |>
#   purrr::compact()

# # Add a per-cluster FDR (BH) computed within each cluster's own p-values,
# # then sort so the most significant genes appear first.
# markers_list <-
#   lapply(
#     markers_list,
#     function(df) {
#       df$p_val_adj_BH <- p.adjust(df$p_val, method = "BH")
#       df[order(df$p_val), ]
#     }
#   )

# # Write each cluster's markers to its own tab-separated file in a dedicated
# # output directory; the "cluster_" prefix keeps file names from starting with
# # a numeric cluster label.
# out_dir <- "DE_results_by_cluster"
# if (!dir.exists(out_dir)) {
#   dir.create(
#     out_dir,
#     showWarnings = FALSE,
#     recursive = TRUE
#   )
# }

# purrr::iwalk(
#   markers_list,
#   function(df, cl) {
#     write.table(
#       df,
#       file = file.path(out_dir, paste0("cluster_", cl, "_markers.tsv")),
#       sep = "\t",
#       quote = FALSE,
#       row.names = FALSE
#     )
#   }
# )

# qs_save(
#   merged_integrated_liver_obj,
#   "merged_integrated_seurat_obj.qs2",
#   nthreads = 4
# )
# qs_save(
#   markers_list,
#   "IH_vs_RA_markers_list.qs2",
#   nthreads = 4
# )

# markers_list <-
#   qs_read(
#     "IH_vs_RA_markers_list.qs2",
#     nthreads = 4
#   )
