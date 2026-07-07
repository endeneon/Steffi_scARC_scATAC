#! /usr/bin/env Rscript

# Siwei 24 Mar 2026
# init ####
{
  library(BiocParallel)
  library(glmGamPoi)
  library(future)
  library(parallel)

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
  # library(loomR)
  # library(anndata)
  library(hdf5r)
  library(arrow)
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
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/murphygrp/SPATIAL/HB_Hypoxia_spatial"
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

# OpenMP-backed code (e.g. Rtsne) reads OMP_NUM_THREADS rather than future's
# worker count, so set it explicitly to use the full LSF allocation.
Sys.setenv(OMP_NUM_THREADS = available_cores - 1)

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
        20
      )
    )
}

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  set.seed(42)
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
    GEX_seurat <-
      GEX_seurat %>%
      SCTransform(verbose = T, method = "glmGamPoi") %>%
      NormalizeData(verbose = T, normalization.method = "LogNormalize") %>%
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
    ## run Harmony separately ####
    DefaultAssay(GEX_seurat) <- "RNA"
    GEX_seurat <-
      RunHarmony(
        GEX_seurat,
        group.by.vars = "orig.ident",
        reduction.use = "pca",
        dims.use = 1:n_dims,
        reduction.save = "integrated.harmony",
        verbose = T
      )

    DefaultAssay(GEX_seurat) <- normalization_method
    GEX_seurat <-
      GEX_seurat %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.cca",
        reduction.name = "tsne.cca",
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.rpca",
        reduction.name = "tsne.rpca",
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.jointpca",
        reduction.name = "tsne.jointpca",
        check_duplicates = F,
        verbose = T,
        seed.use = 42
      ) %>%
      RunTSNE(
        dims = 1:n_dims,
        reduction = "integrated.harmony",
        reduction.name = "tsne.harmony",
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
transformed_visium_list <-
  qs_read("transformed_visium_list.qs2", nthreads = workers_2_use / 2)

# subset the query object
query_liver_obj <-
  transformed_visium_list[c(
    "X1_1",
    'X5_1',
    "X4_1",
    "X11_1"
  )]

liver_seurat_list <-
  lapply(
    query_liver_obj,
    FUN = extract_raw_count_mtx,
    assay = "Spatial"
  )
names(liver_seurat_list) <-
  names(query_liver_obj)

for (i in seq_along(liver_seurat_list)) {
  liver_seurat_list[[i]] <-
    CreateSeuratObject(
      counts = liver_seurat_list[[i]],
      assay = "RNA",
      project = names(liver_seurat_list)[i]
    )
}

merged_liver_obj <-
  merge(
    x = liver_seurat_list[[1]],
    y = liver_seurat_list[2:length(liver_seurat_list)],
    add.cell.ids = names(liver_seurat_list)
  )

n_dims <-
  findPcsElbow(merged_liver_obj)

merged_integrated_liver_obj <-
  IterateIntegrateLayers(
    object = merged_liver_obj,
    normalization_method = "SCT",
    n_dims = n_dims
  )

# qs_save(
#   merged_integrated_liver_obj,
#   "merged_integrated_liver_obj.qs2",
#   nthreads = 4
# )

# merged_integrated_liver_obj <-
#   qs_read(
#     "merged_integrated_liver_obj.qs2",
#     nthreads = 4
#   )

merged_integrated_liver_obj$hypoxia.ident <- NA
merged_integrated_liver_obj$hypoxia.ident[
  merged_integrated_liver_obj$orig.ident %in% c("X1_1", "X5_1")
] <-
  "RA"
merged_integrated_liver_obj$hypoxia.ident[
  merged_integrated_liver_obj$orig.ident %in% c("X4_1", "X11_1")
] <-
  "IH"
merged_integrated_liver_obj$hypoxia.ident <-
  factor(
    merged_integrated_liver_obj$hypoxia.ident,
    levels = c("RA", "IH")
  )
merged_integrated_liver_obj$hypoxia.ident <-
  relevel(
    merged_integrated_liver_obj$hypoxia.ident,
    ref = "RA"
  )
merged_integrated_liver_obj$sex <- "NA"
merged_integrated_liver_obj$sex[
  merged_integrated_liver_obj$orig.ident %in% c("X1_1", "X4_1")
] <- "Male"
merged_integrated_liver_obj$sex[
  merged_integrated_liver_obj$orig.ident %in% c("X5_1", "X11_1")
] <- "Female"

# scCustomize::DimPlot_scCustom(
#   merged_integrated_liver_obj,
#   reduction = "umap.cca",
#   group.by = "orig.ident",
#   label = F,
#   label.size = 4,
#   pt.size = 0.2
# )

merged_integrated_liver_obj <-
  merged_integrated_liver_obj %>%
  FindNeighbors(
    reduction = "integrated.cca",
    dims = 1:n_dims,
    verbose = T
  ) %>%
  FindClusters(
    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5),
    verbose = T
  )
merged_integrated_liver_obj <-
  FindClusters(
    merged_integrated_liver_obj,
    resolution = c(
      0.25,
      0.26,
      0.27,
      0.28,
      0.29
    ),
    verbose = T
  )

merged_integrated_liver_obj <-
  IntegrateKmeansClustering(
    object = merged_integrated_liver_obj,
    clusters = 6,
    integrationMethod = "cca",
    seed.use = 42
  )

# sapply(merged_integrated_liver_obj@meta.data, dplyr::n_distinct)

# scCustomize::DimPlot_scCustom(
#   merged_integrated_liver_obj,
#   reduction = "umap.cca",
#   group.by = "SCT_snn_res.0.26",
#   split.by = "orig.ident",
#   colors_use = DiscretePalette_scCustomize(
#     num_colors = length(unique(merged_integrated_liver_obj$SCT_snn_res.0.26)),
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
merged_integrated_liver_obj[["RNA"]] <-
  JoinLayers(merged_integrated_liver_obj[["RNA"]])
merged_integrated_liver_obj <-
  merged_integrated_liver_obj %>%
  NormalizeData(
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = T
  )
cluster_levels <-
  levels(factor(merged_integrated_liver_obj$SCT_snn_res.0.26))

markers_list <-
  lapply(
    cluster_levels,
    function(cl) {
      obj_sub <-
        subset(
          merged_integrated_liver_obj,
          subset = SCT_snn_res.0.26 == cl
        )
      # Skip clusters that lack both conditions or have too few cells
      cond_counts <- table(obj_sub$hypoxia.ident)
      if (
        !all(c("IH", "RA") %in% names(cond_counts)) ||
          any(cond_counts[c("IH", "RA")] < 3)
      ) {
        return(NULL)
      }
      res <-
        tryCatch(
          FindMarkers(
            obj_sub,
            assay = "RNA",
            slot = "data",
            test.use = "MAST",
            ident.1 = "IH",
            ident.2 = "RA", # reference level, so that the logFC is positive for genes upregulated in IH
            group.by = "hypoxia.ident",
            min.cells.group = 3,
            latent.vars = c("sex"),
            logfc.threshold = 0.05,
            min.pct = 0.05,
            only.pos = F,
            verbose = T
          ),
          error = function(e) {
            warning(
              "FindMarkers failed for cluster ",
              cl,
              ": ",
              conditionMessage(e)
            )
            NULL
          }
        )
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }
      res$cluster <- cl
      tibble::rownames_to_column(res, "gene")
    }
  ) |>
  # Name each list item by its cluster and drop skipped (NULL) clusters
  setNames(cluster_levels) |>
  purrr::compact()

# Add a per-cluster FDR (BH) computed within each cluster's own p-values,
# then sort so the most significant genes appear first.
markers_list <-
  lapply(
    markers_list,
    function(df) {
      df$p_val_adj_BH <- p.adjust(df$p_val, method = "BH")
      df[order(df$p_val), ]
    }
  )

# Write each cluster's markers to its own tab-separated file in a dedicated
# output directory; the "cluster_" prefix keeps file names from starting with
# a numeric cluster label.
out_dir <- "DE_results_by_cluster"
if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
}

purrr::iwalk(
  markers_list,
  function(df, cl) {
    write.table(
      df,
      file = file.path(out_dir, paste0("cluster_", cl, "_markers.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
)

qs_save(
  merged_integrated_liver_obj,
  "merged_integrated_liver_obj.qs2",
  nthreads = 4
)
qs_save(
  markers_list,
  "IH_vs_RA_markers_list.qs2",
  nthreads = 4
)

# markers_list <-
#   qs_read(
#     "IH_vs_RA_markers_list.qs2",
#     nthreads = 4
#   )
