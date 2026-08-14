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
  # library(Signac)
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
if (
  interactive() &&
    (Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode"))
) {
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

# qs_save(
#   merged_integrated_liver_obj,
#   "merged_all_samples_integrated_seurat_obj.qs2",
#   nthreads = 4
# )
# raw_gex_seurat_list <-
#   qs_read(
#     "raw_all_samples_seurat_list.qs2",
#     nthreads = 8
#   )

# merged_liver_obj <-
#   merge(
#     x = raw_gex_seurat_list[[1]],
#     y = raw_gex_seurat_list[2:length(raw_gex_seurat_list)],
#     add.cell.ids = names(raw_gex_seurat_list)
#   )
# merged_liver_obj$orig.ident <- as.factor(merged_liver_obj$orig.ident)
# merged_liver_obj$preparation <- as.factor(merged_liver_obj$preparation)
merged_liver_obj <-
  qs_read("./annotation_script_package/marker_based/step1.qs2", nthreads = 8)
merged_liver_obj <-
  qs_read("./annotation_script_package/marker_based/annotated_seurat.qs2", nthreads = 8)
# merged_liver_obj <-
#   qs_read("merged_all_samples_integrated_seurat_obj.qs2", nthreads = 8)
# colnames(merged_liver_obj)
unique(merged_liver_obj$orig.ident)
VlnPlot(
  merged_liver_obj,
  features = c("percent.mt"),
  ncol = 1,
  group.by = "orig.ident",
  pt.size = 0  
)

# qs_save(
#   merged_liver_obj,
#   "merged_all_samples_integrated_seurat_obj.qs2",
#   nthreads = 8
# )
merged_liver_obj$category <- "Primary"
merged_liver_obj$category[str_detect(
  merged_liver_obj$orig.ident,
  "pre"
)] <- "Resistant"
merged_liver_obj$category[
  merged_liver_obj$orig.ident %in%
    c(
      "BC010",
      "BC004",
      "BC001",
      "BC009",
      "BC014",
      "BC015",
      "BC002",
      "BC003",
      "BC005",
      "BC007",
      "BC011",
      "BC006",
      "BC016",
      "BC008"
    )
] <- "Resistant"


Reductions(merged_liver_obj)
#  [1] "pca"                 "integrated.harmony"  "tsne.unintegrated"   "umap.unintegrated"   "integrated.cca"      "integrated.rpca"     "integrated.jointpca" "tsne.cca"            "tsne.rpca"
# [10] "tsne.jointpca"       "tsne.harmony"        "umap.cca"            "umap.rpca"           "umap.jointpca"       "umap.harmony"
colnames(merged_liver_obj@meta.data)
scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "celltype_broad",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = T,
  colors_use = "glasbey",
  shuffle = TRUE
)
scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "celltype_fine",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = TRUE,
  colors_use = "glasbey"
)
scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "preparation",
  label = F,
  label.size = 4,
  alpha = 1,
  pt.size = 0.2,
  colors_use = "alphabet2",
  shuffle = T
)
#####
merged_hepatocyte_obj <-
  subset(merged_liver_obj, subset = celltype_broad == "Hepatocyte")
df_agg <-
  AggregateExpression(
    merged_hepatocyte_obj,
    assays = "RNA",
    slot = "counts",
    group.by = c("orig.ident"),
    verbose = T,
    return_seurat = F
  )
df_agg <- df_agg$RNA
rownames(df_agg)
colnames(df_agg)
df_meta <-
  merged_liver_obj@meta.data
df_meta <-
  df_meta[!duplicated(df_meta$orig.ident), ]
df_meta_ordered <-
  df_meta[
    match(
      colnames(df_agg),
      gsub("_", "-", as.character(df_meta$orig.ident))
    ),
  ]
df_meta_ordered$category <- "Primary"
df_meta_ordered$category[str_detect(
  df_meta_ordered$orig.ident,
  "pre"
)] <- "Resistant"
df_meta_ordered$category[
  df_meta_ordered$orig.ident %in%
    c(
      "BC010",
      "BC004",
      "BC001",
      "BC009",
      "BC014",
      "BC015",
      "BC002",
      "BC003",
      "BC005",
      "BC007",
      "BC011",
      "BC006",
      "BC016",
      "BC008"
    )
] <- "Resistant"


library(edgeR)
DGE_agg <-
  DGEList(
    counts = df_agg,
    samples = df_meta_ordered$orig.ident,
    group = df_meta_ordered$category
  )
DGE_agg <-
  calcNormFactors(
    DGE_agg,
    method = "TMM"
  )
DGE_cpm <-
  cpm(
    DGE_agg,
    log = TRUE,
    prior.count = 1
  )

genes_2_plot <-
  read.table(
    "sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated_HCC_crossref_DisGeNET_yes.tsv",
    header = T,
    sep = "\t"
  )
genes_2_plot <- genes_2_plot$SYMBOL
genes_2_plot <- genes_2_plot[!duplicated(genes_2_plot)]

df_2_plot <-
  DGE_cpm[rownames(DGE_cpm) %in% genes_2_plot, ]
df_2_plot <-
  as.data.frame(t(df_2_plot))
rownames(df_2_plot)
colnames(df_2_plot)

df_long <-
  df_2_plot |>
  tibble::rownames_to_column("sample_name") |>
  data.table::as.data.table() |>
  data.table::melt(
    id.vars = "sample_name",
    variable.name = "gene_name",
    value.name = "exp_level"
  )

df_long$category <- "Primary"
df_long$category[str_detect(
  df_long$sample_name,
  "pre"
)] <- "Resistant"
df_long$category[
  df_long$sample_name %in%
    c(
      "BC010",
      "BC004",
      "BC001",
      "BC009",
      "BC014",
      "BC015",
      "BC002",
      "BC003",
      "BC005",
      "BC007",
      "BC011",
      "BC006",
      "BC016",
      "BC008"
    )
] <- "Resistant"

df_long$preparation <- "10x"
df_long$preparation[str_detect(
  df_long$sample_name,
  "^BC"
)] <- "FLEX"
df_long$preparation <- as.factor(df_long$preparation)
genes_2_plot <- unique(df_long$gene_name)


DE_plot_list <-
  lapply(
    genes_2_plot,
    function(gene) {
      df_sub <-
        df_long[df_long$gene_name == gene, ]
      p_box <-
        ggplot(
          df_sub,
          aes(
            x = category,
            y = exp_level,
            fill = category
            # colour = as.factor(preparation)
          )
        ) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.7,
          outlier.size = 0,
          alpha = 0.5
        ) +
        geom_point(
          aes(colour = preparation),
          position = position_jitter(
            width = 0.05
          ),
          size = 1,
          alpha = 1
        ) +
        # dashed horizontal line at the mean of each box
        stat_summary(
          fun = mean,
          geom = "errorbar",
          aes(ymax = after_stat(y), ymin = after_stat(y)),
          position = position_dodge(width = 0.8),
          width = 1,
          # linetype = "dashed",
          colour = "darkred"
        ) +
        scale_fill_manual(values = c("orange4", "darkblue")) +
        scale_colour_manual(values = c("black", "darkred")) +
        ggpubr::stat_compare_means(
          # method = "wilcox.test",
          method = "t.test",
          comparisons = list(c("Primary", "Resistant")),
          label = "p.signif" # use "p.format" to show the actual p-value instead of significance level
        ) +
        expand_limits(y = 0) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        ylab("log2(CPM + 1)") +
        # xlab("") +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 315, hjust = 0),
          axis.title.x = element_blank()
        ) +
        labs(title = gene)
      return(p_box)
    }
  )

print(
  patchwork::wrap_plots(DE_plot_list, ncol = 4, guides = "collect") &
    theme(legend.position = "right")
)

table(multiome_obj$orig.ident)
# ggsave()

# merged_integrated_liver_obj <-
#   qs_read(
#     "merged_all_samples_integrated_seurat_obj.qs2",
#     nthreads = 8
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
