#! /usr/bin/env Rscript

# init ####

{
  # library(devtools)
  # library(reticulate)
  #
  # library(SingCellaR)
  library(Seurat)
  library(scCustomize)
  # library(loupeR)
  # library(scDblFinder)
  # library(monocle)

  library(BiocParallel)
  library(glmGamPoi)
  # library(harmony)
  library(parallel)
  library(future)

  library(Matrix)

  library(stringr)
  library(magrittr)
  library(qs2)

  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(ggrastr)

  library(SingleR)
  library(scuttle)

  if (Sys.getenv("VSCODE") == "1") {
    library(languageserver)
    library(showtext)
    library(httpgd)

    showtext::showtext_auto()
  }
}

setwd(
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/murphygrp/SNUCRNASEQ/D6_Mice_Hypoxia/Priya_additional_modules"
)

# showtext::showtext_auto()
# setwd("/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/projects/szhang_dev/STEREO_seq/Human_liver/R_liver")
# determine if R is running in RSTUDIO
if (Sys.getenv("RSTUDIO") == "1" || Sys.getenv("POSITRON") == "1") {
  print("Running under RStudio IDE, use plan(multisession)")
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

## nthreads, need to take care that POSITRON cannot fully figure out the threads to use, so we need to set it manually here.
## We will use all but one thread for the analysis, and leave one thread for the system to use.
workers_2_use <-
  ifelse(
    test = session_plan == "multisession",
    yes = 16,
    no = parallel::detectCores() - 1
  )

{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  set.seed(42)
  options(future.globals.maxSize = workers_2_use * 20 * 1024^3) # 40 G per thread
  future::plan(
    session_plan, # Do NOT use "multisession" here if use LSF, use "multicore" instead
    workers = workers_2_use
  )
}

# load the sketched dataset, Male ####
seurat_obj_2_plot <-
  qs_read(
    file = "Hypoxia_2_obj_D6Male_IH_RAS_seurat_obj_integrated_label_singleR.qs2",
    nthreads = 8
  )

# plot module scores ####
module_list <-
  list(
    "ROS-Hnf4a" = c(
      "Hnf4a",
      "Mapk8",
      "Mapk9",
      "Elk1",
      "Slu7",
      "Tbx3",
      "E2f7",
      "E2f8",
      "Myc",
      "Ctnnb1"
    ),
    "Feed-loop" = c(
      "Lin28b",
      "Mirlet7a-1",
      "Mirlet7a-2",
      "Hmga2",
      "Igf2bp1",
      "Igf2bp2",
      "Igf2bp3",
      "Yap1",
      "Hif1a",
      "Nfkb1",
      "Nfkb2",
      "Myc",
      "Ctnnb1"
    ),
    "Redox-PI3K" = c(
      "Pik3ca",
      "Pik3c2a",
      "Pik3c2b",
      "Pik3r6",
      "Akt1",
      "Akt2",
      "Pten",
      "Foxo1",
      "Gsk3b",
      "G6pc1",
      "G6pc2",
      "G6pc3",
      "Pck1",
      "Pck2"
    )
  )

scCustomize::DimPlot_scCustom(seurat_obj_2_plot,
  reduction = "umap.cca",
  group.by = "RNA_snn_res.0.1",
  colors_use = DiscretePalette(
    n = length(unique(seurat_obj_2_plot$RNA_snn_res.0.1)),
    palette = "glasbey"
  ),
  add_prop_plot = F
)

for (module_name in names(module_list)) {
  print(str_c("Calculating module score for ", module_name, "..."))
  module_genes <- module_list[[module_name]]
  module_score_name <- str_c(module_name, "ModuleScore", sep = "_")
  seurat_obj_2_plot <-
    AddModuleScore(
      object = seurat_obj_2_plot,
      features = list(module_genes),
      name = module_score_name,
      ctrl = 100
    )
}
# Warning: The following features are not present in the object: Mirlet7a-1, Mirlet7a-2, not searching for symbol synonyms
# Warning: The following features are not present in the object: G6pc1, G6pc2, not searching for symbol synonyms

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[1]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[1]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[1]][1:5],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident"
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(text = element_text(family = "sans"))
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[1]][6:10],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident"
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(text = element_text(family = "sans"))
}

## Vln + DotPlot ####

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[11]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[11]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[1]],
        colnames(seurat_obj_2_plot@meta.data)[[11]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}

## module list 2 ####
{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[2]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[2]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[2]][c(1, 4:7)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(text = element_text(family = "sans"))
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[2]][8:13],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(text = element_text(family = "sans"))
}

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[12]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[12]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[2]],
        colnames(seurat_obj_2_plot@meta.data)[[12]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}

## module list 3 ####
{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[3]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[3]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.cca",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(
      text = element_text(family = "sans")
    )
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[3]][c(1:6)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

# "Pten"    "Foxo1"   "Gsk3b"   "G6pc1"   "G6pc2"   "G6pc3"   "Pck1"    "Pck2"
{
  p_vlnlist <-
    sapply(
      X = module_list[[3]][c(7:9, 12:14)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[13]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[13]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[3]],
        colnames(seurat_obj_2_plot@meta.data)[[13]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}


# load the sketched dataset, Female ####
seurat_obj_2_plot <-
  qs_read(
    file = "D6Female_IH_RAS_seurat.qs2",
    nthreads = 8
  )
seurat_obj_2_plot$orig.ident <-
  str_remove_all(
    seurat_obj_2_plot$orig.ident,
    "_Liver"
  )

for (module_name in names(module_list)) {
  print(str_c("Calculating module score for ", module_name, "..."))
  module_genes <- module_list[[module_name]]
  module_score_name <- str_c(module_name, "ModuleScore", sep = "_")
  seurat_obj_2_plot <-
    AddModuleScore(
      object = seurat_obj_2_plot,
      features = list(module_genes),
      name = module_score_name,
      ctrl = 100
    )
}
length(unique(seurat_obj_2_plot$RNA_snn_res.0.1))
colnames(seurat_obj_2_plot@meta.data)

## module list 1 ####
{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[1]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[1]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(
      text = element_text(family = "sans")
    )
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[1]][c(1:5)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

# "Pten"    "Foxo1"   "Gsk3b"   "G6pc1"   "G6pc2"   "G6pc3"   "Pck1"    "Pck2"
{
  p_vlnlist <-
    sapply(
      X = module_list[[1]][c(6:10)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[69]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[69]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[1]],
        colnames(seurat_obj_2_plot@meta.data)[[69]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}

## module list 2 ####
{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[2]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[2]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(
      text = element_text(family = "sans")
    )
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[2]][c(1, 4:7)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

# "Pten"    "Foxo1"   "Gsk3b"   "G6pc1"   "G6pc2"   "G6pc3"   "Pck1"    "Pck2"
{
  p_vlnlist <-
    sapply(
      X = module_list[[2]][c(8:13)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[70]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[70]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[2]],
        colnames(seurat_obj_2_plot@meta.data)[[70]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}


## module list 3 ####
{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[1]],
      features = module_list[[3]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(text = element_text(family = "sans"))
}

{
  p_featureplot_genes <-
    scCustomize::FeaturePlot_scCustom(
      seurat_obj_2_plot[, seurat_obj_2_plot$orig.ident == unique(seurat_obj_2_plot$orig.ident)[2]],
      features = module_list[[3]], # SCF
      colors_use = viridisLite::viridis(
        n = 200,
        direction = -1
      ),
      pt.size = 0.01,
      reduction = "umap.harmony",
      # na_cutoff = 500,
      alpha_exp = 0.5,
      split.by = "orig.ident",
      # num_columns = 4,
      combine = F
    )

  list_plots <-
    vector(
      mode = "list",
      length = length(p_featureplot_genes)
    )

  for (i in seq_along(list_plots)) {
    list_plots[[i]] <-
      p_featureplot_genes[[i]]
  }

  patchwork::wrap_plots(list_plots, ncol = 3) &
    theme(
      text = element_text(family = "sans")
    )
}

{
  p_vlnlist <-
    sapply(
      X = module_list[[3]][c(1:6)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

# "Pten"    "Foxo1"   "Gsk3b"   "G6pc1"   "G6pc2"   "G6pc3"   "Pck1"    "Pck2"
{
  p_vlnlist <-
    sapply(
      X = module_list[[3]][c(7:9, 12:14)],
      FUN = function(gene) {
        scCustomize::VlnPlot_scCustom(
          seurat_obj_2_plot,
          features = gene,
          group.by = "RNA_snn_res.0.1",
          split.by = "orig.ident",
          pt.size = 0.1
        )
      }
    )
  patchwork::wrap_plots(
    p_vlnlist,
    ncol = 1
  ) +
    patchwork::plot_layout(guides = "collect") &
    theme(
      text = element_text(family = "sans"),
      axis.title.x = element_blank()
    )
}

{
  p_module_vlnplot <-
    scCustomize::VlnPlot_scCustom(
      seurat_obj_2_plot,
      features = colnames(seurat_obj_2_plot@meta.data)[[71]],
      group.by = "RNA_snn_res.0.1",
      pt.size = 0,
      split.by = "orig.ident"
    ) &
      labs(
        title = str_remove(
          colnames(seurat_obj_2_plot@meta.data)[[71]],
          "_ModuleScore_?\\d+"
        )
      ) &
      theme(
        text = element_text(family = "sans"),
        axis.text.x = element_text(angle = 0, hjust = 0)
      )
  seurat_obj_2_plot_tmp <- seurat_obj_2_plot
  seurat_obj_2_plot_tmp$cluster_ident <-
    str_c(
      seurat_obj_2_plot_tmp$RNA_snn_res.0.1,
      "_",
      seurat_obj_2_plot_tmp$orig.ident
    )
  p_dotplot <-
    scCustomize::DotPlot_scCustom(
      seurat_obj_2_plot_tmp,
      features = c(
        module_list[[3]],
        colnames(seurat_obj_2_plot@meta.data)[[71]]
      ),
      colors_use = viridisLite::plasma(
        n = 200,
        direction = 1
      ),
      group.by = "cluster_ident",
      flip_axes = T,
      facet_label_rotate = T
    ) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 315, hjust = 0),
      text = element_text(family = "sans")
    ) +
    scale_x_discrete(
      labels = function(x) str_remove(x, "_ModuleScore_?\\d+"),
      limits = rev
    )

  patchwork::wrap_plots(
    p_module_vlnplot,
    p_dotplot,
    ncol = 2,
    widths = c(1, 2)
  ) &
    theme(text = element_text(family = "sans"))
}

## singleR annotation ####
## load mouse liver ref, Fabres et al. 2020
{
  mouse_liver_ref <-
    qs_read("../D6_Male_R_analysis/ref_mouse_liver_final_version.qs2",
      nthreads = 8
    )
  seurat_mouse_liver_ref <-
    CreateSeuratObject(
      counts = mouse_liver_ref@assays$RNA@layers$counts,
      assay = "RNA",
      meta.data = mouse_liver_ref@meta.data
    )
  rownames(seurat_mouse_liver_ref) <-
    rownames(mouse_liver_ref)
  colnames(seurat_mouse_liver_ref) <-
    colnames(mouse_liver_ref)

  sce_mouse_liver_ref <-
    as.SingleCellExperiment(seurat_mouse_liver_ref)
  sce_mouse_liver_ref <-
    logNormCounts(sce_mouse_liver_ref)
  unique(sce_mouse_liver_ref$cell_type__ontology_label)
}

singleR_return_annot_from_scratch <-
  function(seurat_obj,
           ref_sce = NULL,
           transfer_label = "",
           de_method = "wilcox") {
    if (!(class(seurat_obj[["RNA"]]) == "Assay5")) {
      seurat_obj <-
        scCustomize::Convert_Assay(
          seurat_obj = seurat_obj,
          convert_to = "V5"
        )
    }
    raw_seurat <-
      CreateSeuratObject(
        counts = seurat_obj@assays$RNA@layers$counts,
        assay = "RNA",
        meta.data = seurat_obj@meta.data
      )
    colnames(raw_seurat) <- colnames(seurat_obj)
    rownames(raw_seurat) <- rownames(seurat_obj)
    raw_seurat$barcodes <- colnames(seurat_obj)
    raw_sce <-
      as.SingleCellExperiment(raw_seurat)
    raw_sce <-
      logNormCounts(raw_sce)

    pred.label <-
      SingleR::SingleR(
        test = raw_sce,
        ref = ref_sce,
        labels = ref_sce[[transfer_label]],
        de.method = de_method,
        num.threads = 10
      )
    return(pred.label)
  }

sce_predicted_label <-
  singleR_return_annot_from_scratch(
    seurat_obj = seurat_obj_2_plot,
    ref_sce = sce_mouse_liver_ref,
    transfer_label = "cell_type__ontology_label",
    de_method = "wilcox"
  )

colnames(seurat_obj_2_plot)
sum(duplicated(str_split(colnames(seurat_obj_2_plot),
  pattern = "_", simplify = T
)[, 1]))
ncol(seurat_obj_2_plot)
rownames(sce_predicted_label)
df_barcode_lookup <-
  data.frame(
    barcode = rownames(sce_predicted_label),
    predicted_label = sce_predicted_label$pruned.labels
  )

all(df_barcode_lookup$barcode == colnames(seurat_obj_2_plot))
seurat_obj_2_plot$singleR.predicted.ident <- df_barcode_lookup$predicted_label

# singleR labels
cell_ontology_palette <-
  DiscretePalette(
    n = length(unique(mouse_liver_ref$cell_type__ontology_label)),
    palette = "alphabet2"
  )
names(cell_ontology_palette) <-
  unique(mouse_liver_ref$cell_type__ontology_label)

predicted.cell_type__ontology_label <-
  unique(seurat_obj_2_plot$singleR.predicted.ident)
predicted.cell_type__ontology_palette <-
  cell_ontology_palette[names(cell_ontology_palette) %in% predicted.cell_type__ontology_label]
names(predicted.cell_type__ontology_palette)
predicted.cell_type__ontology_palette <-
  predicted.cell_type__ontology_palette[order(names(predicted.cell_type__ontology_palette))]

scCustomize::DimPlot_scCustom(seurat_obj_2_plot,
  reduction = "umap.harmony",
  group.by = "singleR.predicted.ident",
  shuffle = T,
  alpha = 0.5,
  colors_use = predicted.cell_type__ontology_palette,
  split_seurat = F,
  add_prop_plot = F
)
