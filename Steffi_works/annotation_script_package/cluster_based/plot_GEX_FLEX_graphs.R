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
  "/research_jude/rgs01_jude/groups/cab/projects/automapper/common/szhang37/pulled_git_repos/Multiome_main/Steffi_works/annotation_script_package/cluster_based"
)
# determine if R is running in RSTUDIO/VSCode/Positron
if (
  interactive() &&
    ((Sys.getenv("RSTUDIO") == "1" || (Sys.getenv("TERM_PROGRAM") == "vscode")))
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
  systemfonts::register_font(
    "sans",
    plain = systemfonts::match_fonts("DejaVu Sans")$path
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
# merged_liver_obj <-
#   qs_read("./annotation_script_package/marker_based/step1.qs2", nthreads = 8)
merged_liver_obj <-
  qs_read(
    "../marker_based/annotated_seurat_marker_based.qs2",
    nthreads = 8
  )

scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "celltype_broad",
  split.by = "celltype_broad",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = T,
  raster = F,
  # DimPlot_scCustom's split.by branch does `names(colors_overall) <- levels_overall`
  # on `colors_use` as-is, without expanding a palette *name* into colors first
  # (that expansion only happens when colors_use is left NULL). Passing the
  # string "glasbey" therefore fails with the length mismatch error; pass an
  # explicit color vector sized to the number of celltype_broad levels instead.
  colors_use = DiscretePalette_scCustomize(
    num_colors = length(levels(as.factor(merged_liver_obj$celltype_broad))),
    palette = "glasbey"
  )
)

scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "SCT_snn_res.0.25",
  split.by = "SCT_snn_res.0.25",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = T,
  raster = F,
  # DimPlot_scCustom's split.by branch does `names(colors_overall) <- levels_overall`
  # on `colors_use` as-is, without expanding a palette *name* into colors first
  # (that expansion only happens when colors_use is left NULL). Passing the
  # string "glasbey" therefore fails with the length mismatch error; pass an
  # explicit color vector sized to the number of SCT_snn_res.0.25 levels instead.
  colors_use = DiscretePalette_scCustomize(
    num_colors = length(levels(as.factor(merged_liver_obj$SCT_snn_res.0.25))),
    palette = "varibow"
  )
)

# `CellSelector()`/`HoverLocator()` require an interactive graphics device
# (X11/RStudio viewer). No R build on this cluster (conda or module, incl.
# R/4.3.0-4.4.0) is compiled with X11 support, so `X11()` always fails with
# "unable to open connection to X11 display" here. Replace the interactive
# workflow with:
#   1. An interactive `plotly` UMAP (an HTML widget rendered in the
#      Positron/RStudio Viewer pane -- no X11 device needed) to visually
#      identify cells/regions and read off their umap coordinates by hovering.
#   2. A non-interactive polygon gate (`select_cells_by_polygon()`) that
#      reproduces what CellSelector's lasso/box selection did: it flags cells
#      inside a user-supplied polygon and writes the label into
#      `manual_annotation`, exactly like `CellSelector(..., ident = "manual_annotation")`.
# umap_coords <-
#   as.data.frame(Embeddings(merged_liver_obj, reduction = "umap.harmony")[, 1:2])
# colnames(umap_coords) <- c("umap_1", "umap_2")
# umap_coords$cell_barcode <- rownames(umap_coords)
# umap_coords$celltype_broad <- merged_liver_obj$celltype_broad

# umap_plotly <-
#   plotly::plot_ly(
#     umap_coords,
#     x = ~umap_1,
#     y = ~umap_2,
#     color = ~celltype_broad,
#     colors = DiscretePalette_scCustomize(
#       num_colors = length(levels(as.factor(umap_coords$celltype_broad))),
#       palette = "glasbey"
#     ),
#     text = ~cell_barcode,
#     hoverinfo = "text",
#     type = "scatter",
#     mode = "markers",
#     marker = list(size = 4, opacity = 0.7)
#   )
# umap_plotly

# `plotly::event_data("plotly_selected")` only returns data inside a reactive
# context, so capturing a lasso/box selection straight back into R requires a
# small Shiny gadget -- the same mechanism Seurat's own `CellSelector()` uses
# (miniUI + shiny::runGadget()), just with a `plotly` widget instead of a base
# R plot. Since it renders as an HTML widget (Viewer pane/browser), this needs
# no X11 device. Run this yourself in an interactive R console/Positron
# console (do not source/execute non-interactively -- it blocks until you
# click "Done").
lasso_select_cells <-
  function(
    object,
    reduction = "umap.harmony",
    group.by = NULL,
    ident_label = "SelectedCells",
    ident_col = "manual_annotation"
  ) {
    coords <- as.data.frame(Embeddings(object, reduction = reduction)[, 1:2])
    colnames(coords) <- c("umap_1", "umap_2")
    coords$cell_barcode <- rownames(coords)
    if (!is.null(group.by)) {
      coords$group_var <- as.factor(object@meta.data[[group.by]])
    }

    ui <-
      miniUI::miniPage(
        miniUI::gadgetTitleBar("Lasso- or box-select cells, then click Done"),
        miniUI::miniContentPanel(
          plotly::plotlyOutput("umap_plot", height = "100%")
        )
      )

    server <- function(input, output, session) {
      output$umap_plot <- plotly::renderPlotly({
        p <-
          plotly::plot_ly(
            coords,
            x = ~umap_1,
            y = ~umap_2,
            customdata = ~cell_barcode,
            color = if (!is.null(group.by)) ~group_var else I("steelblue"),
            colors = if (!is.null(group.by)) {
              DiscretePalette_scCustomize(
                num_colors = nlevels(coords$group_var),
                palette = "glasbey"
              )
            } else {
              NULL
            },
            type = "scatter",
            mode = "markers",
            marker = list(size = 4, opacity = 0.7),
            source = "lasso_umap"
          )
        plotly::layout(p, dragmode = "lasso")
      })

      selected_barcodes <- shiny::reactiveVal(character(0))

      shiny::observeEvent(
        plotly::event_data("plotly_selected", source = "lasso_umap"),
        {
          ed <- plotly::event_data("plotly_selected", source = "lasso_umap")
          if (!is.null(ed)) {
            selected_barcodes(as.character(ed$customdata))
          }
        }
      )

      shiny::observeEvent(input$done, {
        bc <- selected_barcodes()
        result <- object
        if (length(bc) > 0) {
          if (!ident_col %in% colnames(result@meta.data)) {
            result@meta.data[[ident_col]] <- "Unselected"
          }
          result@meta.data[bc, ident_col] <- ident_label
        }
        shiny::stopApp(
          returnValue = list(
            object = result,
            selected_cells = bc,
            selected_coords = coords[coords$cell_barcode %in% bc, ]
          )
        )
      })
    }

    shiny::runGadget(ui, server, viewer = shiny::paneViewer())
  }

# Example usage (run interactively; drag the lasso/box-select tool in the
# plotly toolbar over the region you want, then click "Done"):
lasso_result <-
  lasso_select_cells(
    merged_liver_obj,
    reduction = "umap.harmony",
    group.by = "celltype_broad",
    ident_label = "SelectedCells",
    ident_col = "manual_annotation"
  )
selected_cells_obj <- lasso_result$object
selected_cells_obj <-
  subset(
    selected_cells_obj,
    subset = manual_annotation == "SelectedCells"
  )
ncol(selected_cells_obj) # number of selected cells
# lasso_result$selected_cells # character vector of selected cell barcodes
# lasso_result$selected_coords # umap_1/umap_2/cell_barcode for selected cells

qs_save(
  selected_cells_obj,
  file = "manual_selection_hepatocytes.qs2",
  nthreads = 8
)
selected_cells_obj <-
  qs_read(
    "manual_selection_hepatocytes.qs2",
    nthreads = 8
  )
# qs_save(
#   selected_cells_obj,
#   file = "manual_selection_macrophages.qs2",
#   nthreads = 8
# )
# Alternative: if you already know the polygon vertices (e.g. read off by
# hovering on `umap_plotly` above), gate cells non-interactively instead.
# select_cells_by_polygon <-
#   function(
#     object,
#     reduction,
#     poly_x,
#     poly_y,
#     ident_label = "SelectedCells",
#     ident_col = "manual_annotation"
#   ) {
#     coords <- Embeddings(object, reduction = reduction)[, 1:2]
#     inside <-
#       sp::point.in.polygon(
#         point.x = coords[, 1],
#         point.y = coords[, 2],
#         pol.x = poly_x,
#         pol.y = poly_y
#       ) >
#         0
#     if (!ident_col %in% colnames(object@meta.data)) {
#       object@meta.data[[ident_col]] <- "Unselected"
#     }
#     object@meta.data[[ident_col]][inside] <- ident_label
#     object
#   }

# # Example usage (fill in vertices traced from `umap_plotly`):
# selected_obj <-
#   select_cells_by_polygon(
#     merged_liver_obj,
#     reduction = "umap.harmony",
#     poly_x = c(-5, -2, -2, -5),
#     poly_y = c(2, 2, 6, 6),
#     ident_label = "SelectedCells",
#     ident_col = "manual_annotation"
#   )
merged_liver_obj$hepatocyte <- FALSE
merged_liver_obj$hepatocyte[
  colnames(merged_liver_obj) %in% colnames(selected_cells_obj)
] <- TRUE
qs_save(
  merged_liver_obj, 
  file = "merged_liver_obj_with_hepatocyte.qs2", 
  nthreads = 8)
# merged_subsetted_obj <- selected_cells_obj
scCustomize::DimPlot_scCustom(
  merged_liver_obj,
  reduction = "umap.harmony",
  group.by = "hepatocyte",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = T,
  raster = F,
  colors_use = c("steelblue", "orange")
)

stesstmerged_subsetted_obj <- selected_cells_obj
scCustomize::DimPlot_scCustom(
  merged_subsetted_obj,
  reduction = "umap.harmony",
  group.by = "manual_annotation",
  label = F,
  label.size = 4,
  pt.size = 0.2,
  alpha = 1,
  shuffle = T,
  raster = F,
  colors_use = c("steelblue", "orange")
)
df_agg <-
  AggregateExpression(
    merged_subsetted_obj,
    assays = "RNA",
    slot = "counts",
    group.by = c("orig.ident"),
    verbose = T,
    return_seurat = F
  )
df_agg <- df_agg$RNA
# rownames(df_agg)
# colnames(df_agg)
df_meta <-
  merged_liver_obj@meta.data
# df_meta <-
#   df_meta[!duplicated(df_meta$orig.ident), ]
df_meta_ordered <-
  df_meta[
    match(
      colnames(df_agg),
      gsub("_", "-", as.character(df_meta$orig.ident))
    ),
  ]
# df_meta_ordered$category <- "Primary"
# df_meta_ordered$category[str_detect(
#   df_meta_ordered$orig.ident,
#   "pre"
# )] <- "Resistant"
# df_meta_ordered$category[
#   df_meta_ordered$orig.ident %in%
#     c(
#       "BC010",
#       "BC004",
#       "BC001",
#       "BC009",
#       "BC014",
#       "BC015",
#       "BC002",
#       "BC003",
#       "BC005",
#       "BC007",
#       "BC011",
#       "BC006",
#       "BC016",
#       "BC008"
#     )
# ] <- "Resistant"

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

sum(rowSums(DGE_cpm > 1) > 10)
DGE_agg <-
  DGEList(
    counts = df_agg[rowSums(DGE_cpm > 1) > 10, ],
    samples = df_meta_ordered$orig.ident,
    group = df_meta_ordered$category
  )
DGE_agg <-
  calcNormFactors(
    DGE_agg,
    method = "TMM"
  )

DGE_agg$samples$category <- df_meta_ordered$category
DGE_agg$samples$preparation <- df_meta_ordered$preparation

design_matrix <-
  model.matrix(
    ~ 0 + category + preparation,
    data = DGE_agg$samples
  )
# colnames(design_matrix) <-
#   gsub(
#     "group",
#     "",
#     colnames(design_matrix)
#   )

DGE_agg <-
  estimateDisp(
    DGE_agg,
    robust = TRUE,
    design = design_matrix
  )
DGE_agg <-
  glmQLFit(
    DGE_agg,
    design = design_matrix,
    robust = TRUE
  )
qlf_test <-
  glmQLFTest(
    DGE_agg,
    contrast = c(-1, 1, 0, 0) # use Primary as reference
  )

genes_2_plot <-
  read.table(
    "../../sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated_HCC_crossref_DisGeNET_yes.tsv",
    header = T,
    sep = "\t"
  )
genes_2_plot <- genes_2_plot$SYMBOL
genes_2_plot <- genes_2_plot[!duplicated(genes_2_plot)]

selected_gene_exp_table <-
  qlf_test$table
selected_gene_exp_table$FDR <-
  p.adjust(
    selected_gene_exp_table$PValue,
    method = "BH"
  )
selected_gene_exp_table$gene_name <-
  rownames(selected_gene_exp_table)
selected_gene_exp_table <-
  selected_gene_exp_table[order(selected_gene_exp_table$PValue), ]
selected_gene_exp_table <-
  selected_gene_exp_table[selected_gene_exp_table$gene_name %in% genes_2_plot, ]

if (!dir.exists("DE_results_by_celltype")) {
  dir.create("DE_results_by_celltype")
}
write.table(
  selected_gene_exp_table,
  file = "DE_results_by_celltype/DE_logCPM_t_test_reference_based_Hepatocytes.tsv",
  sep = "\t",
  quote = F,
  row.names = F
)
# write.table(
#   selected_gene_exp_table,
#   file = "DE_results_by_celltype/DE_logCPM_t_test_reference_based_Macrophages.tsv",
#   sep = "\t",
#   quote = F,
#   row.names = F
# )

df_2_plot <-
  DGE_cpm[rownames(DGE_cpm) %in% genes_2_plot, ]
df_2_plot <-
  as.data.frame(t(df_2_plot))
# rownames(df_2_plot)
# colnames(df_2_plot)

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

df_long$preparation <- "10xGEX"
df_long$preparation[str_detect(
  df_long$sample_name,
  "^BC"
)] <- "10xFLEX"
df_long$preparation <- as.factor(df_long$preparation)
genes_2_plot <- unique(df_long$gene_name)

DE_plot_list <-
  lapply(
    genes_2_plot,
    function(gene) {
      df_sub <-
        df_long[df_long$gene_name == gene, ]

      # Use the PValue already computed by the edgeR glmQLFTest contrast
      # (selected_gene_exp_table), matched by gene_name, instead of letting
      # ggpubr::stat_compare_means() silently recompute its own t-test.
      gene_pvalue <-
        selected_gene_exp_table$PValue[
          selected_gene_exp_table$gene_name == gene
        ]
      pvalue_stars <-
        as.character(
          stats::symnum(
            gene_pvalue,
            corr = FALSE,
            na = FALSE,
            cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
            symbols = c("****", "***", "**", "*", "ns")
          )
        )
      stat_df <-
        data.frame(
          group1 = "Primary",
          group2 = "Resistant",
          label = pvalue_stars,
          y.position = max(df_sub$exp_level) * 1.1
        )

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
        ggpubr::stat_pvalue_manual(
          stat_df,
          label = "label",
          y.position = "y.position"
        ) +
        expand_limits(y = 0) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
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

# if (!requireNamespace("httpgd", quietly = TRUE)) install.packages("httpgd")
# httpgd::hgd()          # opens a cairo/Skia device VS Code can display
# getOption("device"); names(dev.cur())
# options(device = function(...) X11(type = "cairo"))   # force cairo screen device

# print(
#   patchwork::wrap_plots(DE_plot_list, ncol = 4, guides = "collect") &
#     theme(legend.position = "right")
# )
# httpgd::hgd_browse()  # or open the printed URL if the pane doesn't auto-show

# `patchwork::plot_layout(axes = "collect_x")` does not merely hide the
# duplicated tick *labels* -- it removes the entire axis line/ticks from every
# panel except the outermost row of each column (patchwork assumes a
# facet-like shared grid). That leaves interior panels (e.g. HDAC11) with no
# visible x-axis at all. Instead, keep each panel's own axis line/ticks and
# manually blank just the tick *text* on every panel except the true
# bottom-most plot in each column, so "Primary"/"Resistant" labels are shown
# only once per column, at the bottom.
n_de_genes <- length(DE_plot_list)
de_ncol <- 4
de_col_idx <- ((seq_len(n_de_genes) - 1) %% de_ncol) + 1
de_bottom_indices <- sapply(seq_len(de_ncol), function(col) {
  max(which(de_col_idx == col))
})
for (i in seq_along(DE_plot_list)) {
  if (!(i %in% de_bottom_indices)) {
    DE_plot_list[[i]] <-
      DE_plot_list[[i]] + theme(axis.text.x = element_blank())
  }
}

print(
  patchwork::wrap_plots(DE_plot_list, ncol = de_ncol) +
    patchwork::plot_layout(
      guides = "collect",
      axis_titles = "collect"
    ) &
    theme(legend.position = "right")
)
ggsave(
  "DE_results_by_celltype/DE_logCPM_edgeR_test_cluster_based_Hepatocytes.pdf",
  width = 6,
  height = 12
)
# ggsave(
#   "DE_results_by_celltype/DE_logCPM_edgeR_test_cluster_based_Macrophages.pdf",
#   width = 6,
#   height = 12
# )
# what device is active?
# then open a fresh device and re-print

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
