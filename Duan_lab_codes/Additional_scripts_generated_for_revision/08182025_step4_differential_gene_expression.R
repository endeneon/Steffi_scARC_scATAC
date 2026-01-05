library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(scCustomize)

setwd("/data/Christina_Thapa/ASD_gene_network_scRNA/Use_log_normalization_and_harmony/")

save_path <- "./Use_log_normalization_and_harmony/Using_Cellline_As_covar/"

merged_seurat <- qread("merged_integrated_labeled_cleaned.qs")

Idents(merged_seurat) <- "lof_gene"
Idents(merged_seurat) <- "celltype"
new.cluster.id <- c("glut", "GABA", "glut")
length(new.cluster.id)
length(unique(merged_seurat$celltype))
names(new.cluster.id) <- levels(merged_seurat)
merged_seurat_labeled <- RenameIdents(merged_seurat, new.cluster.id)
merged_seurat_labeled$broad_celltype <- merged_seurat_labeled@active.ident

time_points <- unique(merged_seurat_labeled$time_point)
KD_gene <- unique(merged_seurat_labeled$lof_gene)
KD_gene <- KD_gene[KD_gene != 'WT']
line <- unique(merged_seurat_labeled$cell_line)
celltypes <- unique(merged_seurat_labeled$broad_celltype)
libraries <- unique(merged_seurat_labeled$library)

Idents(merged_seurat_labeled) <- "lof_gene"

all_results <- list()


for(ct in celltypes){
  for (gene in KD_gene){
    for (tp in time_points){
      
      subset_obj <- subset(merged_seurat_labeled, 
                           broad_celltype == ct  & time_point == tp)
      subset_obj <- subset(subset_obj, lof_gene %in% c('WT', gene))
      control_cells <- WhichCells(
        subset_obj,
        expression = lof_gene == "WT" 
      )
      gene_cells <- WhichCells(
        subset_obj, 
        expression = lof_gene == gene 
      )
      cells_to_keep <- subset_obj@meta.data %>%
        tibble::rownames_to_column("cell_id") %>%
        group_by(cell_line, library, lof_gene) %>%
        sample_n(size = 630, replace = FALSE) %>%
        pull(cell_id)
      subset_obj <- subset(subset_obj, cells = cells_to_keep)
      
      markers <- FindMarkers(
        subset_obj, 
        ident.1 = gene,
        ident.2 = 'WT',
        logfc.threshold = 0,
        min.pct = 0.01,
        test.use = "MAST",
        latent.vars = "cell_line"
      )
      all_results[[paste(ct, tp, gene, sep = "_")]] <- markers
    }
  }
}




filtered_results <- list()
for (name in names(all_results)){
  df <- all_results[[name]]
  df_filter <- df[abs(df$avg_log2FC) > 0 & df$p_val_adj < 0.05, ]
  df$regulation <- ifelse(df$p_val_adj < 0.05 &df$avg_log2FC > 0, "up", 
                          ifelse(df$p_val_adj < 0.05 &df$avg_log2FC < 0, "down", "ns"))
  all_results[[name]] <- df
  df_filter <- df[abs(df$avg_log2FC) > 0 & df$p_val_adj < 0.05, ]
  filtered_results[[name]] <- df_filter
}

for (name in names(all_results)){
  df <- all_results[[name]]
  write.csv(df, paste0(save_path, name,".csv"), row.names = TRUE)
}

for (name in names(filtered_results)){
  df <- filtered_results[[name]]
  write.csv(df, paste0(save_path, name,"_significant.csv"), row.names = TRUE)
}

save(all_results, file = paste0(save_path, "All_combined_comparison.RData"))
save(filtered_results, file = paste0(save_path, "All_combined_comparison_significant.RData"))
