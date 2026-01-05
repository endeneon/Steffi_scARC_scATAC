library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(scCustomize)


save_path <- "./Use_log_normalization_and_harmony/Prep_for_paper/"

merged_seurat <- qread("merged_integrated_labeled_cleaned.qs")


Idents(merged_seuratl) <- "lof_gene"
Idents(merged_seurat) <- "celltype"
new.cluster.id <- c("Glut", "GABA", "Glut")
length(new.cluster.id)
length(unique(merged_seurat_cpt_null$celltype))
names(new.cluster.id) <- levels(merged_seurat_cpt_null)
merged_seurat_labeled <- RenameIdents(merged_seurat_cpt_null, new.cluster.id)
merged_seurat_labeled$broad_celltype <- merged_seurat_labeled@active.ident

UMAPPlot(merged_seurat_labeled, reduction = "umap")

UMAPPlot(merged_seurat_labeled, reduction = "umap")

UMAPPlot(merged_seurat_labeled, reduction = "umap", group.by = "celltype")
UMAPPlot(merged_seurat_labeled, reduction = "umap", split.by = "lof_gene")
UMAPPlot(merged_seurat_labeled, reduction = "umap", split.by = "cell_line")
meta <- merged_seurat_labeled@meta.data

cell_counts <- meta %>%
  group_by(celltype, broad_celltype, library, lof_gene, cell_line) %>%
  summarise(n_cells = n()) %>%
  arrange(broad_celltype, library)

save(cell_counts, file = paste0(save_path, "Cell_counts_for_plotting.Rdata"))
write.csv(cell_counts, file = paste0(save_path, "Cell_counts_for_plotting.csv"))

celltypes <- unique(merged_seurat_labeled$broad_celltype)
line <- unique(merged_seurat_labeled$cell_line)
merged_seurat_labeled$library <- recode(
  merged_seurat_labeled$library, 
  "MEF0" = "MEF2C_0hr",
  "MEF1" = "MEF2C_1hr", 
  "MEF6" = "MEF2C_6hr", 
  "ROR0" = "RORB_0hr", 
  "ROR1" = "RORB_1hr", 
  "ROR6" = "RORB_6hr", 
  "TCF0" = "TCF4_0hr", 
  "TCF1" = "TCF4_1hr", 
  "TCF6" = "TCF4_6hr", 
  "WT0" = "WT_0hr", 
  "WT1"  = "WT_1hr", 
  "WT6" = "WT_6hr"
)

for (cl in line){
  for(ct in celltypes){
    subset_obj <- subset(merged_seurat_labeled, 
                         cell_line == cl &
                           broad_celltype == ct )
    lrg_erg_plot <- DotPlot(merged_seurat_labeled, 
                            features = c("BDNF", "VGF", 
                                         "FOS" , "NPAS4"),
                            group.by = "library",
                            cols = c("white", "red")) + 
      coord_flip() + 
      RotatedAxis() 
    
    ggsave(lrg_erg_plot, 
           filename = paste0(save_path, "Late_and_early_Response_genes", cl, "_", ct, ".pdf"),
           width = 6, height = 4, dpi = 300)
    ggsave(lrg_erg_plot, 
           filename = paste0(save_path, "Late_and_early_Response_genes", cl, "_", ct, ".png"),
           width = 6, height = 4, dpi = 300, bg = "white")
  }
}


