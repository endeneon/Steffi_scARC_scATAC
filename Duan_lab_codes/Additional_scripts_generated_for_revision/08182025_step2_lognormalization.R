###Open the required libraries
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)


ElbowPlot <- "ElbowPlot"
FeaturePlot <- "FeaturePlot"
UMAPPlot <- "UMAPPlot"
folders_new <- c(ElbowPlot, 
                 FeaturePlot,
                 UMAPPlot)
parent_dir <- "./Use_log_normalization_and_harmony/"
load("All_qc_RNA_files.RData") ##qc removed files
transformed_lst <- list()
for (name in names(RNA.list)) {
  obj <- RNA.list[[name]]
  obj <- PercentageFeatureSet(obj, pattern = c("^MT-"), col.name = "percent.mt")
  obj <- NormalizeData(obj, 
                       normalization.method = "LogNormalize", 
                       scale.factor = 1000)
  obj <- FindVariableFeatures(obj, nfeatures = 8000)
  all.genes <- rownames(obj)
  obj <- ScaleData(obj,features = all.genes)
  obj <- RunPCA(obj) 
  elbow_plot <- ElbowPlot(obj, ndims = 50)
  ggsave(filename = file.path(ElbowPlot, 
                              paste0(name, 
                                     "_elbowPlot.jpg")),
         width = 10, 
         height=10, 
         plot=elbow_plot, bg = "white")
  obj <- FindNeighbors( obj,
                        reduction = 'pca',
                        dims = 1:30)
  obj <- FindClusters( obj,
                       resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                       random.seed = 42, algorithm = 4)
  obj <- RunUMAP(object =  obj,
                 reduction = 'pca',
                 dims = 1:30,
                 seed.use = 42)
  resolutions <- seq(0.1, 1.0, by = 0.1)
  for (res in resolutions) {
    Idents(obj) <- paste0("RNA_snn_res.", res)
    Res_plot <- DimPlot(obj, label=T)
    ggsave(filename = file.path(UMAPPlot,
                                paste0(name, 
                                       "_UMAPPlot_res", res, ".jpg")), 
           width = 10, 
           height=10, 
           plot= Res_plot)
  }
  Idents(obj) <- "RNA_snn_res.0.5"
  dim2 <- FeaturePlot(obj,
                      features = c("SLC17A6", "SLC17A7", 
                                   "GFAP", "S100B", 
                                   "NEFM", "GAD1"))
  ggsave(filename = file.path(FeaturePlot,
                              paste0(name, 
                                     "_FeaturePlot_Res0pt5.jpg")), 
         width = 10, 
         height=10, 
         plot=dim2, bg = "white") 
  transformed_lst[[name]] <- obj
}

for (name in names(transformed_lst)) {
  dim2 <- FeaturePlot(obj,
                      features = c("SLC17A6", "SLC17A7", 
                                   "GFAP", "S100B", 
                                   "NEFM", "GAD1"))
  ggsave(filename = file.path(FeaturePlot,
                              paste0(name, 
                                     "_FeaturePlot_Res0pt5.jpg")), 
         width = 10, 
         height=10, 
         plot=dim2, bg = "white") 
}

save(transformed_lst, file = "log_transformed_data.RData")
