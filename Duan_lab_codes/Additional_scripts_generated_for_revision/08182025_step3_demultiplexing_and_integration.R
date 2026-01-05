###Open the required libraries
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(DoubletFinder)
library(openxlsx)
library(qs)

setwd("./Use_log_normalization_and_harmony/")
load("log_transformed_data.RData")

files_path <- "./barcode_collection"
pathlist <- sort(list.files(path = files_path, full.names = T, pattern = "SNG.barcodes", recursive = T))
barcode_list <- vector(mode = "list", length = 16)
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
  split_paths <- strsplit(basename(pathlist[[i]]), "_")[[1]]
  names(barcode_list)[i] <- split_paths[[1]]
}

barcode_list <- barcode_list[names(barcode_list) != "Undetermined"]

for (name in names(transformed_lst)) {
  obj <- transformed_lst[[name]]
  seurat_barcodes <- colnames(obj@assays$RNA$counts)
  obj$cell_line <- "unassigned"
  barcode_df <- barcode_list[[name]]
  matched_indices <- match(seurat_barcodes, barcode_df$barcode)
  obj$cell_line <- barcode_df$line[matched_indices]
  obj$cell_line[obj$cell_line == "V4500_1"] <- "KOLF2.2J"
  gene <- str_extract(name, "^[A-Za-z]+")
  time <- str_extract(name, "[0-9]+$")
  obj$library <- name
  obj$lof_gene <- gene
  obj$time_point <- time
  transformed_lst[[name]] <- obj
}

save(transformed_lst, file =  "Demultiplexed_files.RData")

load("Demultiplexed_files.RData")


doublet_results <- list()

for (name in names(transformed_lst))  {
  tryCatch({
    obj <- transformed_lst[[name]]
    
    sweep.res <- paramSweep(obj, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    
    nExp <- round(0.08 * nrow(obj@meta.data))
    
    obj <- doubletFinder(obj,
                         PCs  = 1:30, pN = 0.25,
                         pK = as.numeric(pK),
                         nExp = nExp,
                         reuse.pANN = NULL, sct = FALSE)
    transformed_lst[[name]] <- obj
    
    DF.name <- colnames(obj@meta.data)[grepl("DF.classification", colnames(obj@meta.data))]
    
    # Calculate doublet counts and percentage
    counts <- table(obj@meta.data[[DF.name]])
    percent_doublets <- if ("Doublet" %in% names(counts)) {
      counts["Doublet"] / sum(counts) * 100
    } else {
      0 # If no doublets are detected
    }
    
    # Store results
    doublet_results[[name]] <- list(
      counts = counts,
      percent_doublets = percent_doublets
    )
    
    cat("Processed sample:", name, "\n")
    
  }, error = function(e) {
    cat("Error in sample:", name, "\nMessage:", conditionMessage(e), "\n")
    doublet_results[[name]] <- list(error = conditionMessage(e))
  })
}
qsave(transformed_lst, "Doublet_finder_demultiplexed.qs")

doublet_summary <- do.call(rbind, lapply(names(doublet_results), function(name){
  res <- doublet_results[[name]]
  
  if (!is.null(res$error)){
    return(data.frame(Sample = name,
                      Singlets = NA,
                      Doublets = NA,
                      PercentDoublets = NA,
                      Status = paste("Error:", res$error)))
  }
  
  counts <- res$counts
  singlets <- if ("Singlet" %in% names(counts)) counts[["Singlet"]] else 0
  doublets <- if ("Doublet" %in% names(counts)) counts[["Doublet"]] else 0
  
  data.frame(Sample = name,
             Singlets = singlets,
             Doublets = doublets,
             PercentDoublets = round(res$percent_doublets, 2),
             Status = "OK")
}))



save_path <- "./Use_log_normalization_and_harmony/Doubletfinder/"
summary_list <- list()
transformed_lst_filtered <- list()
for (name in names(transformed_lst)){
  res <- transformed_lst[[name]]
  cls_col <- grep("^DF\\.classifications", colnames(res@meta.data), value = TRUE)[1]
  colnames(res@meta.data)[colnames(res@meta.data) == cls_col] <- "doublet_finder"
  plot_umap <- UMAPPlot(res, label = TRUE, group.by = "doublet_finder")
  print(plot_umap)
  ggsave(plot= plot_umap, filename = paste0(save_path,"Inital_umap_group_by_doublets_", name, ".pdf"),
         width = 10, height=10)
  
  plot_umap2 <- UMAPPlot(res, label = TRUE, group.by = "cell_line")
  print(plot_umap2)
  ggsave(plot= plot_umap2, filename = paste0(save_path, "Inital_umap_group_by_cell_line_", name, ".pdf"),
         width = 10, height=10)
  prep_summary <- res@meta.data %>% 
    group_by(cell_line, doublet_finder) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(sample = name, stage = "before")
  filtered_obj <-
    subset(res, subset = doublet_finder == "Singlet" &
             cell_line %in% c("CW20107", "KOLF2.2J")
    )
  
  post_summary <- filtered_obj@meta.data %>% 
    group_by(cell_line, doublet_finder) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(sample = name, stage = "before")
  
  plot_umap <- UMAPPlot(filtered_obj, label = TRUE, group.by = "doublet_finder")
  print(plot_umap)
  ggsave(plot= plot_umap, filename = paste0(save_path,"Filtered_umap_group_by_doublets_", name, ".pdf"),
         width = 10, height=10)
  
  plot_umap2 <- UMAPPlot(filtered_obj, label = TRUE, group.by = "cell_line")
  print(plot_umap2)
  ggsave(plot= plot_umap2, filename = paste0(save_path,"Filtered_umap_group_by_cell_line_", name, ".pdf"),
         width = 10, height=10)
  summary_list[[name]] <- bind_rows(prep_summary, post_summary)
  transformed_lst_filtered[[name]] <- filtered_obj
}
write.csv(doublet_summary, paste0(save_path, "Doublet_summary.csv"))
gc()

qsave(transformed_lst_filtered, "Singlets_subset_cellline_unintegrated.qs")


# wb <- createWorkbook()
# for (name in names(summary_list)){
#   addWorksheet(wb , name)
#   writeData(wb, name, summary_list[[name]])
# }
# 
# saveWorkbook(wb, paste0(save_path, "Doublet_summary_by_sample.xlsx"), overwrite = TRUE)
# save(doublet_results, file = paste0(save_path, "Doublet_results.RData"))
# 




merged_data2 <- merge(x = transformed_lst_filtered[[1]],
                      y = transformed_lst_filtered[2:length(transformed_lst_filtered)],)
merged_seurat <- merged_data2
print(merged_seurat)
qsave(merged_seurat, "merged_non-integrated.qs")
merged_seurat <- qread("merged_non-integrated.qs")

####Analyze without integrating
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent.mt")
merged_seurat <- NormalizeData(merged_seurat) 

merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 3000) 
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, reduction = "pca")
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5, cluster.name ="unintegrated_clusters")
merged_seurat <- RunUMAP(merged_seurat, dims=1:30, reduction = "pca", reduction.name = "umap.unintegrated")

UninUMAP <- DimPlot(merged_seurat, reduction = "umap.unintegrated", group.by = c("seurat_clusters", "library"), alpha = 0.5, raster = F)
ggsave(filename = "Unintegrated_UMAP_library.jpg", width = 20, height = 10, plot= UninUMAP)
UninUMAP2 <- DimPlot(merged_seurat, reduction = "umap.unintegrated", group.by = "library", raster = F)
ggsave(filename = "Unintegrated_UMAP_groupbylibrary.jpg", width = 20, height = 10, plot= UninUMAP2)
qsave(merged_seurat, "merged_non-integrated_2.qs")
#merged_seurat <- readRDS("merged_non-integrated_2.Rds")
#######Performing integration

merged_seurat <- IntegrateLayers(
  merged_seurat, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose=F,
  
)

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])

merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims=1:30)

merged_seurat <- FindClusters(merged_seurat, resolution = 0.1, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.2, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.3, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.4, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.6, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.7, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.8, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.9, algorithm = 4, random.seed = 1)
merged_seurat <- FindClusters(merged_seurat, resolution = 1, algorithm = 4, random.seed = 1)
qsave(merged_seurat, "Seurat_object_integrated.qs")

Idents(merged_seurat) <- "RNA_snn_res.1"
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
InUMAP_01 <- DimPlot(merged_seurat, reduction = "umap", group.by = c("library"), alpha = 0.5)
ggsave(filename = "Integrated_UMAP_groupbylibrary.jpg", width = 10, height = 10, plot= InUMAP_01)
InUMAP_01 <- DimPlot(merged_seurat, reduction = "umap", alpha = 0.5)
ggsave(filename = "Integrated_UMAP.jpg", width = 10, height = 10, plot= InUMAP_01)
sp <- Stacked_VlnPlot(merged_seurat, 
                      features=c("MAP2", "GAD1", "GAD2", "SLC17A6",
                                 "NEFM", "VIM", "OLIG2",
                                 "GFAP", "S100B", 
                                 "EBF1", "SEMA3E",
                                 "BCL11B", "SST", "SATB2", "SOX2",
                                 "SLC17A7", "SERTAD4", "FOXG1", "POU3F2", "LHX2",
                                 "ADCYAP1", 
                                 "CUX1", "CUX2", "DCX"))
ggsave(filename = paste0("StackedVlnPlot_integrated.jpg"),
       width = 10, 
       height = 10,
       plot = sp)
FeaturePlot(merged_seurat, features = c("MAP2", "GAD1", "GAD2", "SLC17A6",
                                        "NEFM", 
                                        "GFAP"), reduction = "umap")

Markers <- FindAllMarkers(merged_seurat, only.pos = TRUE)
Markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
new.cluster.id <- c("nmglut", "GABA", "GABA", "npglut", "nmglut", 
                    "nmglut", "npglut", "nmglut", "GABA", "nmglut",
                    "nmglut", "nmglut", "Astrocyte", "GABA", "npglut", 
                    "npglut", "nmglut","GABA", "GABA", "GABA", 
                    "glut ?", "Astrocyte", "unidentified", "unidentified", "glut ?",
                    "doublet", "glut ?")
length(new.cluster.id)
length(unique(merged_seurat$seurat_clusters))
names(new.cluster.id) <- levels(merged_seurat)
merged_seurat_labeled <- RenameIdents(merged_seurat, new.cluster.id)


DimPlot(merged_seurat_labeled, reduction = "umap", alpha = 0.5, label = T, raster = FALSE)


qsave(merged_seurat_labeled, "merged_integrated_labeled.qs")
