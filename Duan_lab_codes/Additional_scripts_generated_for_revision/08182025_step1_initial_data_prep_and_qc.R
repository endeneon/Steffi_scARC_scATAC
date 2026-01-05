############8/18/2025#####
# scRNA_10X_ASD_Science_08112025
# There are 15 datasets: for 4 Genes (CPT1C, MEF2C, RORB) and One wildtype * 3 stimulation conditions (0hr, 1hr, 6hrs)
# Each dataset is pooled. 


###Open the required libraries
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)

# read h5 ####
h5.list <- 
  list.files('./',
             'filtered_feature_bc_matrix.h5',
             recursive = TRUE,
             include.dirs = FALSE,
             full.names = TRUE)

h5.list <- h5.list[!grepl("Undetermined", h5.list)]


RNA.list <- list()
for (i in 1:length(h5.list)){
  h5.file <- Read10X_h5(filename = h5.list[i])
  fileextract <- strsplit(h5.list[i], "/")[[1]][7]
  pattern <- str_split(fileextract, pattern = "_")[[1]][1]
  cat("h5: ", pattern)
  rna.dat <- CreateSeuratObject(counts = h5.file,
                                assay = 'RNA')
  RNA.list[[pattern]] <- rna.dat
}
rm(rna.dat)

for (i in 1:length(RNA.list)) {
  library <- names(RNA.list)[i]
  RNA.list[[i]] <- PercentageFeatureSet(RNA.list[[i]],
                                        pattern = "^MT-",
                                        col.name = 'percent.mt', assay = 'RNA')
  q <- VlnPlot(RNA.list[[i]],features = c("nCount_RNA", "nFeature_RNA", 
                                          "percent.mt"), ncol=3)
  ggsave(filename = paste0(library, "_count_feature_count_percentmt_RNA.jpg"), 
         width = 10, height=10, plot=q) 
  
}

qc_summary <- data.frame()

for (name in names(RNA.list)) {
  obj <- RNA.list[[name]]
  
  n_cells_before <- ncol(obj)
  n_feature_before <- nrow(obj)
  
  obj <- subset(obj,
                subset = nCount_RNA > 300 &
                  nCount_RNA < 20000 &
                  nFeature_RNA> 1000 &
                  nFeature_RNA < 7500 & 
                  percent.mt < 20)
  
  n_cells_after <- ncol(obj)
  n_feature_after <- nrow(obj)
  
  qc_summary <- rbind(
    qc_summary, 
    data.frame(
      library = name, 
      cells_before = n_cells_before, 
      features_before = n_feature_before,
      cells_after = n_cells_after, 
      features_after = n_feature_after
    )
  )
  
  RNA.list[[name]] <- obj
}
# 
# write.csv(qc_summary, 
#           "Number_of_cells-and_features_before_and_after_qc.csv")

save(RNA.list, file = "All_qc_RNA_files.RData")
