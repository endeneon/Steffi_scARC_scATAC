library(readxl)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

setwd("~/GREAT/")
excel_file <- "Redo_depth_cutoff_Table_S24. ASoC SNPs in each cell type and time point.xlsx"
sheets <- excel_sheets(excel_file)
data_list <- lapply(sheets, function(s) {
  read_excel(excel_file, sheet = s)
})

names(data_list) <- sheets

filtered_list <- lapply(data_list, function(df){
  df %>% 
    mutate(FDR = as.numeric(FDR)) %>%  
    filter(FDR < 0.05) 
})


ASoC_intersect_sum_dir <-
  "ASoC_25_ways_sum_dir"
dir.create(path = ASoC_intersect_sum_dir,
           recursive = T)

## make a list that holds all ASoC SNPs
names(filtered_list)
cell_type_list <-
  c("GABA", "nmglut", "npglut")
cell_time_list <-
  c("0h", "1h", "6h")
## make a list that holds all ASoC SNPs
cell_type_list <- c("GABA", "nmglut", "npglut")
cell_time_list <- c("0h", "1h", "6h")

for (i in seq_along(cell_type_list)) {
  
  cell_type <- cell_type_list[i]
  
  ## collect ASoC SNPs per time
  all_ASoC_SNPs <- vector("list", length(cell_time_list))
  names(all_ASoC_SNPs) <- cell_time_list
  
  for (time in cell_time_list) {
    key <- paste(cell_type, time, "ASoC", sep = "_")
    all_ASoC_SNPs[[time]] <- filtered_list[[key]]
    print(key)
  }
  
  ## -----------------------------
  ## time-specific SNPs (0h / 1h / 6h)
  ## -----------------------------
  for (time in cell_time_list) {
    
    other_times <- setdiff(cell_time_list, time)
    
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[[time]][
        !all_ASoC_SNPs[[time]]$ID %in%
          unique(unlist(lapply(other_times, function(t)
            all_ASoC_SNPs[[t]]$ID))),
      ]
    
    print(paste(cell_type, time,
                "specific SNPs:", nrow(ASoC_stage_specific_output)))
    
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- "."
    
    write.table(
      table_writeout,
      file = paste0(
        ASoC_intersect_sum_dir, "/",
        cell_type, "_", time, "_specific.txt"
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  ## -----------------------------
  ## shared across all times (0h ∩ 1h ∩ 6h)
  ## -----------------------------
  {
    shared_ids <- Reduce(
      intersect,
      lapply(cell_time_list, function(t)
        all_ASoC_SNPs[[t]]$ID)
    )
    
    ASoC_shared_output <-
      all_ASoC_SNPs[[cell_time_list[1]]][
        all_ASoC_SNPs[[cell_time_list[1]]]$ID %in% shared_ids,
      ]
    
    print(paste(cell_type,
                "shared across times:",
                nrow(ASoC_shared_output)))
    
    table_writeout <-
      ASoC_shared_output[, c(1, 2, 2, 6)]
    
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- "."
    
    write.table(
      table_writeout,
      file = paste0(
        ASoC_intersect_sum_dir, "/",
        cell_type, "_shared.txt"
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  ## -----------------------------
  ## at least one time (0h ∪ 1h ∪ 6h)
  ## -----------------------------
  {
    ASoC_at_least_one_output <-
      do.call(rbind, all_ASoC_SNPs)
    
    ASoC_at_least_one_output <-
      ASoC_at_least_one_output[
        !duplicated(ASoC_at_least_one_output$ID),
      ]
    
    print(paste(cell_type,
                "at least one time:",
                nrow(ASoC_at_least_one_output)))
    
    table_writeout <-
      ASoC_at_least_one_output[, c(1, 2, 2, 6)]
    
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- "."
    
    write.table(
      table_writeout,
      file = paste0(
        ASoC_intersect_sum_dir, "/",
        cell_type, "_at_least_one.txt"
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

## =========================================================
## across cell types (per time)
## =========================================================

for (time in cell_time_list) {
  
  df_GABA   <- filtered_list[[paste("GABA",   time, "ASoC", sep = "_")]]
  df_nmglut <- filtered_list[[paste("nmglut", time, "ASoC", sep = "_")]]
  df_npglut <- filtered_list[[paste("npglut", time, "ASoC", sep = "_")]]
  
  ## shared across all cell types
  shared_ids <-
    Reduce(intersect,
           list(df_GABA$ID,
                df_nmglut$ID,
                df_npglut$ID))
  
  ASoC_shared_cells_output <-
    df_GABA[df_GABA$ID %in% shared_ids, ]
  
  print(paste("shared across cell types", time,
              nrow(ASoC_shared_cells_output)))
  
  table_writeout <-
    ASoC_shared_cells_output[, c(1, 2, 2, 6)]
  
  table_writeout[[2]] <-
    as.numeric(table_writeout[[2]]) - 1
  
  table_writeout$SCORE <- 100
  table_writeout$STRAND <- "."
  
  write.table(
    table_writeout,
    file = paste0(
      ASoC_intersect_sum_dir, "/",
      "sharedByAllCellTypes_", time, ".txt"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  ## at least one cell type
  ASoC_any_cell_output <-
    rbind(df_GABA, df_nmglut, df_npglut)
  
  ASoC_any_cell_output <-
    ASoC_any_cell_output[
      !duplicated(ASoC_any_cell_output$ID),
    ]
  
  print(paste("at least one cell type", time,
              nrow(ASoC_any_cell_output)))
  
  table_writeout <-
    ASoC_any_cell_output[, c(1, 2, 2, 6)]
  
  table_writeout[[2]] <-
    as.numeric(table_writeout[[2]]) - 1
  
  table_writeout$SCORE <- 100
  table_writeout$STRAND <- "."
  
  write.table(
    table_writeout,
    file = paste0(
      ASoC_intersect_sum_dir, "/",
      "atLeastInOneCellType_", time, ".txt"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}
