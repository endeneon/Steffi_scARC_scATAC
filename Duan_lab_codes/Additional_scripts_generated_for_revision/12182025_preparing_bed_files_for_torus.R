library(readxl)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

setwd("~/TORUS/Calculate_new_fdr/")
chain <- import.chain("~/Data/Databases/Genomes/hg38/hg38ToHg19.over.chain")
excel_file <- "../Redo_depth_cutoff_Table_S24. ASoC SNPs in each cell type and time point.xlsx"
sheets <- excel_sheets(excel_file)
data_list <- lapply(sheets, function(s) {
  read_excel(excel_file, sheet = s)
})

names(data_list) <- sheets


filtered_list <- lapply(data_list, function(df){
  df %>% 
    mutate(
      pVal = as.numeric(pVal),
      FDR_new = p.adjust(pVal, method = "BH")) %>%  
    filter(FDR_new < 0.05) 
})

df_to_bed <- function(df) {
  
  # Create GRanges (hg38)
  gr_hg38 <- GRanges(
    seqnames = df$CHROM,
    ranges = IRanges(
      names = df$ID,
      start = as.numeric(df$POS),
      end   = as.numeric(df$POS)
    )
  )
  
  # LiftOver to hg19
  gr_hg19 <- liftOver(gr_hg38, chain)
  
  # Convert to data.frame
  result_hg19 <- as.data.frame(gr_hg19)
  result_hg19 <- result_hg19[, -1]
  colnames(result_hg19)[1:2] <- c("ID", "CHR")
  
  # Extract chr number
  index <- seq(2, nrow(result_hg19) * 2, by = 2)
  chr <- unlist(strsplit(as.character(result_hg19$CHR), "chr"))[index]
  
  # Create BED (0-based start)
  bed <- data.frame(
    chr,
    start = as.numeric(result_hg19$start) - 1,
    end   = as.numeric(result_hg19$end)
  )
  
  return(bed)
}

bed_list <- lapply(filtered_list, df_to_bed)
output_dir <- "./bed_files"
dir.create(output_dir, showWarnings = FALSE)

mapply(function(bed, nm) {
  write.table(
    bed,
    file = file.path(output_dir, paste0(nm, ".bed")),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}, bed_list, names(bed_list))
