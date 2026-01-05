library(dplyr)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(readxl)
library(qvalue)
setwd(".pi1_estimation/MAST_leading/")

npglut_files <- list.files(path = ".",
                           pattern = "^revision_.*npglut_*", full.names = T)
nmglut_files <- list.files(path = ".",
                           pattern = "^revision_.*nmglut*", full.names = T)
get_timepoint <- function(filename) {
  # Look for 0hr, 1hr, or 6hr in the filename
  if (grepl("0hr", filename, ignore.case = TRUE)) return("0hr_limma")
  if (grepl("1hr", filename, ignore.case = TRUE)) return("1hr_limma")
  if (grepl("6hr", filename, ignore.case = TRUE)) return("6hr_limma")
  return(NA)  # If none found
}

npglut_list <- lapply(npglut_files, read.csv)
names(npglut_list) <- sapply(npglut_files, get_timepoint)
nmglut_list <- lapply(nmglut_files, read.csv)
names(nmglut_list) <- sapply(nmglut_files, get_timepoint)

# npglut_list <- lapply(npglut_list, function(df) {
#   df %>% dplyr::filter(P.Value < 0.05)
# })
# 
# nmglut_list <- lapply(nmglut_list, function(df) {
#   df %>% dplyr::filter(P.Value < 0.05)
# })



MAST_result <- read_excel("./MAST result for copmparison with Lexi's Pseudobulk.xlsx", 
                          col_types = c("text", "numeric", "skip", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "text", "skip", "skip", 
                                        "skip", "skip", "skip"))
MAST_result <- MAST_result[, !names(MAST_result) %in% "padj"]
names(MAST_result)[names(MAST_result) == "geneIDs"] <- "genes"
names(MAST_result)[names(MAST_result) == "log2FC"]  <- "logFC"
names(MAST_result)[names(MAST_result) == "BH_FDR"]    <- "P.Value"
# MAST_result <- MAST_result %>% 
#   filter(P.Value < 0.05)

nmglut_only <- subset(MAST_result, grepl("^nmglut", category))
nmglut_MAST_split <- split(nmglut_only,
                           sub(".*(0hr|1hr|6hr)$", "\\1", nmglut_only$category))
nmglut_list[paste0(names(nmglut_MAST_split), "_MAST")] <- nmglut_MAST_split

npglut_only <- subset(MAST_result, grepl("^npglut", category))
npglut_MAST_split <- split(npglut_only,
                           sub(".*(0hr|1hr|6hr)$", "\\1", npglut_only$category))
npglut_list[paste0(names(npglut_MAST_split), "_MAST")] <- npglut_MAST_split


# estimate_pi1 <- function(pvals) {
#   if (length(pvals) == 0) {
#     return(0)  # Or NA, depending on preference; handles no leading genes
#   }
#   pvals <- pmin(pmax(pvals, 1e-15), 1 - 1e-15)  # Still clip to avoid exact 0/1, though qvalue can handle
#   tryCatch({
#     qobj <- qvalue(p = pvals)  # Default uses smoother method for pi0
#     pi0 <- qobj$pi0
#     pi1 <- 1 - pi0
#     return(pi1)
#   }, error = function(e) {
#     warning("qvalue estimation failed; falling back to simple estimate")
#     lambda <- 0.5  # Fixed lambda as a robust simple alternative
#     pi0 <- mean(pvals > lambda) / (1 - lambda)
#     pi1 <- 1 - min(pi0, 1)
#     return(pi1)
#   })
# }
estimate_pi1 <- function(pvals) {
  if (length(pvals) == 0) return(0)
  
  lambda <- 0.5
  pi0 <- mean(pvals > lambda) / (1 - lambda)
  pi0 <- min(max(pi0, 0), 1)  # keep in [0,1]
  pi1 <- 1 - pi0
  return(pi1)
}

pi1_heatmap_matrix <- function(limma_list, mast_list, direction = c("up", "down")) {
  direction <- match.arg(direction)
  
  limma_timepoints <- c("0hr_limma", "1hr_limma", "6hr_limma")
  mast_timepoints <- c("0hr_MAST", "1hr_MAST", "6hr_MAST")
  
  pi1_mat <- matrix(NA, nrow = length(limma_timepoints), ncol = length(mast_timepoints),
                    dimnames = list(limma_timepoints, mast_timepoints))
  
  for (i in seq_along(limma_timepoints)) {
    for (j in seq_along(mast_timepoints)) {
      limma_df <- limma_list[[limma_timepoints[i]]]
      mast_df <- mast_list[[mast_timepoints[j]]]
      
      if (direction == "up") {
        lead_genes <- mast_df %>% filter(logFC > 0, P.Value < 0.05) %>% pull(genes)
      } else {
        lead_genes <- mast_df %>% filter(logFC < 0, P.Value < 0.05) %>% pull(genes)
      }
      
      pvals_matched <- limma_df %>% filter(genes %in% lead_genes) %>% pull(P.Value)
      pi1_mat[i, j] <- estimate_pi1(pvals_matched)
    }
  }
  
  return(pi1_mat)
}

# Example: npglut
pi1_up <- pi1_heatmap_matrix(npglut_list[1:3], npglut_list[4:6], direction = "up")
pi1_down <- pi1_heatmap_matrix(npglut_list[1:3], npglut_list[4:6], direction = "down")
pi1_nm_up <- pi1_heatmap_matrix(nmglut_list[1:3], nmglut_list[4:6], direction = "up")
pi1_nm_down <- pi1_heatmap_matrix(nmglut_list[1:3], nmglut_list[4:6], direction = "down")
# Convert to long format for ggplot
pi1_up_df <- as.data.frame(as.table(pi1_up)) %>% rename(Limma = Var1, MAST = Var2, pi1 = Freq)
pi1_down_df <- as.data.frame(as.table(pi1_down)) %>% rename(Limma = Var1, MAST = Var2, pi1 = Freq)
pi1_nm_up_df <- as.data.frame(as.table(pi1_nm_up)) %>% rename(Limma = Var1, MAST = Var2, pi1 = Freq)
pi1_nm_down_df <- as.data.frame(as.table(pi1_nm_down)) %>% rename(Limma = Var1, MAST = Var2, pi1 = Freq)
# Plot heatmap function
plot_pi1_heatmap <- function(df, title){
  ggplot(df, aes(x = MAST, y = Limma, fill = pi1)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(pi1, 2)), color = "black", size = 4)+
    scale_fill_gradient(low = "white", high = "maroon", limits = c(0,1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16)) +
    labs(title = title, x = "Leading (MAST)", y = "Matched (Limma)", fill = expression(pi[1]))
}

npglut_up <- plot_pi1_heatmap(pi1_up_df, "npglut Upregulated Genes")
npglut_down <- plot_pi1_heatmap(pi1_down_df, "npglut Downregulated Genes")
nmglut_up <- plot_pi1_heatmap(pi1_nm_up_df, "nmglut Upregulated Genes")
nmglut_down <-plot_pi1_heatmap(pi1_nm_down_df, "nmglut Downregulated Genes")


ggsave(npglut_up, filename = "Pi1_npglutup.png", width = 4.7, height = 4, bg = "white")
ggsave(npglut_up, filename ="Pi1_npglutup.pdf", width = 4.7, height = 4, bg = "white")
ggsave(npglut_down, filename ="Pi1_npglutdown.png", width = 4.7, height = 4, bg = "white")
ggsave(npglut_down,filename = "Pi1_npglutdown.pdf", width = 4.7, height = 4, bg = "white")
ggsave(nmglut_up, filename ="Pi1_nmglutup.png", width = 4.7, height = 4, bg = "white")
ggsave(nmglut_up,filename = "Pi1_nmglutup.pdf", width = 4.7, height = 4, bg = "white")
ggsave(nmglut_down,filename = "Pi1_nmglutdown.png", width = 4.7, height = 4, bg = "white")
ggsave(nmglut_down, filename ="Pi1_nmglutdown.pdf", width = 4.7, height = 4, bg = "white")
pi1_np_up_df <- pi1_up_df
pi1_np_down_df <- pi1_down_df
save(pi1_up_df, pi1_down_df, pi1_nm_up_df, pi1_nm_down_df, 
     file = "all_pi1_values.RData")
