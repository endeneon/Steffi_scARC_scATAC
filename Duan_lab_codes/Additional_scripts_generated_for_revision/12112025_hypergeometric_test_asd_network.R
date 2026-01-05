library(dplyr)
library(stats)

load("./Use_log_normalization_and_harmony/Using_Cellline_As_covar/All_combined_comparison.RData")
asd_targets <- read_csv("PT_GRN_Fu_et_al_ASD_Netowrk.csv")
colnames(asd_targets)

unique_tf <- unique(asd_targets$tf)
unique_tf <- setdiff(unique_tf, "MLXIP")
tf_list <- list()
for(tf_name in unique_tf){
  df <- asd_targets[asd_targets$tf == tf_name, ]
  celltypes <- setdiff(unique(df$ObservedIn), "Both")
  
  for(types in celltypes){
    df_subset <- df[df$ObservedIn == types | df$ObservedIn == "Both", ]
    types_clean <- sub(" only$", "", types)
    list_name <- paste0(tf_name, "_", types_clean)
    tf_list[[list_name]] <- df_subset$gene
  }
}

name_map <- c( "GABA_0_MEF" = "GABA_MEF2C_0hr",
               "GABA_1_MEF" ="GABA_MEF2C_1hr",
               "GABA_6_MEF"="GABA_MEF2C_6hr",
               "GABA_0_ROR" ="GABA_RORB_0hr",
               "GABA_1_ROR" ="GABA_RORB_1hr",
               "GABA_6_ROR" ="GABA_RORB_6hr",
               "GABA_0_TCF"="GABA_TCF4_0hr",
               "GABA_1_TCF" ="GABA_TCF4_1hr",
               "GABA_6_TCF"="GABA_TCF4_6hr",
               "glut_0_MEF"="Glut_MEF2C_0hr",
               "glut_1_MEF" ="Glut_MEF2C_1hr",
               "glut_6_MEF"="Glut_MEF2C_6hr",
               "glut_0_ROR" ="Glut_RORB_0hr",
               "glut_1_ROR"="Glut_RORB_1hr",
               "glut_6_ROR" ="Glut_RORB_6hr",
               "glut_0_TCF" ="Glut_TCF4_0hr",
               "glut_1_TCF" ="Glut_TCF4_1hr",
               "glut_6_TCF"="Glut_TCF4_6hr"
)
names(all_results) <- name_map[names(all_results)]
results <- lapply(names(all_results), function(deg_name){
  
  parts <- strsplit(deg_name, "_")[[1]]
  celltype <- parts[1]
  gene <- parts[2]
  
  direct_name <- paste0(gene, "_", celltype)
  all_genes <- rownames(all_results[[deg_name]])
  sig_genes <- rownames(all_results[[deg_name]][all_results[[deg_name]]$p_val_adj < 0.05, ,drop= FALSE])
  
  direct_genes <- tf_list[[direct_name]]
  
  direct_in_background <- intersect(direct_genes, all_genes)
  k <- sum(sig_genes %in% direct_in_background)
  
  N <- length(all_genes) ##total tested genes in that run
  M <- length(direct_in_background) # direct in universe
  K <- length(sig_genes) #DEGs
  pval <- if (M > 0 && K > 0){
    phyper(q = k-1, 
           m = M, 
           n = N-M, 
           k = K, lower.tail = FALSE)
  } else {
    NA
  }
  fold_enrichment <- if (M > 0 && K > 0){
    k / (K * (M / N))
  } else {
    NA
  }
  
  overlap_genes <- intersect(sig_genes, direct_in_background)
  data.frame(
    deg_name = deg_name,
    direct_name = direct_name,
    overlap = k,
    overlap_genes = paste(overlap_genes, collapse = ","),
    deg_total = K,
    direct_total = M,
    universe = N,
    fold_enrichment = fold_enrichment,
    pval = pval
  )
  
})


results_df <- bind_rows(results)
results_df$padj <- p.adjust(results_df$pval, method = 'BH')
save(results_df, file ="Overlap_direct_targets_test_for_random.RData")


write.csv(results_df, file ="Overlap_direct_targets_test_for_random.csv")

##To make heatmap
library(ggplot2)
library(tidyr)
library(ggplot2)
library(stringr)
heatmap_df <- results_df %>%
  mutate(number = str_extract(deg_name, "\\d+(?=hr)")) %>%
  mutate(sig = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ ""
  )) %>%
  select(number, direct_name, fold_enrichment, sig)

heatmap_df$number <- factor(
  heatmap_df$number,
  levels = sort(as.numeric(unique(heatmap_df$number)))
)
heatmap_df$number <- paste0(heatmap_df$number, "hr")


desired_order <- c(
  "RORB_Glut",
  "RORB_GABA",
  "MEF2C_Glut",
  "MEF2C_GABA",
  "TCF4_Glut",
  "TCF4_GABA"
)

heatmap_df$direct_name <- factor(heatmap_df$direct_name, levels = desired_order)

plot <- ggplot(heatmap_df, aes(x = direct_name, y = number, fill = fold_enrichment)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sig), size = 5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Fold Enrichment"
  ) +
  labs(
    x = "ASD direct targets",
    y = "DEG comparison",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
plot
ggsave(plot, filename = "enrichment_heatmap.pdf", dpi = 600,
       height = 3, width =6)
ggsave(plot, filename = "enrichment_heatmap.png", dpi = 600,
       height = 3, width =6, bg = "white")
