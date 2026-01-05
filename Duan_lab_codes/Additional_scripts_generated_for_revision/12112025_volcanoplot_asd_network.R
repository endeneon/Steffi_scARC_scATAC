library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggrastr)

load("./Use_log_normalization_and_harmony/Using_Cellline_As_covar/All_combined_comparison.RData")
asd_targets <- read_csv("PT_GRN_Fu_et_al_ASD_Network.csv")
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
all_results <- all_results[!grepl("*CPT", names(all_results))]
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
all_results <- lapply(all_results, function(df){
  colnames(df)[colnames(df) == "regulation"] <- "significance"
  df$significance <- dplyr::recode(df$significance, 
                                   "ns" = "nonsig",
                                   "down" = "neg",
                                   "up" = "pos",
  )
  return(df)
})
save_path <- "./Volcano_plots/"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
for (name in names(all_results)){
  res <- all_results[[name]]
  res$gene <- rownames(res)
  res$neg_log_p_val <- -log10(res$p_val)
  res$significant <- (abs(res$avg_log2FC) > 0 & res$p_val_adj < 0.05)
  parts <- strsplit(name, "_")[[1]]
  celltype <- parts[1]
  gene <- parts[2]
  time <- parts[3]
  direct_name <- paste0(gene, "_", celltype)
  targets <- tf_list[[direct_name]]
  sig_counts <- table(res$significance)
  volcano_plot <- ggplot(data = as.data.frame(res),
                         aes(x = avg_log2FC, 
                             y = neg_log_p_val, 
                             color = significance)) + 
    rasterise(geom_point(size = 1)) +
    scale_color_manual(values = c("pos" ="red3", 
                                  "nonsig" ="grey50", 
                                  "neg" = "blue"),
                       labels = c(
                         "pos" = paste0("pos (", sig_counts["pos"], ")"),
                         "nonsig" = paste0("nonsig (", sig_counts["nonsig"], ")"),
                         "neg" = paste0("neg (", sig_counts["neg"], ")")
                       )) +
    theme_minimal() +
    geom_text_repel(data = subset(res, significant & gene %in% targets),
                    aes(label = gene),
                    box.padding = unit(0.5, 'lines'),
                    point.padding = unit(0.5, 'lines'),
                    segment.size = 0.3,
                    # size = 3, 
                    #min.segment.length = 0.5,
                    max.overlaps = Inf,
                    force_pull = 10,
                    show.legend = F,
                    nudge_y = 20, 
                    nudge_x = ifelse(subset(res, significant & gene %in% targets)$avg_log2FC > 0, 1.5, -1.5))+ 
    xlim(c(-5, 5)) +
    theme(
      plot.title = element_text(hjust = 0.5, size =),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      # axis.line = element_line(size =1),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    ggtitle(paste0(gene, "_KO", " vs. WT", " (", celltype, ", ", time, ")")) +
    labs(
      x = "logFC",
      y = "-log10(P-value)")
  volcano_plot
  ggsave(paste0(save_path, name,
                "_volcano_plot.png"),
         plot = volcano_plot,
         width = 7,
         height = 5,
         dpi = 1200,
         bg = "white")
  ggsave(paste0(save_path, name,
                "_volcano_plot.pdf"),
         plot = volcano_plot,
         width = 7,
         height = 5,
         dpi = 1200)
}
