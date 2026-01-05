library(stringr)
library(ggplot2)
library(RColorBrewer)
setwd("~/GREAT")
df_raw_data <-
  read.table(file = "ASoC_25_ways_sum_dir/intersect_table_EnhvsProm.tsv",
             header = F, sep = "\t",
             quote = "")
colnames(df_raw_data) <-
  c("Sample_SNP", "SNP_coverage", "Interval_region", "Interval_coverage", "Overlap_coverage")
remove <- c("sharedByAllCellTypes_0h", "sharedByAllCellTypes_0h",
            "sharedByAllCellTypes_1h", "sharedByAllCellTypes_1h",
            "sharedByAllCellTypes_6h" ,"sharedByAllCellTypes_6h",
            "atLeastInOneCellType_0h", "atLeastInOneCellType_1h",
            "atLeastInOneCellType_6h",
            "GABA_at_least_one", "nmglut_at_least_one", "npglut_at_least_one")
df_raw_data <- df_raw_data[! df_raw_data$Sample_SNP %in% remove, ]
df_2_plot <- df_raw_data
df_2_plot$Percentage <-
  df_2_plot$Overlap_coverage / df_2_plot$SNP_coverage

df_2_plot$Cell_type <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 1]
df_2_plot$Cell_time <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 2]

df_2_plot$Region <-
  str_split(df_2_plot$Interval_region,
            pattern = "_",
            n = 2,
            simplify = T)[, 2]
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12F",
                  replacement = "Me_F")
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12O",
                  replacement = "Me_O")

df_2_plot$Region[df_2_plot$Region %in% "zEnh_all"] <- "Enhancers"
df_2_plot$Region[df_2_plot$Region %in% "zProm"] <- "Promoters"


unique(df_2_plot$Cell_type)
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type == "atLeastInOneCellType"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type == "sharedByAllCellTypes"), ]
# df_2_plot <-
#   df_2_plot[(df_2_plot$Cell_type == "npglut"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_time == "at_least_one"), ]

# df_2_plot

plot <- ggplot(df_2_plot,
               aes(x = Cell_time,
                   y = Percentage)) +
  geom_bar(aes(fill = Region,
               colour = Cell_type),
           stat = "identity",
           position = "dodge",
           color = NA,
           # position = position_dodge(width = 1),
           # width = 1,
           linewidth = 0) +
  geom_text(aes(label = paste(df_2_plot$Overlap_coverage,
                              df_2_plot$SNP_coverage,
                              sep = ' / ')),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0.5,
            position = position_dodge2(width = 0.9),
            colour = "black") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1"))  +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Set2"),
                      guide = "none") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1.0),
                     labels = scales::percent) +
  labs(x = "") +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 315,
                                          vjust = 0.5,
                                          hjust = 0,
                                          colour = "black")) +
  facet_grid(. ~ Cell_type)

ggsave(plot = plot,
      filename = "FigS17A.pdf", width = 7, 
       height = 5, dpi = 300)
ggsave(plot = plot,
       filename = "FigS17A.png", width = 7, 
       height = 5, dpi = 300, bg = "white")
