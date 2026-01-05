library(ggplot2)

enrich_files <- list.files(path = "./rdata/",
                           pattern = "^enrichment_.*\\.rdata$", full.names = T)

enrich_all <- list()

for (f in enrich_files) {
  load(f)  # loads object called 'enrichment'
  
  enrich_all[[f]] <- enrichment
}

enrichment_all <- do.call(rbind, enrich_all)
se <- (enrichment_all$high - enrichment_all$low) / (2*1.96)
z <- enrichment_all$estimate/se
p <- exp(-0.717*2 - 0.416*z^2)

enrichment_all$p <- p
enrichment_all$lp <- -log10(p)
class_snp <- levels(as.factor(enrichment_all$snp))
order <- unique(enrichment_all$snp)
enrichment_all$snp <- factor(enrichment_all$snp, levels = rev(order))

rename_traits <-  c("SCZ2021" = "SCZ",
                    "Bipolar2024" = "Bipolar",
                    "MDD2025" = "Depression",
                    "neuroticism2018" = "Neuroticism", 
                    "ASD2019"="Autism Spectrum Disorder", 
                    "ADHD"="ADHD",
                    "Intelligence2018"="Intelligence", 
                    "ALcoholic2018"="Alcohol Dependence" , 
                    "PTSD"="PTSD", 
                    "Alz2022"="Alzheimers's Disease", 
                    "Parkinson2019"= "Parkinson's",
                    "Height2014"="Height", 
                    "CD"="Crohn's Disease", 
                    "UC" ="Ulcerative Colitis",
                    "T2D2018"="Type 2 Diabetes")

enrichment_all$trait <- ifelse(enrichment_all$trait %in% names(rename_traits),
                               rename_traits[enrichment_all$trait],
                               enrichment_all$trait) 
enrichment_all$trait <- factor(enrichment_all$trait, 
                               levels = c("SCZ", "Bipolar","Depression", "Neuroticism", 
                                          "Autism Spectrum Disorder", "ADHD",
                                          "Intelligence", "Alcohol Dependence" , 
                                          "PTSD", "Alzheimers's Disease", 
                                          "Parkinson's","Height", "Crohn's Disease", 
                                          "Ulcerative Colitis",
                                          "Type 2 Diabetes"))
limit <- max(abs(enrichment_all$estimate)) * c(-1, 1)
enrichment_all$size_capped <- pmin(enrichment_all$lp, 9)

torus_plot <- ggplot(enrichment_all, aes(x = trait, y = snp, 
                                         size = size_capped, fill = estimate)) +
  geom_point(
    shape = 21,
    color = "black",     # thin black outline
    stroke = 0.3         # very thin outline like your image
  ) +
  scale_radius(limits = c(0, 9),
               breaks = c(0, 3, 6, 9),
               range = c(1, 8))+
  scale_fill_gradient2(
    low  = "blue",
    mid = "white",
    high = "#a50f15",
    limits = c(-5, 5),
    breaks = seq(-5, 5, 2.5)  
  ) +
  labs(fill = "log2OR", size = "-log10P", x = "", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.2),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1)
  ) 
print(torus_plot)

ggsave(torus_plot,
       filename = "Torus_dotplot_output.jpg", 
       width = 7, height = 5, dpi = 600)

ggsave(torus_plot,
       filename = "Torus_dotplot_output.pdf", 
       width = 7, height = 5, dpi = 600)
