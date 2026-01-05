# Skip the processing and directly load the results
library(locuszoomr)
library(ctwas)
library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(ggrepel)
library(patchwork)
library(data.table)

geneloc <- fread("/path/to/geneloc.txt.gz")
mapping_table <- data.frame(molecular_id=finemap_res$molecular_id, gene_type="protein_coding",
                            gene_name=sapply(strsplit(finemap_res$molecular_id,"\\."),function(x){x[2]}))
mapping1 <- merge(mapping_table, geneloc, by.x="gene_name",by.y="geneid")
colnames(mapping1)[4:6] <- c("chrom","start","end")
snp_map <- readRDS("/path/to/snp_map.RDS")
weights <- readRDS("/path/to/SCZ.preprocessed.weights.RDS")
finemap_res_multi <- anno_finemap_res(finemap_res,
                                      snp_map = snp_map,
                                      mapping_table = mapping1,
                                      add_gene_annot = TRUE,
                                      map_by = "molecular_id",
                                      drop_unmapped = TRUE,
                                      add_position = TRUE,
                                      use_gene_pos = "mid")

### Plot for ASOC_rs72986630
data(SLE_gwas_sub)
temp <- list(chr=19, pos=11738921)
loc <- locus(xrange=c(temp$pos[1]-2.5e5,temp$pos[1]+2.5e5), seqname=temp$chr, 
             data=SLE_gwas_sub, ens_db = "EnsDb.Hsapiens.v86")

datAll <- finemap_res_multi[finemap_res_multi$chrom==loc$seqname & finemap_res_multi$pos>loc$xrange[1] & finemap_res_multi$pos<loc$xrange[2],
                            c("type","molecular_id","context","susie_pip","chrom","pos")]
datAll$timepoint <- sapply(strsplit(datAll$context,"_"), function(x){x[1]})
datAll$celltype <- sapply(strsplit(datAll$context,"_"), function(x){x[2]})
datAll[is.na(datAll)] <- "SNP" ## Fill in time point and cell types for SNP
datAll$pos <- datAll$pos/1e6

datAll$label <- ""
sel <- (datAll$susie_pip>0.1) & (datAll$type!="SNP")
datAll$label[sel] <- sapply(strsplit(datAll$molecular_id[sel], "\\."), function(x){paste0(x[2],",",x[1])})
datAll$timepoint[datAll$timepoint!="SNP"] <- paste0(datAll$timepoint[datAll$timepoint!="SNP"],"_gene")
datAll$celltype[datAll$celltype!="SNP"] <- paste0(datAll$celltype[datAll$celltype!="SNP"],"_gene")
colnames(datAll)[7:8] <- c("color","shape")

alpha_vec <- rep(1,length(unique(datAll$type)))
names(alpha_vec) <- unique(datAll$type)
alpha_vec["SNP"] <- 0.5

size_vec <- rep(2,length(unique(datAll$type)))
names(size_vec) <- unique(datAll$type)
size_vec["SNP"] <- 1

shape_vec <- c(17,15,18,16)
names(shape_vec) <- c("GABA_gene","nmglut_gene","npglut_gene","SNP")

color_vec <- c("orange","red","purple","blue")
names(color_vec) <- c("0hr_gene","1hr_gene","6hr_gene","SNP")

p1 <- ggplot(datAll, aes(x=pos, y=susie_pip, shape=shape, color=color, size=type, alpha=type, label=label)) + 
  geom_point() + #geom_vline(xintercept = temp$pos[1]/1e6, linetype="dashed", color = "black", size=1, alpha=0.7) + 
  #  facet_grid(celltype~., scales = "free_y") + 
  #  annotate("text", temp$pos[1]/1e6, 1, hjust = -.25, 
  #           label = paste(near.gene,"eQTL")) +
  scale_alpha_manual(values=alpha_vec, guide="none") +
  scale_size_manual(values=size_vec, guide="none") +
  scale_shape_manual(values=shape_vec, guide="none") + 
  scale_color_manual(values=color_vec, name=NULL) + 
  ylab("cTWAS PIP") + xlim(loc$xrange/1e6) + ylim(0,1) + 
  xlab("") + theme_bw() +
  geom_label_repel(data=datAll[sel,], size=2) + 
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))

datP <- finemap_res_multi[finemap_res_multi$chrom==loc$seqname & finemap_res_multi$pos>loc$xrange[1] & finemap_res_multi$pos<loc$xrange[2],
                          c("type","molecular_id","context","z","chrom","pos")]
datP$logP <- -log10(2*pnorm(-abs(datP$z)))
datP$timepoint <- sapply(strsplit(datP$context,"_"), function(x){x[1]})
datP$celltype <- sapply(strsplit(datP$context,"_"), function(x){x[2]})
datP[is.na(datP)] <- "SNP" ## Fill in the time point and cell types for SNP
datP$pos <- datP$pos/1e6
datP$label <- ""
sel <- (datP$type!="SNP") & (datP$logP>6)
datP$label[sel] <- sapply(strsplit(datP$molecular_id[sel], "\\."), function(x){paste0(x[2],",",x[1])})
datP$timepoint[datP$timepoint!="SNP"] <- paste0(datP$timepoint[datP$timepoint!="SNP"],"_gene")
datP$celltype[datP$celltype!="SNP"] <- paste0(datP$celltype[datP$celltype!="SNP"],"_gene")
colnames(datP)[8:9] <- c("color","shape")

ymax <- max(datP$logP)
p2 <- ggplot(datP, aes(x=pos, y=logP, shape=shape, color=color, size=type, alpha=type, label=label)) + 
  geom_point() + #geom_vline(xintercept = temp$pos[1]/1e6, linetype="dashed", color = "black", size=1, alpha=0.7) + 
  #  facet_grid(celltype~., scales = "free_y") + 
  #  annotate("text", temp$pos[1]/1e6, 1, hjust = -.25, 
  #           label = paste(near.gene,"eQTL")) +
  scale_alpha_manual(values=alpha_vec, guide="none") +
  scale_size_manual(values=size_vec, guide="none") +
  scale_shape_manual(values=shape_vec, name=NULL) + 
  scale_color_manual(values=color_vec, guide="none") + 
  ylab("GWAS/TWAS -log10 P") + xlim(loc$xrange/1e6) + ylim(0,ymax) + 
  xlab("") + theme_bw() +
  geom_label_repel(data=datP[sel,], size=2) + 
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))

p3 <- gg_genetracks(loc, filter_gene_biotype = "protein_coding")

pw <- (p2 + theme(legend.position = "right",legend.box.spacing = unit(0, "pt"),
                  axis.title.x=element_blank(), axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(), plot.margin=margin(0,6,0,6))) /
  (p1 + theme(legend.position = "right",legend.box.spacing = unit(1, "pt"),
              axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks.x=element_blank(), plot.margin=margin(0,6,0,6))) /
  (p3 + theme(plot.margin=margin(0,6,0,6))) +
  plot_layout(ncol = 1) & theme(axis.title.y = element_text(face = "bold"))
pw
ggsave("ASoC_rs7298663_locus_plot.pdf", pw, width = 8, height = 4.5)