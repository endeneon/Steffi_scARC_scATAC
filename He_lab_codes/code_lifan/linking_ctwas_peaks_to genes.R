## cTWAS version 0.5.35
library(ctwas)

finemap_res <- readRDS("/path/to/MDD_2025.finemap_regions_res.RDS")
finemap_res <- finemap_res$susie_alpha_res
finemap_res$gene_name <- sapply(strsplit(finemap_res$molecular_id, "\\."), function(x){x[2]})
# finemap_res1 <- anno_susie_alpha_res(finemap_res,
#                                      mapping_table = mapping_table,
#                                      map_by = "molecular_id",
#                                      drop_unmapped = TRUE)

combined_pip <- combine_gene_pips(finemap_res,
                                  group_by = "gene_name",
                                  by = "context",
                                  method = "combine_cs",
                                  filter_cs = F,
                                  include_cs_id = FALSE)


weight.store <- "ctwas_SCZ/weights/Joint_caQTL_ASoC_SCZ_09112024/"
weightfs <- list.files(weight.store)[1:9]
ws <- list()
for(w in weightfs) {
  n <- strsplit(w,"\\.")[[1]][1]
  mydb <- dbConnect(RSQLite::SQLite(), paste0(weight.store,w))
  ws[[n]] <- dbGetQuery(mydb, 'SELECT * FROM weights')
}
res4$condition <- colnames(res4)[3:11][apply(res4[,3:11],1,which.max)]
res4$condition <- substr(res4$condition,1,nchar(res4$condition)-4)
res4$snp <- "0"
for(i in 1:nrow(res4)) {
  if(!(res4$gene_name[i] %in% ws[[res4$condition[i]]]$gene)) next
  res4$snp[i] <- ws[[res4$condition[i]]][ws[[res4$condition[i]]]$gene== res4$gene_name[i],"rsid"]
}


### Mapping genes with eQTL
eqtl.store <- "/path/to/eQTL/store"
eqtl.files <- list.files(eqtl.store,pattern="eqtl.rds")
eqtls <- list()
for(f in eqtl.files[-10]) {
  n <- substr(f, 1, nchar(f)-47)
  eqtls[[n]] <- readRDS(paste0(eqtl.store,f))$cis$eqtls
}
eqtl1 <- lapply(eqtls,function(x){x[x$snps %in% res4$snp,]})
#saveRDS(eqtl1,"/project/xinhe/lifan/neuron_stim/ctwas_SCZ/eqtl_overlap_ctwas_peaksPIP5.rds")

eqtl2 <- lapply(eqtl1, function(x){x[x$FDR<0.2,]})
eqtl3 <- do.call(rbind, eqtl2)
res4$eqtl.gene <- "NA"
for(s in unique(eqtl3$snps)) {
  res4[res4$snp==s,"eqtl.gene"] <- paste(unique(eqtl3[eqtl3$snps==s,"gene"]),collapse=",")
}

### Mapping genes with co-activation

p2g <- readRDS("/path/to/peak_target_gene_results.rds")
p2g1 <- lapply(p2g, function(x){x[x$fdr<0.2,]})
p2g2 <- p2g1
for(n in names(p2g1)) {
  p2g2[[n]]$celltype <- n
}
#p2g2 <- lapply(p2g1, function(x){x[x$Peak %in% res4$range,]})
p2g3 <- do.call(rbind, p2g2)

### Construct peak ranges
res4$range <- "NA"
cpeaks <- strsplit(res4$gene_name[startsWith(res4$gene_name,"cPeak")],"_")
cpeaks <- sapply(cpeaks,function(x){paste0(x[2],":",x[3],"-",x[4])})
res4$range[startsWith(res4$gene_name,"cPeak")] <- cpeaks
res4$p2g.gene <- "NA"
for(p in unique(p2g3$Peak)) {
  res4[res4$range==p,"p2g.gene"] <- paste(unique(p2g3[p2g3$Peak==p,"Gene"]),collapse=",")
}

asoc.store <- "/path/to/ASoC_results_folder"
asocfs <- list.files(asoc.store)
asocs <- list()
for(f in asocfs) {
  n <- substr(f,1,nchar(f)-18)
  asocs[[n]] <- read.table(paste0(asoc.store,f),header=T)
}

asoc.snp <- sapply(strsplit(res4$gene_name[startsWith(res4$gene_name,"ASOC")],"_"),function(x){x[2]})
asoc.cond <- res4$condition[startsWith(res4$gene_name,"ASOC")]
asoc.range <- rep("",length(asoc.snp))
for(i in seq_along(asoc.snp)) {
  temp <- asocs[[asoc.cond[i]]][asocs[[asoc.cond[i]]]$ID==asoc.snp[i],]
  asoc.range[i] <- paste0(temp$CHROM,":",temp$POS-250, "-", temp$POS+250)
}
res4$range[startsWith(res4$gene_name,"ASOC")] <- asoc.range
asoc.range1 <- GRanges(asoc.range)
#p2g1.range <- lapply(p2g1, function(x){GRanges(x$Peak)})
p2g2 <- lapply(p2g1, function(x){x[findOverlaps(asoc.range1,GRanges(x$Peak))@to,]})
p2g3 <- do.call(rbind,p2g2)
match <- findOverlaps(asoc.range1,GRanges(p2g3$Peak))
for(i in unique(match@from)) {
  res4[res4$snp==asoc.snp[i],"p2g.gene"] <- paste(unique(p2g3[match[match@from==i]@to,"Gene"]),collapse=",")
  res4[res4$snp==asoc.snp[i],"p2g.peak"] <- paste(unique(p2g3[match[match@from==i]@to,"Gene"]),collapse=",")
}
saveRDS(res4[order(res4$PIP_sum,decreasing=T),], "SCZ_cTWAS_peak_target_gene.rds")
res4 <- readRDS("SCZ_cTWAS_peak_target_gene.rds")



res4.p2g <- res4[res4$p2g.gene!="NA",c("genename","snp","range","PIP_sum","condition","p2g.gene")]
N.p2g <- sapply(strsplit(res4.p2g$p2g.gene,","),length)
res5.p2g <- res4.p2g[rep(rownames(res4.p2g),times=N.p2g),-5]
res5.p2g$p2g.gene <- unlist(strsplit(res4.p2g$p2g.gene,","))
lapply(p2g1, function(x){x[(x$Gene==res5.p2g$p2g.gene[1]) & (x$Peak==res5.p2g$genename[1]),]})


### Compute distance to the nearest genes
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

res4 <- readRDS("ctwas_SCZ/SCZ_cTWAS_peak_target_gene_withABC.rds")
pcg <- read.table("symbol_entrez_autosome_protein_coding.txt")
prom <- promoters(TxDb.Hsapiens.UCSC.hg38.knownGene,upstream = 0,downstream = 0)
g = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
g1 <- g[g$gene_id %in% pcg$entrezgene_id]
tss <- resize(g1, width=1, fix='start')
dist <- distanceToNearest(GRanges(res4$range), tss)

#gene.match <- findOverlaps(GRanges(res4$range), tss)
entrez2symbol <- select(org.Hs.eg.db, key=tss$gene_id, columns="SYMBOL", keytype="ENTREZID")
res4$tss.gene <- "NA"
res4$tss.dist <- Inf
for(i in unique(dist@from)) {
  geneset <- tss[dist[dist@from==i]@to]$gene_id
  res4$tss.gene[i] <- paste(entrez2symbol[entrez2symbol$ENTREZID %in% geneset,"SYMBOL"],collapse=",")
  res4$tss.dist[i] <- dist@elementMetadata$distance[dist@from==i]
}
saveRDS(res4[order(res4$PIP_sum,decreasing=T),], "SCZ_cTWAS_peak_target_gene_withABC_filterTSSgenes.rds")

### Linking to target gene with ABC score
library(data.table)
abc <- fread("ctwas_SCZ/All_ABC_EnhancerPredictions_threshold0.021_self_promoter.tsv")
abc.range <- makeGRangesFromDataFrame(abc,keep.extra.columns = T)
abc.overlap <- findOverlaps(GRanges(res4$range), abc.range)

saveRDS(res4[order(res4$PIP_sum,decreasing=T),], "SCZ_cTWAS_peak_target_gene_withABC.rds")
res4$ABC.gene <- "NA"
res4$ABC.max <- "NA"
for(i in unique(abc.overlap@from)) {
  abc.sub <- abc.range[abc.overlap[abc.overlap@from==i]@to,]
  geneset <- unique(abc.sub$TargetGene)
  res4$ABC.gene[i] <- paste(geneset,collapse=",")
  res4$ABC.max[i] <- abc.sub$TargetGene[which.max(abc.sub$ABC.Score)]
}

saveRDS(res4[order(res4$PIP_sum,decreasing=T),], "SCZ_cTWAS_peak_target_gene_withABCmax.rds")
write.table(res4[order(res4$combined_pip,decreasing=T),],"ctwas_SCZ/new_SCZ_cTWAS_peak_target_gene_withABCmax.tsv",quote=F,row.names=F,sep="\t")
