library(data.table)
library(rtracklayer)

setwd("path/to/our/neuron_stim/eqtls")

eqtlfs <- list.files()
eqtls <- list()
for(f in eqtlfs) {
  n <- substr(f, 1, nchar(f)-37)
  #eqtls[[n]] <- readRDS(f)$cis$eqtls
  eqtls[[n]] <- fread(f)
}

## Mapping SNPs from Doparminergic to rsid

snpinfo <- fread("path/to/snpinfo_100lines.txt") # effect alleles and genomic positions
snpinfo$variant_id <- paste(snpinfo$chr, snpinfo$pos, snpinfo$alt, sep="_")

ref_pairs <- readRDS("../rb_est/ref_pairs.rds") # Reference gene-SNP pairs for computing rb

## core function to Estimate rb between one set of eQTL to another
rb_est <- function(eqtl_ns, eqtl_jb, ref_pairs) {
  null_ns <- eqtl_ns[gene %in% ref_pairs$gene & `p-value`>0.01,1:3]
  null_jb <- eqtl_jb[hgnc_symbol %in% ref_pairs$gene & `pval`>0.01,]
  null.sumstats <- merge(null_ns, null_jb, by.x=c("SNP","gene"), by.y=c("rsid","hgnc_symbol"), all=F)
  dtCor <- null.sumstats[, .(mCor = cor(beta.x,beta.y)), by=gene]
  re <- mean(dtCor$mCor, na.rm=T)
  
  pair1 <- merge(ref_pairs[,1:2], eqtl_ns, by=c("gene","SNP"), all=F)
  pair1$se <- pair1$beta / pair1[,"t-stat"]
  pair2 <- merge(ref_pairs[,1:2], eqtl_jb, by.x=c("gene","SNP"), by.y=c("hgnc_symbol","rsid"), all=F)
  #pair2$se <- pair2$beta / pair2[,"t-stat"]
  
  var.e1 <- mean(pair1$se^2)
  var.e2 <- mean(pair2$beta_se^2)
  
  var.b1 <- var(pair1$beta)
  var.b2 <- var(pair2$beta)
  
  pairs <- merge(pair1, pair2, by=c("gene","SNP"),all=F)
  covar.bhat <- cov(pairs$beta.x, pairs$beta.y)
  
  rb <- (covar.bhat - re*sqrt(var.e1*var.e2))/sqrt(var.b1-var.e1)/sqrt(var.b2-var.e2)
  rb
}

## Load doparminergic dataset
## lift over SNP to HG38
## map genes to gene symbols
## Overlap with SNP info and obtains rsid

ch <- import.chain("../../jerber_dopamin/hg19ToHg38.over.chain")

da.store <- "../../jerber_dopamin/eqtl_summary_stats_renamed/"
dafs <- list.files(da.store, pattern="gz$")
dans <- substr(dafs, 1, nchar(dafs)-30)

ebs <- matrix(0,nrow=9, ncol=length(dafs))
rownames(ebs) <- names(eqtls)
colnames(ebs) <- dans

for (i in 1:length(dafs)) {
  temp <- fread(paste0(da.store, dafs[i]))
  temp1 <- makeGRangesFromDataFrame(temp, seqnames.field = "snp_chromosome", start.field = "snp_position", 
                                    end.field="snp_position", keep.extra.columns=T)
  seqlevelsStyle(temp1) <- "UCSC"
  temp2 <- unlist(liftOver(temp1, ch))
  temp3 <- data.table(feature = temp2$feature_id, 
                      alt=temp2$assessed_allele, 
                      beta=temp2$beta,
                      beta_se=temp2$beta_se,
                      pval=temp2$p_value,
                      snp_chr=as.character(seqnames(temp2)),
                      snp_pos=start(temp2))
  temp3$variant_id <- paste(temp3$snp_chr, temp3$snp_pos, temp3$alt, sep="_")
  genes <- readRDS("../../geneMapping.rds")
  temp4 <- merge(temp3, genes, by.x="feature", by.y="ensembl_gene_id", all=F)
  
  temp5 <- merge(temp4, snpinfo, by="variant_id", all=F)
  
  for(j in 1:length(eqtls)) {
    cat("running:", names(eqtls)[j], "&", dans[i], "\n")
    ebs[names(eqtls)[j], dans[i]] <- rb_est(eqtls[[j]], temp5, ref_pairs)
  }
}

ebs ## The matrix of rb estimate between every context of our eQTLs and every context of those in DA neurons (Jerber, et al.)

