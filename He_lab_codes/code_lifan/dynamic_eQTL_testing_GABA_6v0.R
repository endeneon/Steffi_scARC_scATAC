library(data.table)

snps = fread('/path/to/genotype_information.txt')

covars = read.table("covariates.txt")
exPC_GABA = read.table("GABA_residual_expression_PC40.txt")

covar_GABA <- merge(covars, exPC_GABA[,1:20], all=F, by=0)
covar_GABA[['offset']] = 1
covar_GABA[["indiv"]] = sapply(strsplit(covar_GABA$Row.names,"\\|"), function(x){x[1]})
covar_GABA[["timepoint"]] = sapply(strsplit(covar_GABA$Row.names,"\\|"), function(x){x[2]})
expr_GABA = read.table("mateqtl_input_100lines/psbulk_typenorm_GABA.txt",sep="\t")

#pstop = read.table("limix_qtl/GABA_candidate_pairs.txt.gz",sep=" ")
rest = readRDS("/path/to/GABA_0hr_main_effect_eqtls.rds")
stim = readRDS("/path/to/GABA_6hr_main_effect_eqtls.rds")
pstop = merge(rest[rest$FDR<0.2,c("gene","snps")], stim[stim$FDR<0.2,c("gene","snps")], by=c("gene","snps"), all=T)
pstop$pval = 1
pstop$beta = 0
pstop[,"t-stat"] = 0
cov1 <- covar_GABA[covar_GABA$timepoint %in% paste0(c(0,6),"hr_GABA"),]

for (i in 1:nrow(pstop)) {
  geno = as.numeric(snps[snps$V1==pstop$snps[i], cov1$indiv, with=F])
  y = as.numeric(expr_GABA[pstop$gene[i],gsub("\\|","\\.",cov1$Row.names)])
  dat = cbind(y,geno,cov1)
  #m0 <- as.formula(paste0("y~",paste(colnames(dat)[-c(1,3,36)],collapse = "+")))
  m1 <- as.formula(paste0("y~",paste(colnames(dat)[-c(1,3,36)],collapse = "+"),"+geno*timepoint"))
  #lm0 <- lm(m0, dat)
  lm1 <- lm(m1, dat)
  # We don't need to compare alternative model with the null model
  # Because the p value is equivalent to t-stat P values in the fixed-effct model
  pstop[i,c("beta","t-stat","pval")]<- summary(lm1)$coefficients[35,c(1,3,4)]
}

pstop1 <- pstop[,c("snps","gene","beta","t-stat","pval")]
colnames(pstop1) <- c("SNP","gene","beta","t-stat","p-value")
fwrite(pstop1, "dynamic_QTL_fixed_pairwise/GABA_6v0_fdr20.txt.gz", quote=F, sep="\t")







