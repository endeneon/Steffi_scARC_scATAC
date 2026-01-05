library(ctwas)
library(ggplot2)
library(tools)
library(logging)
library(doParallel)
library(data.table)
library(optparse)

### Analysis script for running multigroup cTWAS 
option_list <- list(
  make_option(
    "--output_folder", action = "store", default = NA, type = 'character',
    help = "Folder for ouput files. [required]"
  ),
  make_option(
    "--output_prefix", action = "store", default = NA, type = 'character',
    help = "The prefix for output files. [required]"
  ),
  make_option(
    "--sumstats", action = "store", default = NA, type = 'character',
    help = "File location for GWAS summary statistics. File suffix of 'RDS' indicates processed Z score. [required]"
  ),
  make_option(
    "--weight_folder", action = "store", default = NA, type = 'character',
    help = "The folder for all the weights stored in sqlite db format. [required]"
  ),
  make_option(
    "--gwas_n", action = "store", default = NA, type = 'integer',
    help = "Sample size of GWAS. [required]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

ld_R_dir <- "/project/mstephens/causalTWAS/UKB_LDR_0.1"
z_snp_file <- opt$sumstats
outname <- opt$output_prefix
outputdir <- opt$output_folder
weight.store <- opt$weight_folder
weightfs <- list.files(weight.store, pattern=".db", full.names=T)

cat("output_folder:",outputdir,"\n")
cat("output_prefix:",outname,"\n")
cat("weights_folder:",weight.store,"\n")
cat("GWAS sumstats:",z_snp_file,"\n")

dir.create(outputdir, showWarnings=F, recursive=T)
gwas_n <- opt$gwas_n
thin <- 0.1
##### Preprocess region_info #####
cat("##### Preprocess region_info ##### \n")
region_info <- readRDS("/project2/xinhe/shared_data/multigroup_ctwas/LD_region_info/region_info.RDS")
LD_map <- readRDS("/project2/xinhe/shared_data/multigroup_ctwas/LD_region_info/LD_map.RDS")
snp_map <- readRDS("/project2/xinhe/shared_data/multigroup_ctwas/LD_region_info/snp_map.RDS")

##### Preprocess GWAS z-scores #####
cat("##### Preprocess z-scores ##### \n")
processed_z_snp_file <- file.path(outputdir, paste0(outname, ".preprocessed.z_snp.RDS"))
if (file.exists(processed_z_snp_file)){
  z_snp <- readRDS(processed_z_snp_file)
} else {
  if (endsWith(toupper(z_snp_file), ".RDS")) {
    file.copy(z_snp_file, processed_z_snp_file)
    z_snp <- readRDS(processed_z_snp_file)
  } else {
    z_snp <- fread(z_snp_file)
    if(is.null(z_snp$BETA)) {z_snp$BETA <- log(z_snp$OR)}
    z_snp$Z <- z_snp$BETA/z_snp$SE
    z_snp <- z_snp[,c("SNP","A1","A2","Z")]
    colnames(z_snp) <- c("id","A1","A2","z")
    # z_snp <- preprocess_z_snp(na.omit(z_snp), snp_map, 
    #                           drop_multiallelic = TRUE, 
    #                           drop_strand_ambig = FALSE)
    saveRDS(z_snp, file = processed_z_snp_file)
  }
}

##### Preprocess weights #####
cat("##### Preprocess weights ##### \n")
processed_weight_file <- file.path(outputdir, paste0(outname, ".preprocessed.weights.RDS"))
if (file.exists(processed_weight_file)){
  cat(sprintf("Load preprocessed weight: %s\n", processed_weight_file))
  weights <- readRDS(processed_weight_file)
}else{
  runtime <- system.time({
    weights <- list()
    for(f in weightfs){
      n <- substr(basename(f),1,nchar(basename(f))-3)
      weights[[n]] <- preprocess_weights(weight_file = f,
                                         region_info = region_info,
                                         gwas_snp_ids = z_snp$id,
                                         snp_map = snp_map,
                                         LD_map = LD_map,
                                         type = "eQTL",
                                         context = n,
                                         weight_format = "PredictDB",
                                         ncore = 10,
                                         drop_strand_ambig = TRUE,
                                         scale_predictdb_weights = FALSE,
                                         load_predictdb_LD = TRUE,
                                         filter_protein_coding_genes = TRUE)  
    }
    
    #weights <- c(weights1,weights2,weights3,weights4,weights5,weights6,weights7,weights8,weights9)   
    weights <- do.call(c, weights)
    saveRDS(weights, file = processed_weight_file)
  })
  cat(sprintf("Preprocessing weights took %0.2f minutes\n",runtime["elapsed"]/60))
}

##### Impute gene z-scores #####
cat("##### Imputing gene z-scores ##### \n")
z_gene_file <- file.path(outputdir, paste0(outname, ".z_gene.RDS"))
if( file.exists(z_gene_file) ){
  loginfo("Load gene z-scores from %s \n", z_gene_file)
  z_gene <- readRDS(z_gene_file)
}else{
  runtime <- system.time({
    z_snp <- readRDS(processed_z_snp_file)
    z_gene <- compute_gene_z(z_snp, weights, ncore=15)
  })
  saveRDS(z_gene, file = z_gene_file)
  loginfo("Imputing gene z-scores took %0.2f minutes\n",runtime["elapsed"]/60)
}

##### Assemble region_data #####
region_data_file <- file.path(outputdir, paste0(outname, ".region_data.RDS"))
if (file.exists(region_data_file)) {
  region_data <- readRDS(region_data_file)
} else{
  runtime <- system.time({
    loginfo("Assemble region_data with thin = %.2f", thin)
    region_data <- assemble_region_data(region_info, 
                                        z_snp, 
                                        z_gene, 
                                        weights,
                                        snp_map,
                                        thin = thin,
                                        ncore = 10)
  })
  loginfo("Assembling region_data took %0.2f minutes\n",runtime["elapsed"]/60)
  saveRDS(region_data, region_data_file)
}

#### Estimate parameters #####
cat("##### Estimating parameters ##### \n")
param_file <- file.path(outputdir, paste0(outname,".param.RDS"))
if (file.exists(param_file)) {
  param <- readRDS(param_file)
} else{
  runtime <- system.time({
    param <- est_param(region_data, 
                       group_prior_var_structure = "shared_type",
                       null_method = "ctwas",
                       niter_prefit = 3,
                       min_gene = 0,
                       min_var = 2,                          
                       min_p_single_effect = 0.8,
                       niter = 50, 
                       ncore = 15,
                       verbose=TRUE)
  })
  saveRDS(param, param_file)
  loginfo("Parameter estimation took %0.2f minutes\n",runtime["elapsed"]/60)
}

group_prior <- param$group_prior
group_prior_var <- param$group_prior_var

### Assess parameter estimates #####
ctwas_parameters <- summarize_param(param, gwas_n)
saveRDS(ctwas_parameters, paste0(outputdir, "/", outname,".parameters.RDS"))

if (thin < 1){
  region_data <- expand_region_data(region_data,
                                    snp_map,
                                    z_snp,
                                    ncore = 6)
}

screened_region_data_file <- file.path(outputdir, paste0(outname,".screened_region_data.RDS"))
screen_summary_file <- file.path(outputdir, paste0(outname,".screened_summary.RDS"))
if (file.exists(screened_region_data_file)) {
  screened_region_data <- readRDS(screened_region_data_file)
} else{
  runtime <- system.time({
    loginfo("screen region with L=1 without LD")
    screen_res <- screen_regions(region_data,
                                 group_prior = group_prior,
                                 group_prior_var = group_prior_var,
                                 min_nonSNP_PIP = 0.5,
                                 ncore = 6)
  })
  loginfo("screening region took %0.2f minutes\n",runtime["elapsed"]/60)
  screened_region_data <- screen_res$screened_region_data
  screen_summary <- screen_res$screen_summary
  saveRDS(screened_region_data, screened_region_data_file)
  saveRDS(screen_summary, screen_summary_file)
}

finemap_regions_file <- file.path(outputdir, paste0(outname,".finemap_regions_res.RDS"))
if (file.exists(finemap_regions_file)) {
  finemap_res <- readRDS(finemap_regions_file)
} else{
  runtime <- system.time({
    finemap_res <- finemap_regions(screened_region_data,
                                   LD_map = LD_map,
                                   weights = weights,
                                   group_prior = group_prior,
                                   group_prior_var = group_prior_var,
                                   null_method = "ctwas",
                                   L = 5,
                                   ncore = 10,
                                   verbose = TRUE,
                                   save_cor = FALSE,
                                   logfile = NULL)
  })
  saveRDS(finemap_res, finemap_regions_file)
  loginfo("Finemapping took %0.2f minutes\n",runtime["elapsed"]/60)
}