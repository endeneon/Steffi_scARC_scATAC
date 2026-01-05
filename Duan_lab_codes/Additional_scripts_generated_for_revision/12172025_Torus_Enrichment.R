library(mapgen)
library(vroom)
library(bigsnpr)
setwd("~/TORUS/Calculate_new_fdr/")
LD_Blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
head(LD_Blocks, 3)


# bigSNP <- bigsnpr::snp_readBed('~/Data/Databases/GWAS/MAGMA_ref_new/g1k_eur/g1000_eur.bed')
# bigSNP <- snp_attach(bigSNP)
gwas_dir <- "~/TORUS/GWAS_files/"
bed_dir <- "~/TORUS/Calculate_new_fdr/bed_files/"
torus_input_dir <- "./torus_input4"


file_gwas <- list.files(gwas_dir, pattern = "\\.Rdata$", full.names = TRUE)
annotation_bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

for (gwas_file in file_gwas) {
  
  # Load GWAS Rdata
  load(gwas_file)  # assumes it loads object 'gwas'
  trait_name <- tools::file_path_sans_ext(basename(gwas_file))
  
  enrichment <- data.frame()  # initialize empty dataframe
  
  # Loop over BED files
  for (bed_file in annotation_bed_files) {
    
    bed_name <- tools::file_path_sans_ext(basename(bed_file))
    
    # Prepare TORUS input files
    torus.files <- prepare_torus_input_files(gwas.sumstats, bed_file, torus_input_dir = torus_input_dir)
    
    # Run TORUS
    torus.result <- run_torus(
      torus.files$torus_annot_file, 
      torus.files$torus_zscore_file,
      option = "est-prior",
      torus_path = "~/Data/Tools/torus/src/torus"
    )
    
    # Extract enrichment
    torus.enrich <- torus.result$enrich
    tmp <- torus.enrich[torus.enrich$term == paste0(bed_name, ".bed.1"), ]
    tmp$snp <- bed_name
    tmp$trait <- trait_name
    
    enrichment <- rbind(enrichment, tmp)
  }
  
  # Save enrichment for this GWAS
  save_file <- paste0("enrichment_", trait_name, ".rdata")
  save(enrichment, file = save_file)
  cat("Saved:", save_file, "\n")
}

##Torus enrichment for intelligence
load("~/TORUS/GWAS_files/Intelligence2018_re.Rdata")
trait_name <- "Intelligence2018_re"
enrichment <- data.frame() 
for (bed_file in annotation_bed_files) {
  
  bed_name <- tools::file_path_sans_ext(basename(bed_file))
  
  # Prepare TORUS input files
  torus.files <- prepare_torus_input_files(gwas.sumstats, bed_file, torus_input_dir = torus_input_dir)
  
  # Run TORUS
  torus.result <- run_torus(
    torus.files$torus_annot_file, 
    torus.files$torus_zscore_file,
    option = "est-prior",
    torus_path = "~/Data/Tools/torus/src/torus"
  )
  
  # Extract enrichment
  torus.enrich <- torus.result$enrich
  tmp <- torus.enrich[torus.enrich$term == paste0(bed_name, ".bed.1"), ]
  tmp$snp <- bed_name
  tmp$trait <- trait_name
  
  enrichment <- rbind(enrichment, tmp)
}

# Save enrichment for this GWAS
save_file <- paste0("enrichment_", trait_name, ".rdata")
save(enrichment, file = save_file)

##Torus enrichment for T2D
load("~/TORUS/GWAS_files/T2D2018.Rdata")
trait_name <- "T2D2018"
enrichment <- data.frame() 
for (bed_file in annotation_bed_files) {
  
  bed_name <- tools::file_path_sans_ext(basename(bed_file))
  
  # Prepare TORUS input files
  torus.files <- prepare_torus_input_files(gwas.sumstats, bed_file, torus_input_dir = torus_input_dir)
  
  # Run TORUS
  torus.result <- run_torus(
    torus.files$torus_annot_file, 
    torus.files$torus_zscore_file,
    option = "est-prior",
    torus_path = "~/Data/Tools/torus/src/torus"
  )
  
  # Extract enrichment
  torus.enrich <- torus.result$enrich
  tmp <- torus.enrich[torus.enrich$term == paste0(bed_name, ".bed.1"), ]
  tmp$snp <- bed_name
  tmp$trait <- trait_name
  
  enrichment <- rbind(enrichment, tmp)
}

# Save enrichment for this GWAS
save_file <- paste0("enrichment_", trait_name, ".rdata")
save(enrichment, file = save_file)

