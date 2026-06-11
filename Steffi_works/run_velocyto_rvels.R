# Siwei 16 Sept 2025
# Prepare data input for velocyto

# init #####
{
  library(Seurat)
  library(SeuratDisk)
  # library(Signac)
  # library(loupeR)
  # library(scDblFinder)
  # library(monocle)

  library(BiocParallel)
  library(glmGamPoi)
  # library(harmony)
  library(parallel)
  library(future)

  # library(MAST)

  library(stringr)
  library(magrittr)
  library(qs2)

  library(ggplot2)
  library(RColorBrewer)

  # library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(NMF)
  library(rsvd)

  library(velocyto.R)
}

# determine if R is running in RSTUDIO or POSITRON, and set the future plan accordingly
if (Sys.getenv("RSTUDIO") == "1" | Sys.getenv("POSITRON") == "1") {
  print("Running under RStudio IDE, use plan(multisession)")
  session_plan <- "multisession"
} else {
  print("Running under Rscript, use plan(multicore)")
  session_plan <- "multicore"
}

# preload functions ####
get_available_workers <-
  function(x) {
    future::plan(session_plan) # check here!
    return(future::nbrOfFreeWorkers())
  }


## nthreads
workers_2_use <- max(get_available_workers(1) - 1, 1)
workers_2_use <- 16
{
  options(bitmapType = "cairo")
  options(stringsAsFactors = F)
  options(expressions = 20000)
  set.seed(42)
  options(future.globals.maxSize = workers_2_use * 20 * 1024^3) # 40 G per thread
  future::plan(
    session_plan, # Do NOT use "multisession" here if submitting LSF jobs, use "multicore" instead
    workers = workers_2_use
  )
}

# load data #####
gex_merged <-
  qs_read("hashtag_seurat.qs2")
DefaultAssay(gex_merged) <- "RNA"
gex_meta <-
  gex_merged[, sort(as.character(colnames(gex_merged)))]

ldat <-
  read.loom.matrices(
    "outs/velocyto_output/merged_sorted_filtered_bcmatrix_only_5WPDD.loom"
  )

for (i in 1:length(ldat)) {
  print(i)
  colnames(ldat[[i]]) <-
    str_c(
      str_remove_all(
        string = str_split(string = colnames(ldat[[i]]), pattern = ":") %>%
          sapply("[", 2),
        pattern = 'x'
      ),
      "-1"
    )
  ldat[[i]] <-
    ldat[[i]][, sort(as.character(colnames(ldat[[i]])))]
  ldat[[i]] <-
    ldat[[i]][, colnames(gex_meta)]
}

## make a function to use a specific metadata column and return a vector of colours ####
fun_generate_colours <-
  function(
    input_vector = "",
    cell_barcodes = "",
    palette_name = "Dark2"
  ) {
    unique_values <- unique(input_vector)
    palette_colours <-
      brewer.pal(
        n = length(unique_values),
        name = palette_name
      )
    mapped_colour <-
      setNames(palette_colours, nm = unique_values)
    returned_palette <-
      mapped_colour[input_vector]
    # returned_palette <- cell_barcodes
    names(returned_palette) <- cell_barcodes
    return(returned_palette)
  }

cell.colors <-
  fun_generate_colours(
    input_vector = gex_meta$Assignment,
    palette_name = "Dark2",
    cell_barcodes = colnames(gex_meta)
  ) # a vector of many colours, equal to the number of cells

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$ambiguous

# filter expression matrices based on some minimum max-cluster averages
filter_value <- 0.05
emat <-
  filter.genes.by.cluster.expression(
    emat,
    cell.colors,
    min.max.cluster.average = filter_value
  )
nmat <-
  filter.genes.by.cluster.expression(
    nmat,
    cell.colors,
    min.max.cluster.average = filter_value
  )
length(intersect(rownames(emat), rownames(nmat)))
smat <-
  filter.genes.by.cluster.expression(
    smat,
    cell.colors,
    min.max.cluster.average = filter_value
  )
length(intersect(intersect(rownames(emat), rownames(nmat)), rownames(smat))) # 78

# We’ll start with what is perhaps the most robust estimate,
# that combines cell kNN pooling with the gamma fit based on an extreme quantiles:
# Using min/max quantile fit, in which case gene-specific offsets do not
# require spanning read (smat) fit. Here the fit is based on the top/bottom 5%
# of cells (by spliced expression magnitude).
fit.quantile <- 0.05
rvel.qf <-
  gene.relative.velocity.estimates(
    emat = emat,
    nmat = nmat,
    smat = NULL,
    deltaT = 1,
    kCells = 5,
    fit.quantile = fit.quantile
  )

# Calculate velocity use spanning reads (smat) #####
rvel <-
  gene.relative.velocity.estimates(
    emat,
    nmat,
    smat = smat,
    kCells = 5,
    fit.quantile = fit.quantile,
    diagonal.quantiles = TRUE
  )

# calculate the most basic version of velocity estimates, #####
# using relative gamma fit, without cell kNN smoothing
# (i.e. actual single-cell velocity):
rvel1 <-
  gene.relative.velocity.estimates(
    emat,
    nmat,
    deltaT = 1,
    deltaT2 = 1,
    kCells = 1,
    fit.quantile = fit.quantile
  )

qs_savem(
  list(
    rvel.qf = rvel.qf,
    rvel = rvel,
    rvel1 = rvel1
  ),
  "velocity_rvel_estimates.qs2"
)
# emb_tsne <-
#   gex_meta@reductions$tsne@cell.embeddings
