rm(list=ls())
# setwd("Anja/Anja_results/20200115/")
packs = c("mgsub", "tools", "data.table", "optparse", "dplyr")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(packs, require, character.only=TRUE)
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="DIRECTORY with input data", metavar="character"),
  make_option(c("-z", "--zone"), type="character", default=NULL,
              help="Name of the ZONE found in the input filenames [example: CZA]", metavar="character"),
  make_option(c("-v", "--variant"), type="character", default=NULL,
              help="VARIANT TYPE, either SNP or INDEL", metavar="character"),
  make_option(c("-i", "--indv"), type="character", default=NULL,
              help="List of INDIVIDUALS ids", metavar="character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Infer ancestral vs derived alleles.",
                          epilogue = "Example: Rscript scripts/ancestry_allele.R -d geno_matrix -z CZA -v INDEL -i CZCLI01_CZA_INDEL.filt2.012.id")
opt = parse_args(opt_parser)

if (is.null(opt$dir) | is.null(opt$zone) | is.null(opt$variant) | is.null(opt$indv)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# idir <- "geno_matrix"
idir <- opt$dir
# zone = "CZA"
zone = opt$zone
# vartype = "INDEL"
vartype = opt$variant
# indv <- read.csv("geno_matrix/CZCLI01_CZA_INDEL.filt2.012.id")
indv <- read.csv(opt$indv)

gt_fl <- list.files(path = idir, pattern = paste0(zone, '_', vartype, '.*.complete'), full.names = TRUE)
cat('\nSome of the files that will be analysed:\n')
head(gt_fl)

gt_ls <- lapply(gt_fl, read.table, header=TRUE)
cat('\nTotal number of input genomic intervals:', length(gt_ls))
# lapply(gt_ls, head)
gt_ls <- gt_ls[sapply(gt_ls, function(x) dim(x)[2]) > 1]
cat('\nTotal number of input genomic intervals with called genotypes:', length(gt_ls))
# lapply(gt_ls, tail)
cat('\n')

cat('\n... inferring allele ancestry using', as.character(tail(indv$snail_ID, n=2)), 'as outgroup(s) ...\n')
gt_anc <- rbindlist(lapply(gt_ls, FUN = function(x) {
  outg_gt <- data.frame(tail(x[, -1], n=2))
  colnames(outg_gt) <- colnames(x)[-1]
  rownames(outg_gt) <- as.character(tail(indv$snail_ID, n=2))
  
  anc <- t(apply(X = outg_gt, MARGIN = 2, FUN = function(y) {
    mgsub(string = as.character(y), pattern = as.character(c(0, 1, 2)),
          replacement = c("ref_anc", "het", "alt_anc"), recycle = FALSE)
  }))
  
  out_dt <- data.frame(cp=rownames(anc), anc, zone, vartype)
  colnames(out_dt) <- c("cp", rownames(outg_gt), "ZONE", "VTYPE")
  return(out_dt)
}))

cat('\n... writing output', paste0('ANC_GT_', zone, '_', vartype, '.csv'), 'which contains',
    nrow(gt_anc), vartype, 'variants ...\n')
write.csv(x = data.frame(gt_anc), file = paste0('ANC_GT_', zone, '_', vartype, '.csv'), row.names = FALSE)

cat('\nNumber of variants that are homozygous for the reference allele in the outgroup(s) = ref_anc\n')
cat('Number of variants that are homozygous for the alternative allele in the outgroup(s) = alt_anc\n')
cat('Number of variants that are heterozygous in the outgroup(s) = het\n')
sapply(X = as.character(tail(indv$snail_ID, n=2)), FUN = function(i) {
  # i
  table(data.frame(gt_anc)[, i])
})

cat('\nMISSION COMPLETE!\n')
