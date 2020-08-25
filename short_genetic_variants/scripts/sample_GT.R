rm(list=ls())

.packages = c("tools", "tidyr", "data.table", "optparse", "Rmisc", "dplyr")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-g", "--genotypes"), type = "character", default = NULL,
              help = "Genotypes encoded as 0,1,2 in a individuals by variant table.", metavar = "character"),
  make_option(c("-n", "--number"), type = "integer", default = NULL,
              help = "Indicate the number of individuals/genotypes to sample.", metavar = "integer"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Downsample genotypes and count the number of alleles.",
                          epilogue = "Example: Rscript scripts/sample_GT.R -g geno_matrix/CRAB/GM_CZA_CRAB_INDEL.filt2-1022.012 -n 23")
opt = parse_args(opt_parser)

if (is.null(opt$genotypes) | is.null(opt$number)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# infl <- 'test/GM_CZA_CRAB_INDEL.filt2-1022.012'
# gt <- read.table(file = infl, sep = '\t', row.names = 1)
infl <- opt$genotypes
gt = read.table(infl, sep = '\t', row.names = 1)
indv = read.table(paste0(infl, '.indv'), sep = '\t')
pos = read.table(paste0(infl, '.pos'), sep = '\t')
colnames(pos) <- c("CHROM", "POS")
# n_smp = 20
n_smp = opt$number

# gt_sum <- apply(X = gt, MARGIN = 2, FUN = sum)
# gt_sum

file.exists("out_sample_GT.log")
success <- FALSE
i <- 1
while (!success) {
  if (i == 30) {
    mgt <- sum(gt_smp < 0)
    cat(paste(infl, ': there is/are', mgt, 'missing genotype(s).'), file = 'out_sample_GT.log', sep = '\n', append = TRUE)
    success <- TRUE
  } else {
    idx_smp <- replicate(n = ncol(gt), expr = sample(x = 1:nrow(gt), size = n_smp))
    gt_smp <- sapply(X = 1:ncol(gt), FUN = function(x) gt[idx_smp[, x], x])
    success <- sum(gt_smp < 0) == 0
    i = i + 1
  }
  
}
# gt[13,4]
# gt_smp
gtsmp_sum <- apply(X = gt_smp, MARGIN = 2, FUN = sum)
# gtsmp_sum

indv_smp <- data.frame(apply(X = idx_smp, MARGIN = 2, FUN = function(x) indv[x, 1]))
write.table(x = indv_smp, file = paste0(infl, '.samples.csv'), quote = FALSE, sep = ',', row.names = FALSE, col.names = FALSE)

# indv[c(54,33,27,61),1]
# indv[c(62,35,11,4),1]

parts <- strsplit(file_path_sans_ext(basename(infl)), split = "_")[[1]]
parts <- c(parts[-length(parts)], strsplit(parts[length(parts)], split = "[.]")[[1]][1])

out_dt <- data.frame(pos, cp = paste(pos[,1], pos[,2], sep = '_'), AC = gtsmp_sum, ZONE = parts[2])

if (parts[3] == 'WAVE') {
  out_dt$ECOT <- paste(parts[3], parts[4], sep = '_')
  out_dt$VTYPE <- parts[5]
} else {
  out_dt$ECOT <- parts[3]
  out_dt$VTYPE <- parts[4]
}

write.csv(x = out_dt, file = paste0(infl, '.', n_smp, 'N.csv'), quote = FALSE, row.names = FALSE)
