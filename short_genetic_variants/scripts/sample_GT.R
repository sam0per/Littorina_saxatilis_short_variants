rm(list=ls())

pkgs <- c("tools", "data.table", "optparse", "dplyr")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(pkgs, require, character.only=TRUE)
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
# n_smp = 66
n_smp = opt$number

file.exists("out_sample_GT.log")
if (nrow(gt) < n_smp) {
  
  cat(paste(infl, ':', nrow(gt), 'individuals are less than', n_smp, 'filter cutoff.'),
      file = 'out_sample_GT.log', sep = '\n', append = TRUE)

} else if (ncol(gt) == 0) {
  
  cat(paste(infl, ': no information available.'),
      file = 'out_sample_GT.log', sep = '\n', append = TRUE)

} else {
  
  gt_smp <- apply(X = gt, MARGIN = 2, FUN = function(x) {
    pos_idx <- which(x >= 0)
    if (length(pos_idx) >= n_smp) {
      pindv <- as.character(indv[pos_idx,])
      pool <- x[x >= 0]
      idx_smp <- sample(x = 1:length(pool), size = n_smp, replace = FALSE)
      data.frame(GT = pool[idx_smp], ID = pindv[idx_smp])
    }
  })
  
  gt_smp <- Filter(function(x) length(x) > 0, gt_smp)
  pos <- pos[which(colnames(gt) %in% names(gt_smp)), ]
  
  gt_smp <- do.call(what = 'cbind', args = gt_smp)
  
  if (ncol(gt_smp) == 2) {
    gtsmp_sum <- sum(gt_smp[, grepl(pattern = 'GT', x = colnames(gt_smp))])
  } else {
    gtsmp_sum <- apply(X = gt_smp[, grepl(pattern = 'GT', x = colnames(gt_smp))], MARGIN = 2, FUN = sum)
  }
  
  indv_smp <- data.frame(ID = gt_smp[, grepl(pattern = 'ID', x = colnames(gt_smp))])
  write.table(x = indv_smp, file = paste0(infl, '.samples.csv'), quote = FALSE, sep = ',', row.names = FALSE, col.names = FALSE)
  
  parts <- strsplit(file_path_sans_ext(basename(infl)), split = "_")[[1]]
  parts <- c(parts[-length(parts)], strsplit(parts[length(parts)], split = "[.]")[[1]][1])
  
  out_dt <- data.frame(pos, cp = paste(pos[,1], pos[,2], sep = '_'), AC = gtsmp_sum, N = n_smp, ZONE = parts[2])
  
  if (parts[3] == 'WAVE') {
    out_dt$ECOT <- paste(parts[3], parts[4], sep = '_')
    out_dt$VTYPE <- parts[5]
  } else {
    out_dt$ECOT <- parts[3]
    out_dt$VTYPE <- parts[4]
  }
  
  write.csv(x = out_dt, file = paste0(infl, '.', n_smp, 'N.csv'), quote = FALSE, row.names = FALSE)
  
  cat('MISSION COMPLETED.\n')
}
