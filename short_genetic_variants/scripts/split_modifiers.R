rm(list=ls())

pkgs <- c("optparse")
# pkgs <- c("tools", "ggplot2", "data.table", "optparse", "dplyr", "patchwork", "gridExtra", "RColorBrewer")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "table of variants annotated by SnpEff as MODIFIER.",
              metavar = "character"))

opt_parser = OptionParser(usage = "Rscript scripts/split_modifiers.R -f results/marker_density/MD_CZA_CRAB_nongenic_SW_count.txt",
                          option_list=option_list,
                          description = "Generate input file for the software DH.")
opt = parse_args(opt_parser)

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("Missing input file.\n", call.=FALSE)
}

# dt <- read.table(file = 'results/marker_density/MD_CZA_CRAB_nongenic_SW_count.txt', header = TRUE)
dt <- read.table(file = opt$file, header = TRUE)
# head(dt)
# length(unique(dt$cp))

dt$FUN <- NA
for (i in 1:nrow(dt)) {
  stsp <- strsplit(x = as.character(dt$ANN[i]), split = '\\|')[[1]][2]
  # uim <- paste(unique(stsp[stsp %in% imp]), collapse = ':')
  dt$FUN[i] <- stsp
}
dt$FUN <- as.character(dt$FUN)
tb <- data.frame(table(dt$FUN))



lapply(X = as.character(tb$Var1), FUN = function(x) {
  dd <- dt[dt$FUN==x, ]
  an_vt <- paste(x, unique(dd$VTYPE), sep = '_')
  
  # write.table(x = dtp, file = paste('results/marker_density/MD', unique(dtp$ZONE), levels(ic$ECOT),
  #                                   an_vt, 'count.txt', sep = '_'),
  #             quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  dac <- dd[, 'DAC']
  df <- data.frame(matrix(data = 0, nrow = unique(dd$N)*2, ncol = length(dac),
                          dimnames = list(1:(unique(dd$N)*2), as.character(dd[, 'cp']))),
                   stringsAsFactors = FALSE)
  
  for (i in 1:ncol(df)) {
    df[, i] <- c(rep(1, dac[i]), rep(0, (unique(dd$N)*2)-dac[i]))
  }
  # head(df)
  # colSums(df)
  
  fileConn <- file(paste('summary/haplotypes/HAP', levels(dd$ZONE), levels(dd$ECOT),
                         an_vt, 'DH.txt', sep = '_'))
  writeLines(c('# ', an_vt, '\n//\nsegsites: ', length(dac),
               '\npositions: ', paste(round(1:length(dac)/length(dac), 4), collapse = ' '), '\n'),
             con = fileConn, sep = '')
  close(fileConn)
  
  write.table(x = df,
              file = paste('summary/haplotypes/HAP', levels(dd$ZONE), levels(dd$ECOT),
                           an_vt, 'DH.txt', sep = '_'),
              quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE, append = TRUE)
  
  if (!file.exists('data/MD_nongenic_combinations.csv')) {
    file.create('data/MD_nongenic_combinations.csv')
  }
  
  write.table(x = data.frame(paste('summary/haplotypes/HAP', levels(dd$ZONE), levels(dd$ECOT),
                                   an_vt, 'DH.txt', sep = '_'),
                             unique(dd$N)*2), file = 'data/MD_nongenic_combinations.csv', append = TRUE, quote = FALSE,
              sep = ',', row.names = FALSE, col.names = FALSE)
})
# dt[dt$cp=='Contig9943_30730',]
# head(dt[dt$FUN=='3_prime_UTR_variant', ])
# dt[dt$FUN=='3_prime_UTR_variant', ][6, 'DAC']

cat('MISSION COMPLETE!\n')
