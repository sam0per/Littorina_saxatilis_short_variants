rm(list=ls())

pkgs <- c("tools", "ggplot2", "data.table", "optparse", "dplyr", "patchwork", "gridExtra", "RColorBrewer")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-i", "--vone"), type = "character", default = NULL,
              help = "csv file with allele count for a type of variant and population (e.g., CRAB INDELs or WAVE SNPs).",
              metavar = "character"),
  make_option(c("-c", "--csv"), type = "character", default = NULL,
              help = "csv file with other info such as ancestral state.", metavar = "character"),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = paste("Rscript scripts/summary_stats.R -i summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv",
                                        "-c results/Lsax_short_var_czs_daf_inv_findv.csv"),
                          option_list=option_list,
                          description = "Perform Mann-Whitney U tests for a difference between the WW, WS, SW and SS DAF spectra.")
opt = parse_args(opt_parser)

if (is.null(opt$vone) | is.null(opt$csv)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_WAVE_LEFT_SNP_filt2_59N.csv'))
# head(ic)
# sum(duplicated(ic))
ic <- unique(read.csv(file = opt$vone))

# dt <- unique(read.csv(file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv'))
# head(dt)
dt <- unique(read.csv(file = opt$csv))

ann_fl <- list.files(path = 'annotated', pattern = paste(levels(ic$ZONE), 'SNP',
                                                           sep = '_'), full.names = TRUE)
ann <- read.table(file = ann_fl, header = TRUE)
ann$VTYPE <- 'SNP'
  
ann$ZONE <- levels(ic$ZONE)
ann$cp <- paste(ann$CHROM, ann$POS, sep = '_')

dtn <- unique(merge(x = dt, y = ic))
dtn <- unique(merge(x = dtn, y = ann))

# head(dtn)
# table(dtn$ZONE)
# table(dtn$VTYPE)
# table(dtn$ECOT)
# table(dtn$CLASS)
# table(dtn$gBGC) 

if (opt$rminv) {
  cat('1.  Varinats inside inversions will be removed.\n')
  dtn <- dtn[dtn$invRui==FALSE, ]
}

# table(dtn$NE_W_Lcomp)
# nrow(unique(dtn))

dtn$DAC <- ifelse(test = dtn$NE_W_Lcomp=='ref_anc:ref_anc', yes = dtn$AC, no = (unique(dtn$N)*2)-dtn$AC)

dtp <- dtn[!is.na(dtn$av) & dtn$DAC!=0 & dtn$DAC!=(unique(dtn$N)*2), ]
# dtp[duplicated(dtp$cp), ]
dtp <- dtp[!duplicated(dtp$cp), ]

gc <- combn(x = levels(dtp$gBGC), m = 2)

apply(X = gc, MARGIN = 2, FUN = function(x) {
  xr <- length(unique(dtp[dtp$gBGC==x[1], 'cp']))
  xc <- length(unique(dtp[dtp$gBGC==x[2], 'cp']))
  
  sfsr <- merge(x = data.frame(Var1 = 1:((unique(dtp$N)*2)-1)), y = data.frame(table(dtp[dtp$gBGC==x[1], 'DAC'])),
                by = 'Var1', all.x = TRUE)
  sfsa <- merge(x = sfsr, y = data.frame(table(dtp[dtp$gBGC==x[2], 'DAC'])), by = 'Var1', all.x = TRUE)
  sfsa <- apply(X = sfsa[,-1], MARGIN = 2, FUN = function(x) {
    ifelse(test = is.na(x), yes = 0, no = x)
  })
  colnames(sfsa) <- x
  # sfsa[, 1] <- sfsa[, 1] + 1
  # sfsa[, 2] <- sfsa[, 2] + 1
  # assign(x[1], sfsa[, 1] + 1)
  # assign(x[2], sfsa[, 2] + 1)
  
  wt <- wilcox.test(x = sfsa[,1], y = sfsa[,2], alternative = "two.sided", paired = FALSE)
  # normalized U statistic
  # dividing the U statistic by its maximum possible value in the test, the product of the number of SNPs in the two categories (Bamber 1975)
  # apply(X = sfsa, MARGIN = 2, sum)
  nwt <- wt$statistic/(xr * xc)
  cat(x, ': U-norm =', nwt, '\n', sep = ' ')
  
  ze <- paste(unique(dtp$ZONE), unique(dtp$ECOT), sep = '_')
  sfsa <- data.frame(N=rep(1:dim(sfsa)[1], 2), C=c(sfsa[,1], sfsa[,2]),
                     V=c(rep(x[1], dim(sfsa)[1]), rep(x[2], dim(sfsa)[1])),
                     TOT=c(rep(xr, dim(sfsa)[1]), rep(xc, dim(sfsa)[1])),
                     ZE=ze)
  sfsa$P <- sfsa$C/sfsa$TOT
  ggplot(data = sfsa, aes(x = N, y = P, width = .7)) +
    geom_col(aes(fill = V), position = 'dodge') +
    annotate(geom = 'text', x=max(sfsa$N)-15, y=max(sfsa$P)-0.02, label=paste(gsub(pattern = '_', replacement = ' ', x = ze),
                                                                              'U-norm =', round(nwt,5), sep = ' ')) +
    labs(fill = '', x = 'Derived allele frequency', y = 'Proportion') +
    scale_fill_manual(values = brewer.pal(n = 12, name = "Paired")[6:7]) +
    theme(legend.position = 'top',
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))

})
