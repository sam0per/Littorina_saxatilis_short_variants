rm(list=ls())

pkgs <- c("tools", "ggplot2", "data.table", "optparse", "dplyr", "patchwork", "gridExtra")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(pkgs, require, character.only=TRUE)
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-i", "--vone"), type = "character", default = NULL,
              help = "csv file with allele count for a type of variant and population (e.g., CRAB INDELs or WAVE coding INDEL).",
              metavar = "character"),
  make_option(c("-j", "--vtwo"), type = "character", default = NULL,
              help = "csv file with allele count for another type of variant and population (e.g., WAVE INDELs, CRAB SNPs or WAVE non-coding INDELs).",
              metavar = "character"),
  make_option(opt_str = "--by", type = "character", default = NULL,
              help = "column name that defines what to compare (e.g., VTYPE or ECOT).",
              metavar = "character"),
  make_option(c("-t", "--types"), type = "character", default = NULL,
              help = "two types of comparison (e.g., INDEL:SNP, coding_INDEL:noncoding_INDEL or CRAB:WAVE).",
              metavar = "character"),
  make_option(c("-c", "--csv"), type = "character", default = NULL,
              help = "csv file with other info such as ancestral state.", metavar = "character"),
  make_option(opt_str = "--Test1", action = 'store_true', default = FALSE,
              help = "add this flag to perform Test 1 of Nielsen et al 2005."),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = paste("Rscript scripts/SFS_SP.R -i summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv",
                                        "-j summary/allele_count/AC_CZA_CRAB_SNP_filt2_66N.csv",
                                        "--by variant",
                                        "-t INDEL:SNP",
                                        "-c results/Lsax_short_var_czs_daf_inv_findv.csv"),
                          option_list=option_list,
                          description = "Generate SFS and test whether there is a significant difference.")
opt = parse_args(opt_parser)

if (is.null(opt$vone) | is.null(opt$vtwo) | is.null(opt$by) | is.null(opt$types) | is.null(opt$csv)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_LEFT_SNP_filt2_42N.csv'))
# head(ic)
# sum(duplicated(ic))
ic <- unique(read.csv(file = opt$vone))

# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_SNP_filt2_66N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZA_WAVE_LEFT_INDEL_filt2_59N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_RIGHT_SNP_filt2_42N.csv'))
sc <- unique(read.csv(file = opt$vtwo))

if (unique(ic$N) != unique(sc$N)) {
  stop("The total number of individuals N must be the same.\n", call.=FALSE)
}

# cnm <- 'VTYPE'
# cnm <- 'ECOT'
cnm <- opt$by

# tv <- strsplit(c('INDEL:SNP'), split = ":")[[1]]
# tv <- strsplit(c('CRAB:WAVE_LEFT'), split = ":")[[1]]
# tv <- strsplit(c('WAVE_LEFT:WAVE_RIGHT'), split = ":")[[1]]
tv <- strsplit(opt$types, split = ":")[[1]]

# dt <- unique(read.csv(file = 'results/Lsax_short_var_czs_daf_inv_findv.csv'))
# head(dt)
dt <- unique(read.csv(file = opt$csv))

dtn <- merge(x = dt, y = rbind(ic,sc))
if (opt$rminv) {
  cat('1.  Varinats inside inversions will be removed.\n')
  dtn <- dtn[dtn$invRui==FALSE, ]
}
# head(dtn)
# table(dtn$AC)
# table(dtn$NE_W_Lcomp)
# nrow(unique(dtn))

dtn$DAC <- ifelse(test = dtn$NE_W_Lcomp=='ref_anc:ref_anc', yes = dtn$AC, no = (unique(dtn$N)*2)-dtn$AC)
# head(dtn[dtn$NE_W_Lcomp=='alt_anc:alt_anc', ])
# head(dtn[dtn$DAC==0, ])
# head(dtn[dtn$DAC==(unique(dtn$N)*2), ])
dtp <- dtn[!is.na(dtn$av) & dtn$DAC!=0 & dtn$DAC!=(unique(dtn$N)*2), ]
# table(dtp$DAC)
# hist(dtp$DAC, breaks = 100)
# dtp[dtp$DAF==1,]

# ggplot(data = dtp, aes(x = DAC, fill = VTYPE)) +
#   facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
#   geom_histogram(col = 'black', position = 'dodge', binwidth = 2)
# 
# ggplot(data = dtp, aes(x = DAC, fill = VTYPE)) +
#   facet_wrap(~VTYPE, scales = 'free_y') +
#   geom_histogram(col = 'black', binwidth = 2)


# Suppose we sample 5 individuals, then there are 10 chromosomes
# The SFS is usually expressed in terms of counts of the derived allele
# since variants must be polymorphic, this means 9 possible classes: 1 - 9 copies of the derived allele

# simple example of SNP data

# snp <- c(1000,200,100,50,40,30,20,10,20) # counts of SNPs in each of the 9 classes

snp <- merge(x = data.frame(Var1=1:((unique(dtp$N)*2)-1)), y = data.frame(table(dtp[dtp[, cnm]==tv[2], 'DAC'])),
             by = 'Var1', all.x = TRUE)
# ggplot(data = data.frame(N=1:length(snp), C=snp), aes(x = N, y = C)) +
#   geom_col(col = 'black')

# and indels
# indel <- c(180,30,15,10,8,8,8,9,12) # counts of indels in each of the 9 classes
indel <- merge(x = snp, y = data.frame(table(dtp[dtp[, cnm]==tv[1], 'DAC'])), by = 'Var1', all.x = TRUE)
indel <- apply(X = indel[,-1], MARGIN = 2, FUN = function(x) {
  ifelse(test = is.na(x), yes = 0, no = x)
})
# indel <- indel[order(indel$Var1), 3] + 1
snp <- indel[, 1] + 1
indel <- indel[, 2] + 1
# snp <- snp[order(snp$Var1), 2]


# snp <- ifelse(test = is.na(snp), yes = 0.00001, no = snp)
# ggplot(data = data.frame(N=1:length(indel), C=indel), aes(x = N, y = C)) +
#   geom_col(col = 'black')

DAC_h <- ggplot(data = data.frame(N=rep(1:length(indel), 2), C=c(indel,snp), V=c(rep(tv[1], length(indel)),
                                                                        rep(tv[2], length(indel)))), aes(x = N, y = C)) +
  facet_wrap(facets = ~V, nrow = 2, scales = 'free') +
  geom_col(aes(fill = V), col = 'black', position = 'dodge') +
  scale_fill_manual(values = c("#1B9E77", "#666666")) +
  labs(x = '', y = 'Derived allele count') +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
# DAC_h
# ggsave(filename = 'test/figures/dac_hist.pdf', plot = DAC_h)


# Nielsen et al Test 1 compares the likelihood of the indel data, given the frequencies expected from the SNP data
# to the likelihood of the indel data, based on the indel frequencies
# each likelihood is just a multinomial probability

snp_p <- snp/sum(snp)   # ML frequencies of the classes based on SNP data
indel_p <- indel/sum(indel)   # ML frequencies of the classes based on indel data

# 'difference in proportion': (indels in category/all indels) - (SNPs in category/all SNPs).
dDAP_h <- ggplot(data = data.frame(N=1:length(indel), dP=indel_p-snp_p), aes(x = N, y = dP)) +
  geom_point() +
  labs(x = '', y = paste(tv[1], 'prop. -', tv[2], 'prop.')) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
# dDAP_h

class_contrib <- rep(0,((unique(dtn$N)*2)-1))
if (opt$Test1) {
  
  LL0 <- sum(indel*log(snp_p)) # log-likelihood given expectations based on SNPs

  LL1 <- sum(indel*log(indel_p)) # log-likelihood given expectations based on indels

  # chi_app <- 2*(LL1-LL0) # chi-square approximation with n_classes-1 = 8 df
  # chi_app
  # can examine the contributions of different classes to the overall chisquare

  # class_contrib <- rep(0,9)
  # i=1
  for (i in 1:((unique(dtn$N)*2)-1)){class_contrib[i] <- 2*((indel[i]*log(indel_p[i])+sum(indel[-i])*log(sum(indel_p[-i])))-((indel[i]*log(snp_p[i]))+sum(indel[-i])*log(sum(snp_p[-i]))))}

  # plot(1:9,class_contrib)
  
} else {
  
  # alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)
  
  joint_p <- (snp+indel)/sum(snp+indel)
  
  LL0 <- sum(indel*log(joint_p)) + sum(snp*log(joint_p))
  
  LL1 <- sum(indel*log(indel_p)) + sum(snp*log(snp_p))
  
  a <- rep(0,(unique(dtn$N)*2)-1)
  for (i in 1:((unique(dtn$N)*2)-1)){a[i] <- indel[i]*log(indel_p[i]) + snp[i]*log(snp_p[i]) + sum(indel[-i])*log(sum(indel_p[-i])) + sum(snp[-i])*log(sum(snp_p[-i]))}
  b <- rep(0,(unique(dtn$N)*2)-1)
  for (i in 1:((unique(dtn$N)*2)-1)){b[i] <- (indel[i]*log(joint_p[i])) + (snp[i]*log(joint_p[i])) + sum(indel[-i])*log(sum(joint_p[-i])) + sum(snp[-i])*log(sum(joint_p[-i]))}
  
  # class_contrib <- rep(0,((unique(dtn$N)*2)-1))
  for (i in 1:((unique(dtn$N)*2)-1)){class_contrib[i] <- 2*(a[i]-b[i])}
  
}

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)
cat('Chi-square test statistic =', chi_app, '\n')
# qchisq(p = 0.05, df = 260)
ndf <- ((unique(dtn$N)*2)-1-1)*(length(unique(dtn[, cnm]))-1)
cat('Chi-square critical value at 0.01 with', ndf, 'df =', qchisq(p = 0.01, df = ndf), '\n')
cat('Test statistic - critical value =', chi_app-qchisq(p = 0.01, df = ndf), '\n')

chi_tb <- data.frame(chi2_stat=chi_app,
                     chi2_crit=qchisq(p = 0.01, df = ndf),
                     diff=chi_app-qchisq(p = 0.01, df = ndf))


# plot(1:((unique(dtn$N)*2)-1), class_contrib)

mytheme.b <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7, col = 'black')),
  rowhead = list(fg_params=list(cex = 0.7)))
contr_p <- ggplot(data = data.frame(N=1:((unique(dtn$N)*2)-1), CC=class_contrib)) +
  labs(x = 'Derived allele class', y = 'Contribution') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  annotation_custom(tableGrob(round(chi_tb, 3), theme = mytheme.b, rows = NULL),
                    xmin=((unique(dtn$N)*2)-41), xmax=((unique(dtn$N)*2)-21),
                    ymin = max(class_contrib)-5, ymax=max(class_contrib)) +
  geom_point(aes(x = N, y = CC))

da_cc <- DAC_h / dDAP_h / contr_p +
  plot_layout(heights = c(2, 1, 1))
# da_cc

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv')), split = "_")[[1]]
parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]
if (cnm == 'VTYPE') {
  if (parts[3] == 'WAVE') {
    ECOT <- paste(parts[3], parts[4], sep = '_')
    # VTYPE <- parts[5]
  } else {
    ECOT <- parts[3]
    # VTYPE <- parts[4]
  }
  if (opt$rminv) {
    ggsave(filename = paste('figures/SFS', parts[2], ECOT, tv[1], tv[2] ,'chi_noinv.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  } else {
    ggsave(filename = paste('figures/SFS', parts[2], ECOT, tv[1], tv[2] ,'chi.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  }
} else {
  if (parts[3] == 'WAVE') {
    # ECOT <- paste(parts[3], parts[4], sep = '_')
    VTYPE <- parts[5]
  } else {
    # ECOT <- parts[3]
    VTYPE <- parts[4]
  }
  if (opt$rminv) {
    ggsave(filename = paste('figures/SFS', parts[2], tv[1], tv[2], VTYPE,'chi_noinv.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  } else {
    ggsave(filename = paste('figures/SFS', parts[2], tv[1], tv[2], VTYPE,'chi.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  }
}


# 
# 
# 