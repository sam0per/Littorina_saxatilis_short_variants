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
              help = "csv file with allele count for a type of variant and population (e.g., CRAB INDELs or WAVE coding INDEL).",
              metavar = "character"),
  make_option(c("-j", "--vtwo"), type = "character", default = NULL,
              help = "csv file with allele count for another type of variant and population (e.g., WAVE INDELs, CRAB SNPs or WAVE non-coding INDELs).",
              metavar = "character"),
  make_option(c("-t", "--types"), type = "character", default = NULL,
              help = "two types of comparison (e.g., INDEL:SNP, coding_INDEL:noncoding_INDEL or CRAB:WAVE).",
              metavar = "character"),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = paste("Rscript scripts/SFS_test_short.R",
                                        "-i results/marker_density/MD_CZA_CRAB_coding_INDEL_count.txt",
                                        "-j results/marker_density/MD_CZA_CRAB_noncoding_INDEL_count.txt"),
                          option_list=option_list,
                          description = "Generate SFS from derived allele count and test whether there is a significant difference.")
opt = parse_args(opt_parser)

if (is.null(opt$vone) | is.null(opt$vtwo)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}
# ii <- "results/marker_density/MD_CZA_CRAB_coding_INDEL_count.txt"
# jj <- "results/marker_density/MD_CZA_CRAB_noncoding_INDEL_count.txt"
ii <- opt$vone
jj <- opt$vtwo
ic <- unique(read.table(file = ii, header = TRUE))
sc <- unique(read.table(file = jj, header = TRUE))
if (unique(ic$N) != unique(sc$N)) {
  stop("The total number of individuals N must be the same.\n", call.=FALSE)
}

cnm <- data.frame(II=strsplit(x = basename(ii), split = "_")[[1]],
                  JJ=strsplit(x = basename(jj), split = "_")[[1]])
tv <- c(as.character(cnm[(nrow(cnm)-2):(nrow(cnm)-1), 1]), as.character(cnm[(nrow(cnm)-2):(nrow(cnm)-1), 2]))

jsd <- merge(x = data.frame(Var1=1:((unique(sc$N)*2)-1)), y = data.frame(table(sc$DAC)),
             by = 'Var1', all.x = TRUE)
iid <- merge(x = jsd, y = data.frame(table(ic$DAC)), by = 'Var1', all.x = TRUE)

iid <- apply(X = iid[,-1], MARGIN = 2, FUN = function(x) {
  ifelse(test = is.na(x), yes = 0, no = x)
})

# indel <- indel[order(indel$Var1), 3] + 1
# jsd <- iid[, 1] + 1
jsd <- iid[, 1]
# iid <- iid[, 2] + 1
iid <- iid[, 2]

da <- data.frame(N=rep(1:length(iid), 2), C=c(iid,jsd),
                 V=c(rep(tv[2], length(iid)), rep(tv[4], length(jsd))),
                 A=c(rep(tv[1], length(iid)), rep(tv[3], length(jsd))))
da$AV <- paste(da$A, da$V)
DAC_h <- ggplot(data = da, aes(x = N, y = C)) +
  facet_wrap(facets = ~AV, nrow = 2, scales = 'free') +
  geom_col(aes(fill = AV), col = 'black', position = 'dodge') +
  # scale_fill_manual(values = c("#1B9E77", "#666666")) +
  # scale_fill_manual(values = fac_pal) +
  labs(x = '', y = 'Derived allele count') +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
# DAC_h

jsd_p <- jsd/sum(jsd)   # ML frequencies of the classes based on SNP data
iid_p <- iid/sum(iid)   # ML frequencies of the classes based on indel data
# sum(jsd_p[1:round(10/length(jsd_p)*100)])
# sum(iid_p[1:round(10/length(iid_p)*100)])

# 'difference in proportion': (indels in category/all indels) - (SNPs in category/all SNPs).
dDAP_h <- ggplot(data = data.frame(N=1:length(iid_p), dP=iid_p-jsd_p), aes(x = N, y = dP)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  geom_point() +
  # labs(x = '', y = paste(tv[1], 'prop. -', tv[2], 'prop.')) +
  labs(x = '', y = 'Difference in proportion') +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
# dDAP_h

class_contrib <- rep(0,((unique(sc$N)*2)-1))
# TO AVOID INF AFTER LOG TRANSFORMATION
jsd <- jsd + 1
iid <- iid + 1
jsd_p <- jsd/sum(jsd)
iid_p <- iid/sum(iid)
# alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)

joint_p <- (jsd+iid)/sum(jsd+iid)

LL0 <- sum(iid*log(joint_p)) + sum(jsd*log(joint_p))

LL1 <- sum(iid*log(iid_p)) + sum(jsd*log(jsd_p))

a <- rep(0,(unique(sc$N)*2)-1)
for (i in 1:((unique(sc$N)*2)-1)){a[i] <- iid[i]*log(iid_p[i]) + jsd[i]*log(jsd_p[i]) + sum(iid[-i])*log(sum(iid_p[-i])) + sum(jsd[-i])*log(sum(jsd_p[-i]))}
b <- rep(0,(unique(sc$N)*2)-1)
for (i in 1:((unique(sc$N)*2)-1)){b[i] <- (iid[i]*log(joint_p[i])) + (jsd[i]*log(joint_p[i])) + sum(iid[-i])*log(sum(joint_p[-i])) + sum(jsd[-i])*log(sum(joint_p[-i]))}

# class_contrib <- rep(0,((unique(sc$N)*2)-1))
for (i in 1:((unique(sc$N)*2)-1)){class_contrib[i] <- 2*(a[i]-b[i])}

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)
cat('Chi-square test statistic =', chi_app, '\n')
# qchisq(p = 0.05, df = 260)
ndf <- ((unique(sc$N)*2)-1-1)*(length(unique(da$AV))-1)
cat('Chi-square critical value at 0.01 with', ndf, 'df =', qchisq(p = 0.01, df = ndf), '\n')
cat('Test statistic - critical value =', chi_app-qchisq(p = 0.01, df = ndf), '\n')

chi_tb <- data.frame(chi2_stat=chi_app,
                     chi2_crit=qchisq(p = 0.01, df = ndf),
                     diff=chi_app-qchisq(p = 0.01, df = ndf))
# chi_tb
colnames(chi_tb) <- c("x2 stat", "x2 crit", "x2 stat-crit")

# plot(1:((unique(sc$N)*2)-1), class_contrib)

mytheme.b <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7, col = 'black')),
  rowhead = list(fg_params=list(cex = 0.7)))
contr_p <- ggplot(data = data.frame(N=1:((unique(ic$N)*2)-1), CC=class_contrib)) +
  labs(x = 'Derived allele class', y = 'Contribution') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  annotation_custom(tableGrob(round(chi_tb, 3), theme = mytheme.b, rows = NULL),
                    xmin=((unique(ic$N)*2)-(unique(ic$N)*2)*0.4), xmax=((unique(ic$N)*2)-(unique(ic$N)*2)*0.2),
                    ymin = max(class_contrib)-(max(class_contrib)*0.6), ymax=max(class_contrib)) +
  geom_point(aes(x = N, y = CC))


da_cc <- DAC_h / dDAP_h / contr_p +
  plot_layout(heights = c(2, 1, 1))
# da_cc

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]
parts <- cnm[c(-1, (-nrow(cnm)+2):-nrow(cnm)), ]
if (identical(as.character(parts[,1]),as.character(parts[,2]))) {
  parts <- paste(as.character(parts[,1]), collapse = '_')
} else {
  parts <- paste(paste(as.character(parts[,1]), collapse = '_'), paste(as.character(parts[,2]), collapse = '_'), sep = '_')
}
ggsave(filename = paste('figures/SFS', parts, paste(tv, collapse = '_'),
                        'chi_noinv.pdf', sep = "_"), plot = da_cc, scale = 0.75, dpi = "screen")


