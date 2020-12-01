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
  make_option(opt_str = "--test1", action = 'store_true', default = FALSE,
              help = "add this flag to perform Nielsen et al. 2005 Test1-like test."))

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
# ii <- "results/marker_density/MD_CZA_WAVE_RIGHT_noncoding_INDEL_count.txt"
# jj <- "results/marker_density/MD_CZA_WAVE_RIGHT_noncoding_SNP_count.txt"
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

fac_pal <- as.vector(x = rep("NA", 2), mode = "character")
if (grepl(pattern = "IN|DE", x = tv[2])) {
  fac_pal[1] <- brewer.pal(n = 3, name = "Paired")[1]
} else {
  fac_pal[1] <- brewer.pal(n = 3, name = "Paired")[3]
}

if (grepl(pattern = "IN|DE", x = tv[4])) {
  fac_pal[2] <- brewer.pal(n = 4, name = "Paired")[2]
} else if (paste(tv[3:4], collapse = "_")=="coding_SNP") {
  fac_pal[2] <- brewer.pal(n = 4, name = "Paired")[3]
} else {
  fac_pal[2] <- brewer.pal(n = 4, name = "Paired")[4]
}

DAC_h <- ggplot(data = da, aes(x = N, y = C)) +
  facet_wrap(facets = ~AV, nrow = 2, scales = 'free') +
  geom_col(aes(fill = AV), col = 'black', position = 'dodge') +
  # scale_fill_manual(values = c("#1B9E77", "#666666")) +
  scale_fill_manual(values = fac_pal) +
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

parts <- cnm[c(-1, (-nrow(cnm)+2):-nrow(cnm)), ]
if (identical(as.character(parts[,1]),as.character(parts[,2]))) {
  parts <- paste(as.character(parts[,1]), collapse = '_')
} else {
  parts <- paste(paste(as.character(parts[,1]), collapse = '_'), paste(as.character(parts[,2]), collapse = '_'), sep = '_')
}

i0 <- which(iid==0)
j0 <- which(jsd==0)
if (length(i0)!=0 | length(j0)!=0) {
  iid <- iid[-unique(c(i0, j0))]
  jsd <- jsd[-unique(c(i0, j0))]
}

jsd_p <- jsd/sum(jsd)   # ML frequencies of the classes based on SNP data
iid_p <- iid/sum(iid)   # ML frequencies of the classes based on indel data

class_contrib <- rep(0, length(iid))
# alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)

joint_p <- (jsd+iid)/sum(jsd+iid)

LL0 <- sum(iid*log(joint_p)) + sum(jsd*log(joint_p))

LL1 <- sum(iid*log(iid_p)) + sum(jsd*log(jsd_p))

a <- rep(0,length(iid))
for (i in 1:length(iid)){a[i] <- iid[i]*log(iid_p[i]) + jsd[i]*log(jsd_p[i]) + sum(iid[-i])*log(sum(iid_p[-i])) + sum(jsd[-i])*log(sum(jsd_p[-i]))}
b <- rep(0,length(iid))
for (i in 1:length(iid)){b[i] <- (iid[i]*log(joint_p[i])) + (jsd[i]*log(joint_p[i])) + sum(iid[-i])*log(sum(joint_p[-i])) + sum(jsd[-i])*log(sum(joint_p[-i]))}

# class_contrib <- rep(0,((unique(sc$N)*2)-1))
for (i in 1:length(iid)){class_contrib[i] <- 2*(a[i]-b[i])}

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)
cat('Chi-square test statistic =', chi_app, '\n')
# qchisq(p = 0.05, df = 260)
ndf <- (length(iid)-1)*(length(unique(da$AV))-1)
pval <- pchisq(q = chi_app, df = ndf)
# cat('Chi-square critical value at 0.01 with', ndf, 'df =', qchisq(p = 0.01, df = ndf), '\n')
# cat('Test statistic - critical value =', chi_app-qchisq(p = 0.01, df = ndf), '\n')

# chi_tb <- data.frame(chi2_stat=chi_app,
#                      chi2_crit=qchisq(p = 0.01, df = ndf),
#                      diff=chi_app-qchisq(p = 0.01, df = ndf))
# chi_tb
# colnames(chi_tb) <- c("x2 stat", "x2 crit", "x2 stat-crit")

mytheme.b <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7, col = 'black')),
  rowhead = list(fg_params=list(cex = 0.7)))

contr_d <- data.frame(N=1:((unique(ic$N)*2)-1), CC = 0)
contr_d$CC[-unique(c(i0,j0))] <- class_contrib

Nc <- nrow(contr_d)
# tmp_iid <- vector(mode = "integer", length = Nc)
# tmp_jsd <- vector(mode = "integer", length = Nc)

# tmp_iid[-c(i0,j0)] <- iid
# tmp_jsd[-c(i0,j0)] <- jsd
iid <- da$C[1:Nc]
jsd <- da$C[(Nc+1):nrow(da)]

if (grepl(pattern = "noncoding_INDEL", x = jj)) {
  
  if (51<Nc) {
    brk <- c(0:51, length(jsd))
    bn <- cut(x = 1:Nc, breaks = brk, include.lowest = TRUE)
    new_ij <- data.frame(iid, jsd, Bin=bn)
    
    iid <- aggregate(x = new_ij$iid, by = list(new_ij$Bin), sum)$x
    jsd <- aggregate(x = new_ij$jsd, by = list(new_ij$Bin), sum)$x
  }
  
} else if (grepl(pattern = "_coding_SNP", x = jj)) {
  
  brk <- c(0:18, 23, 28, 33, length(jsd))
  bn <- cut(x = 1:Nc, breaks = brk, include.lowest = TRUE)
  new_ij <- data.frame(iid, jsd, Bin=bn)
  
  iid <- aggregate(x = new_ij$iid, by = list(new_ij$Bin), sum)$x
  jsd <- aggregate(x = new_ij$jsd, by = list(new_ij$Bin), sum)$x
  
}

mx <- matrix(nrow = length(iid), ncol = 16)
colnames(mx) <- c("Population", "SFS", "Chi-squared", "p-value", paste("observed", tv[1], tv[2], sep = '_'),
                  paste("observed", tv[3], tv[4], sep = '_'), paste("sum", tv[1], tv[2], sep = '_'),
                  paste("sum", tv[3], tv[4], sep = '_'), paste("expected", tv[1], tv[2], sep = '_'),
                  paste("expected", tv[3], tv[4], sep = '_'), paste("sum_expected", tv[1], tv[2], sep = '_'),
                  paste("sum_expected", tv[3], tv[4], sep = '_'), paste("residuals_obs", tv[1], tv[2], sep = '_'),
                  paste("residuals_obs", tv[3], tv[4], sep = '_'), paste("residuals_sum", tv[1], tv[2], sep = '_'),
                  paste("residuals_sum", tv[3], tv[4], sep = '_'))
for (i in 1:length(iid)) {
  # i <- 41
  tt <- data.frame(c=c(iid[i],jsd[i]), t=c(sum(iid[-i]),sum(jsd[-i])))
  cs <- chisq.test(tt)
  mx[i,1] <- parts
  mx[i,2] <- paste(tv, collapse = '_')
  mx[i,3] <- round(cs$statistic, 4)
  mx[i,4] <- round(cs$p.value, 4)
  mx[i,5] <- iid[i]
  mx[i,6] <- jsd[i]
  mx[i,7] <- sum(iid[-i])
  mx[i,8] <- sum(jsd[-i])
  mx[i,9] <- round(cs$expected[1,1])
  mx[i,10] <- round(cs$expected[2,1])
  mx[i,11] <- round(cs$expected[1,2])
  mx[i,12] <- round(cs$expected[2,2])
  mx[i,13] <- round(cs$residuals[1,1], 4)
  mx[i,14] <- round(cs$residuals[2,1], 4)
  mx[i,15] <- round(cs$residuals[1,2], 4)
  mx[i,16] <- round(cs$residuals[2,2], 4)
}
mx <- as.data.frame(mx)
mx$`p-value` <- as.numeric(as.character(mx$`p-value`))
mx$`Chi-squared` <- as.numeric(as.character(mx$`Chi-squared`))
mx$signif <- ifelse(test = mx$`p-value` < 0.05, yes = TRUE, no = FALSE)
mx$pal <- ifelse(test = mx$signif == TRUE, yes = "red", no = "black")
mx$N <- 1:nrow(mx)
# colnames(contr_d)

if (exists("bn")) {
  
  contr_d$Bin <- as.character(bn)
  mx$Bin <- unique(as.character(bn))
  contr_d <- merge(x = contr_d, y = mx[, c("Bin", "Chi-squared", "pal")], by = "Bin", all.x = TRUE)
  contr_d <- contr_d[order(contr_d$N),]
  
} else {
  
  contr_d <- merge(x = contr_d, y = mx[, c("N", "Chi-squared", "pal")], all.x = TRUE)
  
}

# contr_d$pal <- ifelse(test = is.na(contr_d$pal), yes = mx$pal[nrow(mx)], no = contr_d$pal)
# table(mx$signif)
# mx[mx$signif==TRUE,]

write.table(x = mx, file = paste('results/SFS_comp/SFS', parts, paste(tv, collapse = '_'),
                                 'chi_noinv.csv', sep = "_"), append = FALSE, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = TRUE)

contr_p <- ggplot(data = contr_d) +
  labs(x = 'Derived allele class', y = 'Contribution') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  # annotation_custom(tableGrob(round(chi_tb, 3), theme = mytheme.b, rows = NULL),
  #                   xmin=(unique(ic$N)*2)-(unique(ic$N)*2)*0.2, xmax=((unique(ic$N)*2)-(unique(ic$N)*2)*0.2),
  #                   ymin = max(mx$`Chi-squared`)-(max(mx$`Chi-squared`)*0.6), ymax=max(mx$`Chi-squared`)) +
  annotate(geom="text", x=Nc*0.8, y=max(mx$`Chi-squared`)*0.8, label=paste0("p-value = ", round(pval, 2)), color="black") +
  geom_point(aes(x = N, y = `Chi-squared`), col = contr_d$pal)

# plot(1:((unique(sc$N)*2)-1), class_contrib)

da_cc <- DAC_h / dDAP_h / contr_p +
  plot_layout(heights = c(2, 1, 1))
# da_cc

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]

if (opt$test1) {
  ggsave(filename = paste('figures/SFS', parts, paste(tv, collapse = '_'),
                          'chi_noinv.pdf', sep = "_"), plot = da_cc, scale = 0.75, dpi = "screen")
} else {
  ggsave(filename = paste('figures/SFS_comp/SFS', parts, paste(tv, collapse = '_'),
                          'chi_noinv.pdf', sep = "_"), plot = da_cc, scale = 0.75, dpi = "screen")
}



