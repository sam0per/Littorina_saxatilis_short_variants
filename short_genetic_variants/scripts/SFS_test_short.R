rm(list = setdiff(ls(), c('parts', 'cnm')))

# tv <- strsplit(c('syn_INDEL:syn_SNP'), split = ":")[[1]]
tv <- strsplit(c('nongenic_INDEL:nongenic_SNP'), split = ":")[[1]]

tv2 <- unlist(strsplit(tv, split = "_"))


# sc <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_SNP_syn_count.txt', header = TRUE, sep = '\t')
sc <- read.table(file = 'results/marker_density/MD_CZD_CRAB_syn_WWSS_count.txt', header = TRUE, sep = '\t')
head(sc)
# ic <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_INDEL_syn_count.txt', header = TRUE, sep = '\t')
ic <- read.table(file = 'results/marker_density/MD_CZA_CRAB_nonsyn_DEL_count.txt', header = TRUE, sep = '\t')
head(ic)

snp <- merge(x = data.frame(Var1=1:((unique(sc$N)*2)-1)), y = data.frame(table(sc$DAC)),
             by = 'Var1', all.x = TRUE)
indel <- merge(x = snp, y = data.frame(table(ic$DAC)), by = 'Var1', all.x = TRUE)

indel <- apply(X = indel[,-1], MARGIN = 2, FUN = function(x) {
  ifelse(test = is.na(x), yes = 0, no = x)
})

# indel <- indel[order(indel$Var1), 3] + 1
snp <- indel[, 1] + 1
indel <- indel[, 2] + 1

DAC_h <- ggplot(data = data.frame(N=rep(1:length(indel), 2), C=c(indel,snp), V=c(rep(tv[1], length(indel)),
                                                                        rep(tv[2], length(indel)))),
                aes(x = N, y = C-1)) +
  facet_wrap(facets = ~V, nrow = 2, scales = 'free') +
  geom_col(aes(fill = V), col = 'black', position = 'dodge') +
  # scale_fill_manual(values = c("#1B9E77", "#666666")) +
  # scale_fill_manual(values = fac_pal) +
  labs(x = '', y = 'Derived allele count') +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
DAC_h

snp_p <- snp/sum(snp)   # ML frequencies of the classes based on SNP data
indel_p <- indel/sum(indel)   # ML frequencies of the classes based on indel data
sum(snp_p[1:round(10/length(indel)*100)])
sum(indel_p[1:round(10/length(indel)*100)])
indel_p-snp_p

# 'difference in proportion': (indels in category/all indels) - (SNPs in category/all SNPs).
dDAP_h <- ggplot(data = data.frame(N=1:length(indel), dP=indel_p-snp_p), aes(x = N, y = dP)) +
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
dDAP_h

class_contrib <- rep(0,((unique(sc$N)*2)-1))

# alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)

joint_p <- (snp+indel)/sum(snp+indel)

LL0 <- sum(indel*log(joint_p)) + sum(snp*log(joint_p))

LL1 <- sum(indel*log(indel_p)) + sum(snp*log(snp_p))

a <- rep(0,(unique(sc$N)*2)-1)
for (i in 1:((unique(sc$N)*2)-1)){a[i] <- indel[i]*log(indel_p[i]) + snp[i]*log(snp_p[i]) + sum(indel[-i])*log(sum(indel_p[-i])) + sum(snp[-i])*log(sum(snp_p[-i]))}
b <- rep(0,(unique(sc$N)*2)-1)
for (i in 1:((unique(sc$N)*2)-1)){b[i] <- (indel[i]*log(joint_p[i])) + (snp[i]*log(joint_p[i])) + sum(indel[-i])*log(sum(joint_p[-i])) + sum(snp[-i])*log(sum(joint_p[-i]))}

# class_contrib <- rep(0,((unique(sc$N)*2)-1))
for (i in 1:((unique(sc$N)*2)-1)){class_contrib[i] <- 2*(a[i]-b[i])}

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)
cat('Chi-square test statistic =', chi_app, '\n')
# qchisq(p = 0.05, df = 260)
ndf <- ((unique(sc$N)*2)-1-1)*(length(tv)-1)
cat('Chi-square critical value at 0.01 with', ndf, 'df =', qchisq(p = 0.01, df = ndf), '\n')
cat('Test statistic - critical value =', chi_app-qchisq(p = 0.01, df = ndf), '\n')

chi_tb <- data.frame(chi2_stat=chi_app,
                     chi2_crit=qchisq(p = 0.01, df = ndf),
                     diff=chi_app-qchisq(p = 0.01, df = ndf))
chi_tb

# plot(1:((unique(sc$N)*2)-1), class_contrib)

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
  

da_cc

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]

ggsave(filename = paste('figures/SFS', parts[2], paste(unique(dtn$ECOT), collapse = '_'), tv[1], tv[2],
                        'chi_noinv.pdf', sep = "_"), plot = da_cc,
       width = 10, height = 7)