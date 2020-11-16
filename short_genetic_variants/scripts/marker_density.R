rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "reshape2", "tidyr", "tools", "data.table", "RColorBrewer", "dplyr", "textshape",
              "devtools", "optparse", "boot", "pscl", "lmtest", "gridExtra", "MASS")
# source("https://bioconductor.org/biocLite.R")
# biocLite("snpStats")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

option_list = list(
  make_option(c("-i", "--island"), type = "character", default = NULL,
              help = "Name of the island, CZA, CZB or CZD.", metavar = "character"),
  make_option(c("-e", "--ecotype"), type = "character", default = NULL,
              help = "Name of the ecotype, CRAB, WAVE_LEFT or WAVE_RIGHT.", metavar = "character"),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "Add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = "Rscript scripts/marker_density.R -i CZA -e CRAB --rminv",
                          option_list=option_list,
                          description = "Examine INDEL and SNP marker density per contig.")
opt = parse_args(opt_parser)

if (is.null(opt$island) | is.null(opt$ecotype)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

isl <- opt$island
eco <- opt$ecotype
# isl <- 'CZD'
# eco <- 'CRAB'
  
dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv")
dt$SIZE <- abs(nchar(as.character(dt$REF)) - nchar(as.character(dt$ALT)))
dt <- dt[dt$SIZE <= 50, ]
# range(dt$SIZE)
# rm_inv = FALSE

# head(dt)
if (opt$rminv) {
  cat('Variants inside inversions are being removed ...\n')
  dt <- dt[dt$invRui==FALSE, ]
}

# fxd <- dt[dt$DAF==0 | dt$DAF==1, ]
# fxd <- unique(fxd[, c('cp', 'VTYPE')])
# table(fxd$VTYPE)
# dtp <- dt[dt$DAF!=0 & dt$DAF!=1, ]
# nrow(dtp)
# dtp <- unique(dtp[, c('cp', 'VTYPE')])
# table(dtp$VTYPE)
# (tb <- data.frame(table(dtp$VTYPE), table(fxd$VTYPE)))
# tb$Freq/tb$Freq.1

dt <- dt[!is.na(dt$av) & dt$DAF!=0 & dt$DAF!=1, ]
# dt$ANN <- ifelse(test = is.na(dt$ANN), yes = 'coding', no = as.character(dt$ANN))
# dtu <- unique(dt[, c('cp', 'VTYPE', 'LG', 'av', 'invRui')])
# table(dtu$VTYPE, dtu$invRui)
# table(dtu$VTYPE)
# table(dtu$VTYPE, dtu$LG)
# vLG <- data.frame(table(dtu$VTYPE, dtu$LG))
# sum(vLG[vLG$Var1=='INDEL', 'Freq'])
# sum(vLG[vLG$Var1=='SNP', 'Freq'])
# table(dtu[is.na(dtu$LG), 'VTYPE'])
# nrow(dtu) - sum(data.frame(table(dtu$VTYPE, dtu$LG))$Freq)


# str(table(dt$ZONE, dt$VTYPE, dt$ECOT))
# table(dt$ZONE, dt$VTYPE, dt$ECOT)
# tot_v <- data.frame(rbind(table(dt$ZONE, dt$VTYPE, dt$ECOT)[1,,],
#                           table(dt$ZONE, dt$VTYPE, dt$ECOT)[2,,],
#                           table(dt$ZONE, dt$VTYPE, dt$ECOT)[3,,]),
#                     VTYPE = rep(c('INDEL', 'SNP'), 3),
#                     ZONE = c('CZA', 'CZA', 'CZB', 'CZB', 'CZD', 'CZD'))

dt$class <- paste(dt$ZONE, dt$ECOT, dt$VTYPE, sep = ":")
tot_v <- data.frame(table(dt$class))

# head(dt)
dt$av <- round(dt$av)
dt$LGAV <- paste(dt$LG, dt$av, sep = '_')

pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(LGAV = dt1$LGAV), length))
  # dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
# head(pr)
# pr[pr$Len==329,]

vt <- split(pr, f = pr$VTYPE)
# lapply(vt, head)
# identical(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$SNP$CHROM, vt$INDEL$CHROM)

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'LGAV'), all = TRUE)
# hist(mer$Len)
# range(mer$Len)
# head(mer)

# fai_path <- "data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
# fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# colnames(fai) <- c("CHROM", "Length")
# head(fai)
# fai <- fai[fai$Length >= 250, ]

# mer <- merge(mer, fai, by = 'CHROM')

# mer$pal <- cut(mer$Length, breaks = 10, labels = FALSE)
# table(mer$pal)
# head(mer)
# mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
# mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

mer <- separate(data = mer, col = LGAV, into = c('LG', 'AV'), sep = "_")
mer$RATIO <- mer$x.y / mer$x.x
# table(mer$ECOT)
mer$ECOT <- factor(mer$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))

# len_pal <- colorRampPalette(c("grey", "black"))
# pprop <- ggplot(data = mer, aes(x = prop.x, y = prop.y)) +
#   facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
#   geom_abline(slope = 1, linetype = 'dashed') +
#   # geom_point(aes(col = as.factor(pal))) +
#   geom_point(alpha = 0.4) +
#   # geom_smooth(method='lm', formula= y~x) +
#   # scale_color_manual(values = len_pal(10)) +
#   labs(x = 'INDEL relative proportion', y = 'SNP relative proportion', col = 'bin 50000 bp') +
#   theme(legend.text = element_text(size = 12), legend.position = 'top',
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 12),
#         # legend.position = "top",
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "#91bfdb", color = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2)) +
#   guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
# pprop
# mer[which.max(mer$prop.x), ]
# if (opt$rminv) {
#   ggsave(filename = "figures/MD_snp_vs_indel_prop_noinv.pdf", plot = pprop, width = 10, height = 10)
# } else {
#   ggsave(filename = "figures/MD_snp_vs_indel_prop.pdf", plot = pprop, width = 10, height = 10)
# }

# head(mer)
# mer$IE <- paste(mer$ISL, mer$ECOT, sep = ':')
# get_vratio <- function(data, cm, fac) {
#   dd <- data[data[, cm] == fac, ]
#   rr <- sum(dd$x.x)/sum(dd$x.y)
#   # data[data[, cm] == fac, 'VRATIO'] <- rr
#   return(rr)
# }
# # mer$VRATIO <- 0
# for (i in 1:9) {
#   print(get_vratio(data = mer, cm = 'IE', fac = unique(mer$IE)[i]))
# }

# REMOVE SHORT CONTIGS
# cf <- 500
# mer <- mer[mer$Length >= cf, ]

# REMOVE SNP WITH 0 COUNT
# mer <- mer[mer$x.y!=0, ]
# sum(mer$x.x)/sum(mer$x.y)

# VARIANT COUNT + 1
# mer$x.x <- mer$x.x + 1
# mer$x.y <- mer$x.y + 1

# pcount <- ggplot(data = mer, aes(x = x.y, y = x.x)) +
#   facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
#   # geom_abline(slope = 0.2, linetype = 'dashed') +
#   geom_abline(slope = 0.2, linetype = 'dashed', size = 1.2) +
#   geom_point(alpha = 0.2) +
#   # geom_abline(slope = 0.06, intercept = -0.16, col = 'blue') +
#   # geom_point(aes(col = as.factor(pal))) +
#   geom_smooth(method = 'lm', formula = y~x) +
#   # scale_color_manual(values = len_pal(10)) +
#   labs(x = 'SNP count', y = 'INDEL count', col = 'bin 50000 bp') +
#   theme(legend.text = element_text(size = 12), legend.position = 'top',
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 12),
#         # legend.position = "top",
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "#91bfdb", color = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2)) +
#   guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
# pcount
# ggsave(filename = "figures/MD_indel_vs_snp_count_lm_noinv.pdf", plot = pcount, width = 10, height = 10, dpi = "screen")
# summary(lm(formula = x.x~x.y, data = mer))

# if (opt$rminv) {
#   ggsave(filename = "figures/MD_indel_vs_snp_count_noinv.pdf", plot = pcount, width = 10, height = 10)
# } else {
#   ggsave(filename = "figures/MD_indel_vs_snp_count.pdf", plot = pcount, width = 10, height = 10)
# }
# 
# 
# 
## one island and one ecotype at the time
tt <- as.data.frame(mer[mer$ISL==isl & mer$ECOT==eco, ])
tt$ZONE <- tt$ISL
tt$ISL <- NULL
# head(tt)
# str(tt)
tt$av <- as.integer(tt$AV)
tt$LG <- factor(tt$LG, levels = 1:17)
tt <- merge(x = tt, y = unique(dt[, c(intersect(colnames(tt), colnames(dt)), 'invRui')]), all.x = TRUE)
# nrow(unique(tt))

tt$pal <- ifelse(test = is.infinite(tt$RATIO), yes = 'red', no = 'black')
tt$pal <- ifelse(test = tt$invRui==FALSE, yes = tt$pal, no = 'green')

# pratio <- ggplot(data = tt, aes(x = av, y = RATIO)) +
#   facet_wrap(~ LG) +
#   geom_point(alpha = 0.5, size = 2, col = tt$pal) +
#   labs(x = 'Map position', y = 'SNP to INDEL') +
#   theme(legend.text = element_text(size = 12), legend.position = 'top',
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 12),
#         # legend.position = "top",
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "#91bfdb", color = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2)) +
#   guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
# pratio
# ggsave(filename = "figures/MD_mappos_snp_to_indel_count.pdf", plot = pratio, width = 10, height = 10)

sl <- sum(tt$x.x)/sum(tt$x.y)
pcount <- ggplot(data = tt, aes(x = x.y, y = x.x)) +
  facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_abline(slope = sl, linetype = 'dashed', size = 1.2) +
  geom_point(alpha = 0.2, size = 2, col = tt$pal) +
  geom_smooth(method = 'lm', formula = y~x) +
  labs(x = 'SNP count', y = 'INDEL count') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
ggsave(filename = paste('figures/MD_lm', isl, eco, 'noinv_indel_snp_count.pdf', sep = "_"),
       plot = pcount, scale = 0.8, dpi = "screen")

lms <- summary(lm(formula = x.x ~ x.y, data = tt[tt$invRui==FALSE, ]))
lmsd <- round(data.frame(lms$coefficients), 3)
lmsd$Par <- row.names(lmsd)
lmsd <- lmsd[, c(ncol(lmsd), 1:(ncol(lmsd)-1))]
lmsd <- rbind(lmsd, data.frame(Par='Expected slope', Estimate=round(sl, 3), Std..Error=NA, t.value=NA, Pr...t..=NA))
lmsd$ISL <- isl
lmsd$ECOT <- eco
write.table(x = lmsd, file = paste('results/marker_density/lm/MD_lm', isl, eco, 'noinv_indel_snp_count.txt', sep = "_"),
            append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
# tt <- mer
# tt <- mer[mer$ISL=='CZA' & mer$ECOT=='WAVE_LEFT', ]
# head(tt)
# tt[tt$x.x==0 & tt$x.y==0,]
# tt[tt$prop.x==0 & tt$prop.y==0,]
# table(tt$x.y)
# REMOVE SNP WITH 0 COUNT
# tt <- tt[tt$x.y!=0, ]
# write.table(x = tt, file = 'test/marker_density/CZA_CRAB_marker_density.csv', quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
# tt <- read.csv(file = 'CZA_CRAB_marker_density.csv')

## PERMUTATION LINEAR MODEL
bar <- replicate(n = 10000, expr = rbinom(n = nrow(tt), size = tt$x.y, prob = sum(tt$x.x)/sum(tt$x.y)))
# bar[1:10,1:10]
# dim(bar)
mexp <- apply(X = bar, MARGIN = 2, FUN = function(x) {
  mlin <- lm(x ~ tt$x.y)
  coef(mlin)
  # mpoi <- glm(x ~ tt$x.y, family = poisson(link = 'log'))
  # exp(coef(mpoi))
})
cexp <- as.data.frame(t(mexp))
colnames(cexp) <- c('intercept', 'slope')
# head(cexp)
cexp$cpal <- 'bootstrap'
cexp$ZONE <- isl
cexp$ECOT <- eco
# with(data = cexp, plot(x = intercept, slope))
# tt$SUM.x <- sum(tt$x.x)
# tt.poi <- glm(x.x ~ x.y + offset(log(tt$SUM.x)), data = tt, family=poisson(link = 'log'))
# tt.poi <- glm(x.x ~ x.y + offset(log(tt$Length)), data = tt, family=poisson(link = 'log'))

# M1 <- glm(x.x ~ x.y, data = tt, family = poisson(link = 'log'))
M1 <- lm(x.x ~ x.y, data = tt)
# cm1 <- as.data.frame(t(exp(coef(M1))))
cm1 <- as.data.frame(t(coef(M1)))
colnames(cm1) <- c('intercept', 'slope')
cm1$cpal <- 'observation'

# ttf <- rbind(cm1, cexp)
# ttf$cpal <- factor(x = ttf$cpal, levels = c('observation', 'bootstrap'))
g <- ggplot(data = cexp, aes(x = slope)) +
  facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_histogram(bins = 20, fill = 'black', col = 'white') +
  geom_vline(xintercept = cm1$slope, col = 'green', size = 1.5) +
  # labs(col = '', title = paste(isl, eco)) +
  # scale_color_manual(values = c('black', 'blue')) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
        # legend.text = element_text(size = 18),
        # title = element_text(size = 16)) +
  # guides(col = guide_legend(override.aes = list(size=2)))
# g

ggsave(filename = paste('figures/MD_boot', isl, eco, 'noinv_slope.pdf', sep = "_"),
       plot = g, scale = 0.8, dpi = "screen")

cat('\n...MISSION COMPLETE!...\n')
