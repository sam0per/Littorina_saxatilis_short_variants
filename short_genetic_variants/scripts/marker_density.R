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
# isl <- 'CZA'
# eco <- 'WAVE_LEFT'
  
dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv")

# rm_inv = FALSE

# head(dt)
if (opt$rminv) {
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
pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(CHROM = dt1$CHROM), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
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

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'CHROM'), all = TRUE)
# hist(mer$Len)
# range(mer$Len)
# head(mer)

fai_path <- "data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
colnames(fai) <- c("CHROM", "Length")
# head(fai)

mer <- merge(mer, fai, by = 'CHROM')

mer$pal <- cut(mer$Length, breaks = 10, labels = FALSE)
# table(mer$pal)
# head(mer)
mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

# table(mer$ECOT)
mer$ECOT <- factor(mer$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))

len_pal <- colorRampPalette(c("grey", "black"))
pprop <- ggplot(data = mer, aes(x = prop.x, y = prop.y)) +
  facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
  geom_abline(slope = 1, linetype = 'dashed') +
  # geom_point(aes(col = as.factor(pal))) +
  geom_point(alpha = 0.4) +
  # geom_smooth(method='lm', formula= y~x) +
  # scale_color_manual(values = len_pal(10)) +
  labs(x = 'INDEL relative proportion', y = 'SNP relative proportion', col = 'bin 50000 bp') +
  theme(legend.text = element_text(size = 12), legend.position = 'top',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        # legend.position = "top",
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
# pprop
# mer[which.max(mer$prop.x), ]
if (opt$rminv) {
  ggsave(filename = "figures/MD_snp_vs_indel_prop_noinv.pdf", plot = pprop, width = 10, height = 10)
} else {
  ggsave(filename = "figures/MD_snp_vs_indel_prop.pdf", plot = pprop, width = 10, height = 10)
}

# head(mer)
mer$IE <- paste(mer$ISL, mer$ECOT, sep = ':')
get_vratio <- function(data, cm, fac) {
  dd <- data[data[, cm] == fac, ]
  rr <- sum(dd$x.x)/sum(dd$x.y)
  # data[data[, cm] == fac, 'VRATIO'] <- rr
  return(rr)
}
# mer$VRATIO <- 0
for (i in 1:9) {
  print(get_vratio(data = mer, cm = 'IE', fac = unique(mer$IE)[i]))
}

# REMOVE SNP WITH 0 COUNT
# mer <- mer[mer$x.y!=0, ]
# sum(mer$x.x)/sum(mer$x.y)

# VARIANT COUNT + 1
# mer$x.x <- mer$x.x + 1
# mer$x.y <- mer$x.y + 1

pcount <- ggplot(data = mer, aes(x = x.y, y = x.x)) +
  facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
  # geom_abline(slope = 0.2, linetype = 'dashed') +
  geom_abline(slope = 0.2, linetype = 'dashed', size = 1.2) +
  geom_point(alpha = 0.2) +
  # geom_abline(slope = 0.06, intercept = -0.16, col = 'blue') +
  # geom_point(aes(col = as.factor(pal))) +
  geom_smooth(method = 'lm', formula = y~x) +
  # scale_color_manual(values = len_pal(10)) +
  labs(x = 'SNP count', y = 'INDEL count', col = 'bin 50000 bp') +
  theme(legend.text = element_text(size = 12), legend.position = 'top',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        # legend.position = "top",
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  guides(col = guide_legend(override.aes = list(size=3), nrow = 1))
# pcount
# ggsave(filename = "figures/MD_indel_vs_snp_count_lm_noinv.pdf", plot = pcount, width = 10, height = 10, dpi = "screen")
# summary(lm(formula = x.x~x.y, data = mer))

if (opt$rminv) {
  ggsave(filename = "figures/MD_indel_vs_snp_count_noinv.pdf", plot = pcount, width = 10, height = 10)
} else {
  ggsave(filename = "figures/MD_indel_vs_snp_count.pdf", plot = pcount, width = 10, height = 10)
}
# 
# 
# 
## one island and one ecotype at the time
tt <- mer[mer$ISL==isl & mer$ECOT==eco, ]
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

## PERMUTATION POISSON MODEL
bar <- replicate(n = 10000, expr = rbinom(n = nrow(tt), size = tt$x.y, prob = sum(tt$x.x)/sum(tt$x.y)))
# head(bar)
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

ttf <- rbind(cm1, cexp)

g <- ggplot(data = ttf, aes(x = intercept, y = slope, col = cpal)) +
  geom_point(size = 2) +
  labs(col = '', title = paste(isl, eco)) +
  scale_color_manual(values = c('black', 'blue')) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 18),
        title = element_text(size = 16)) +
  guides(col = guide_legend(override.aes = list(size=2)))
# g

ggsave(filename = paste('figures/MD_boot', isl, eco, 'all_slope_intercept.pdf', sep = "_"),
       plot = g, width = 10, height = 7, dpi = "screen")

rm(list = setdiff(x = ls(), y = c('cm1', 'cexp', 'isl', 'eco')))

fs <- list.files(path = 'results/marker_density', pattern = paste(isl, eco, sep = '_'), full.names = TRUE)
fs <- fs[grepl(pattern = 'nongenic', x = fs)]
ind <- read.table(file = fs[1], header = TRUE, sep = '\t')
snp <- read.table(file = fs[2], header = TRUE, sep = '\t')
cs <- intersect(colnames(ind), colnames(snp))
ind <- ind[, colnames(ind) %in% cs]
snp <- snp[, colnames(snp) %in% cs]
dt <- as.data.frame(rbind(ind, snp))
# dt <- merge(x = as.data.frame(rbind(ind, snp)), y = dt, all.y = TRUE)

dt <- dt[!is.na(dt$av) & dt$DAF!=0 & dt$DAF!=1, ]

dt$class <- paste(dt$ZONE, dt$ECOT, dt$VTYPE, sep = ":")
tot_v <- data.frame(table(dt$class))
tot_v[1,2]/tot_v[2,2]

# head(dt)
pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(CHROM = dt1$CHROM), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
# head(pr)

vt <- split(pr, f = pr$VTYPE)
# lapply(vt, head)
# identical(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$SNP$CHROM, vt$INDEL$CHROM)

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'CHROM'), all = TRUE)
# head(mer)
# sum(is.na(mer$x.x))
# sum(is.na(mer$x.y))

mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

# table(mer$ECOT)
# sum(mer$x.x)/sum(mer$x.y)

# summary(lm(formula = x.x~x.y, data = mer))

# VARIANT COUNT + 1
# mer$x.x <- mer$x.x + 1
# mer$x.y <- mer$x.y + 1

M2 <- lm(x.x ~ x.y, data = mer)
cm2 <- as.data.frame(t(coef(M2)))
colnames(cm2) <- c('intercept', 'slope')
cm2$cpal <- 'nongenic'

## PERMUTATION POISSON MODEL
bar <- replicate(n = 10000, expr = rbinom(n = nrow(mer), size = mer$x.y, prob = sum(mer$x.x)/sum(mer$x.y)))
# head(bar)
mexp <- apply(X = bar, MARGIN = 2, FUN = function(x) {
  mlin <- lm(x ~ mer$x.y)
  coef(mlin)
  # mpoi <- glm(x ~ tt$x.y, family = poisson(link = 'log'))
  # exp(coef(mpoi))
})
cexp <- as.data.frame(t(mexp))
colnames(cexp) <- c('intercept', 'slope')
# head(cexp)
cexp$cpal <- 'bootstrap nongenic'

# ttf <- rbind(cm1, cm2, cexp)
ttf <- rbind(cm2, cexp)

g <- ggplot(data = ttf, aes(x = intercept, y = slope, col = cpal)) +
  geom_point(size = 2) +
  labs(col = '', title = paste(isl, eco)) +
  scale_color_manual(values = c('black', 'red')) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        title = element_text(size = 16)) +
  guides(col = guide_legend(override.aes = list(size=2)))
# g

# ggsave(filename = paste('figures/marker_dens_perm', opt$island, opt$ecotype, 'slope_intercept.pdf', sep = "_"),
#        plot = g, width = 10, height = 7)

rm(list = setdiff(x = ls(), y = c('ttf', 'isl', 'eco')))

fs <- list.files(path = 'results/marker_density', pattern = paste(isl, eco, sep = '_'), full.names = TRUE)
fs <- fs[grepl(pattern = 'nonsyn', x = fs)]
ind <- read.table(file = fs[1], header = TRUE, sep = '\t')
snp <- read.table(file = fs[2], header = TRUE, sep = '\t')
cs <- intersect(colnames(ind), colnames(snp))
ind <- ind[, colnames(ind) %in% cs]
snp <- snp[, colnames(snp) %in% cs]
dt <- as.data.frame(rbind(ind, snp))
# dt <- merge(x = as.data.frame(rbind(ind, snp)), y = dt, all.y = TRUE)

dt <- dt[!is.na(dt$av) & dt$DAF!=0 & dt$DAF!=1, ]

dt$class <- paste(dt$ZONE, dt$ECOT, dt$VTYPE, sep = ":")
tot_v <- data.frame(table(dt$class))
tot_v[1,2]/tot_v[2,2]

# head(dt)
pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(CHROM = dt1$CHROM), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
# head(pr)

vt <- split(pr, f = pr$VTYPE)
# lapply(vt, head)
# identical(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$SNP$CHROM, vt$INDEL$CHROM)

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'CHROM'), all = TRUE)
# head(mer)
# sum(is.na(mer$x.x))
# sum(is.na(mer$x.y))

mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

# table(mer$ECOT)
# sum(mer$x.x)/sum(mer$x.y)

# summary(lm(formula = x.x~x.y, data = mer))

# VARIANT COUNT + 1
# mer$x.x <- mer$x.x + 1
# mer$x.y <- mer$x.y + 1

M3 <- lm(x.x ~ x.y, data = mer)
cm3 <- as.data.frame(t(coef(M3)))
colnames(cm3) <- c('intercept', 'slope')
cm3$cpal <- 'nonsyn'
ttf <- rbind(cm3, ttf)
table(ttf$cpal)

## PERMUTATION POISSON MODEL
bar <- replicate(n = 10000, expr = rbinom(n = nrow(mer), size = mer$x.y, prob = sum(mer$x.x)/sum(mer$x.y)))
# head(bar)
mexp <- apply(X = bar, MARGIN = 2, FUN = function(x) {
  mlin <- lm(x ~ mer$x.y)
  coef(mlin)
  # mpoi <- glm(x ~ tt$x.y, family = poisson(link = 'log'))
  # exp(coef(mpoi))
})
cexp <- as.data.frame(t(mexp))
colnames(cexp) <- c('intercept', 'slope')
# head(cexp)
cexp$cpal <- 'bootstrap nonsyn'
ttf <- rbind(cexp, ttf)
table(ttf$cpal)
ttf$cpal <- factor(ttf$cpal, levels = c('bootstrap nongenic', 'nongenic', 'bootstrap nonsyn', 'nonsyn'))

f <- ggplot(data = ttf, aes(x = intercept, y = slope, col = cpal)) +
  geom_point(size = 2) +
  labs(col = '', title = paste(isl, eco)) +
  scale_color_manual(values = brewer.pal(n = length(levels(ttf$cpal)), name = 'Paired')) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        title = element_text(size = 16)) +
  guides(col = guide_legend(override.aes = list(size=2)))
f

rm(list = setdiff(x = ls(), y = c('ttf', 'isl', 'eco')))

fs <- list.files(path = 'results/marker_density', pattern = paste(isl, eco, sep = '_'), full.names = TRUE)
fs <- fs[grepl(pattern = '_syn', x = fs)]
ind <- read.table(file = fs[1], header = TRUE, sep = '\t')
snp <- read.table(file = fs[2], header = TRUE, sep = '\t')
cs <- intersect(colnames(ind), colnames(snp))
ind <- ind[, colnames(ind) %in% cs]
snp <- snp[, colnames(snp) %in% cs]
dt <- as.data.frame(rbind(ind, snp))
# dt <- merge(x = as.data.frame(rbind(ind, snp)), y = dt, all.y = TRUE)

dt <- dt[!is.na(dt$av) & dt$DAF!=0 & dt$DAF!=1, ]

dt$class <- paste(dt$ZONE, dt$ECOT, dt$VTYPE, sep = ":")
tot_v <- data.frame(table(dt$class))
tot_v[1,2]/tot_v[2,2]

# head(dt)
pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(CHROM = dt1$CHROM), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
# head(pr)

vt <- split(pr, f = pr$VTYPE)
# lapply(vt, head)
# identical(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$INDEL$CHROM, vt$SNP$CHROM)
# setdiff(vt$SNP$CHROM, vt$INDEL$CHROM)

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'CHROM'), all = TRUE)
# head(mer)
# sum(is.na(mer$x.x))
# sum(is.na(mer$x.y))

mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

# table(mer$ECOT)
# sum(mer$x.x)/sum(mer$x.y)

# summary(lm(formula = x.x~x.y, data = mer))

# VARIANT COUNT + 1
# mer$x.x <- mer$x.x + 1
# mer$x.y <- mer$x.y + 1

M4 <- lm(x.x ~ x.y, data = mer)
cm4 <- as.data.frame(t(coef(M4)))
colnames(cm4) <- c('intercept', 'slope')
cm4$cpal <- 'syn'
ttf <- rbind(cm4, ttf)

## PERMUTATION POISSON MODEL
bar <- replicate(n = 10000, expr = rbinom(n = nrow(mer), size = mer$x.y, prob = sum(mer$x.x)/sum(mer$x.y)))
# head(bar)
mexp <- apply(X = bar, MARGIN = 2, FUN = function(x) {
  mlin <- lm(x ~ mer$x.y)
  coef(mlin)
  # mpoi <- glm(x ~ tt$x.y, family = poisson(link = 'log'))
  # exp(coef(mpoi))
})
cexp <- as.data.frame(t(mexp))
colnames(cexp) <- c('intercept', 'slope')
# head(cexp)
cexp$cpal <- 'bootstrap syn'
ttf <- rbind(cexp, ttf)
table(ttf$cpal)
ttf$cpal <- factor(ttf$cpal, levels = c('bootstrap nongenic', 'nongenic', 'bootstrap nonsyn', 'nonsyn',
                                        'bootstrap syn', 'syn'))

f <- ggplot(data = ttf, aes(x = intercept, y = slope, col = cpal)) +
  geom_point(size = 2) +
  labs(col = '', title = paste(isl, eco)) +
  scale_color_manual(values = brewer.pal(n = length(levels(ttf$cpal)), name = 'Paired')) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        title = element_text(size = 16)) +
  guides(col = guide_legend(override.aes = list(size=2)))
f

ggsave(filename = paste('figures/MD_boot', isl, eco, 'ann_slope_intercept.pdf', sep = "_"),
       plot = f, width = 10, height = 7, dpi = "screen")

# est1 <- cbind(Estimate = coef(M1), confint(M1))

## Check for over/underdispersion in the model
# E2 <- resid(M1, type = "pearson")
# N <- nrow(tt)
# p <- length(coef(M1))   
# cat(as.character(M1$call), ': Dispersion =', sum(E2^2) / (N - p), '\n')

## NEGATIVE BINOMIAL
# M2 <- glm.nb(x.x ~ x.y, data = tt)
# summary(M2)

## Check for over/underdispersion in the model
# E2 <- resid(M2, type = "pearson")
# N <- nrow(tt)
# p <- length(coef(M2))   
# cat(as.character(M2$call), ': Dispersion =', sum(E2^2) / (N - p), '\n')

# pchisq(2 * (logLik(M2) - logLik(M1)), df = 1, lower.tail = FALSE)
# est <- cbind(Estimate = coef(M2), confint(M2))
# exp(est)
# newdata2 <- data.frame(x.y = seq(from = min(tt$x.y), to = max(tt$x.y), length.out = 1000))
# 
# newdata2 <- cbind(newdata2,
#                   M1=predict(M1, newdata2, type = "link", se.fit=TRUE),
#                   M2=predict(M2, newdata2, type = "link", se.fit=TRUE))
# head(newdata2)
# newdata2 <- within(newdata2, {
#   INDEL.M1 <- exp(M1.fit)
#   LL.M1 <- exp(M1.fit - 1.96 * M1.se.fit)
#   UL.M1 <- exp(M1.fit + 1.96 * M1.se.fit)
#   INDEL.M2 <- exp(M2.fit)
#   LL.M2 <- exp(M2.fit - 1.96 * M2.se.fit)
#   UL.M2 <- exp(M2.fit + 1.96 * M2.se.fit)
# })

# ggplot(newdata2, aes(x.y, INDELc)) +
#   geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
#   geom_line()

## ZERO-INFLATED MODELS
# library(pscl)
# cat('Percentage of 0 in the y variable =', 100*sum(tt$x.x == 0)/nrow(tt), '\n')

# 100*sum(tt$x.y == 0)/nrow(tt)
# sum(tt$x.y == 0)
# sum(tt$x.x == 0)
# sum(tt$prop.y == 1)
# sum(tt$prop.y == 0)
# sum(tt$prop.x == 1)
# sum(tt$prop.x == 0)

# model.zi = zeroinfl(x.x ~ x.y,
#                     data = tt,
#                     dist = "poisson")
### dist = "negbin" may be used
# summary(model.zi)
# model.zi$coefficients
# fitted(model.zi)[1]
# tt$x.y[1]
## Dispersion statistic
# E2 <- resid(model.zi, type = "pearson")
# N  <- nrow(tt)
# p  <- length(coef(model.zi))
# cat(as.character(model.zi$call), ': Dispersion =', sum(E2^2) / (N - p), '\n')
# 
# M4 <- zeroinfl(x.x ~ x.y, data = tt, dist = "negbin")
# summary(M4)
## Dispersion Statistic
# E2 <- resid(M4, type = "pearson")
# N  <- nrow(tt)
# p  <- length(coef(M4)) + 1 # '+1' is due to theta
# cat(as.character(M4$call), ': Dispersion =', sum(E2^2) / (N - p), '\n')
# 
# lrtest(model.zi, M4)
# vuong(M1,M4)
# vuong(M2,M4)
# vuong(M2,M1)

# glm(x.x ~ prop.y + offset(log(tt$SUM.x)), data = tt, family=poisson(link = 'log'))

# Qx <- quantile(tt$x.x, probs=c(.25, .75), na.rm = FALSE)
# iqrx <- IQR(tt$x.x)
# upx <- Qx[2]+1.5*iqrx # Upper Range  
# lowx <- Qx[1]-1.5*iqrx # Lower Range
# outlx <- subset(tt, tt$x.x < (Qx[1] - 1.5*iqrx) | tt$x.x > (Qx[2]+1.5*iqrx))
# 
# Qy <- quantile(tt$x.y, probs=c(.25, .75), na.rm = FALSE)
# iqry <- IQR(tt$x.y)
# upy <- Qy[2]+1.5*iqry # Upper Range  
# lowy <- Qy[1]-1.5*iqry # Lower Range
# outly <- subset(tt, tt$x.y < (Qy[1] - 1.5*iqry) | tt$x.y > (Qy[2]+1.5*iqry))
# 
# tt_outl <- unique(rbind(outlx, outly))
# cat('Count of outliers per contig length:\n')
# data.frame(table(tt_outl$pal))

# tt1_outl <- tt_outl[, c('CHROM', 'Length')]
# if (file.exists('results/Lsax_short_var_czs_outliers.csv')) {
#   outl <- read.csv('results/Lsax_short_var_czs_outliers.csv')
#   tt1_outl <- merge(tt1_outl, outl)
# }
# 
# write.table(x = tt1_outl, file = 'results/Lsax_short_var_czs_outliers.csv', quote = FALSE,
#             sep = ',', row.names = FALSE)

# mytheme.r <- gridExtra::ttheme_default(
#   core = list(fg_params=list(cex = 1)),
#   colhead = list(fg_params=list(cex = 1, col = 'red')),
#   rowhead = list(fg_params=list(cex = 1)))
# mytheme.b <- gridExtra::ttheme_default(
#   core = list(fg_params=list(cex = 1)),
#   colhead = list(fg_params=list(cex = 1, col = 'blue')),
#   rowhead = list(fg_params=list(cex = 1)))
# g <- ggplot(data = tt, aes(x = x.y, y = x.x)) +
#   geom_point(col = 'black') +
#   # geom_text(data = tt_outl, aes(x = x.y, y = x.x, label = CHROM),hjust = 0, vjust = 0) +
#   geom_ribbon(data = newdata2, aes(x = x.y, y = INDEL.M1, ymin = LL.M1, ymax = UL.M1), fill = 'blue', alpha = .25) +
#   geom_ribbon(data = newdata2, aes(x = x.y, y = INDEL.M2, ymin = LL.M2, ymax = UL.M2), fill = 'red', alpha = .25) +
#   geom_line(data = newdata2, aes(x = x.y, y = INDEL.M1), col = 'blue') +
#   geom_line(data = newdata2, aes(x = x.y, y = INDEL.M2), col = 'red') +
#   annotation_custom(tableGrob(round(exp(est), 3), theme = mytheme.r),
#                     xmin=5, xmax=25, ymin = max(newdata2$INDEL.M2)-10, ymax=max(newdata2$INDEL.M2)) +
#   annotation_custom(tableGrob(round(exp(est1), 3), theme = mytheme.b),
#                     xmin=30, xmax=60, ymin = max(newdata2$INDEL.M2)-10, ymax=max(newdata2$INDEL.M2)) +
#   labs(x = 'SNP count', y = 'INDEL count', title = paste(opt$island, opt$ecotype))
# ggsave(filename = paste('figures/model', opt$island, opt$ecotype, 'SNP_INDEL_count.pdf', sep = "_"),
#        plot = g, width = 8, height = 7)
# g <- ggplot(data = tt, aes(x = prop.x, y = prop.y)) + geom_point(col = 'black')
# geom_abline(slope = M4$coefficients$count[2], intercept = M4$coefficients$count[1], col = 'red') +
# geom_abline(slope = tt.poi$coefficients[2], intercept = tt.poi$coefficients[1], col = 'red') +
# geom_abline(slope = 1, linetype = 'dashed')
# geom_smooth(method = 'lm', formula = y~x)

# library(lmtest)
# library(sandwich)

# tt.ols <- lm(prop.x ~ prop.y, data = tt)
# tt.ols <- lm(x.x ~ x.y, data = tt)
# round(summary(tt.ols)$coefficients, 3)
# ggplot(data = tt, aes(x = x.y, y = x.x)) +
#   geom_point() +
#   geom_smooth(method = 'lm', formula = y~x) +
#   geom_abline(slope = 1, linetype = 'dashed')
# tt.ols <- lm(prop.y ~ prop.x, data = tt)
# tt.ols <- lm(x.y ~ x.x, data = tt)
# summary(tt.ols)
# round(summary(tt.ols)$coefficients, 3)
# ggplot(data = tt, aes(x = x.x, y = x.y)) +
#   geom_point() +
#   geom_smooth(method = 'lm', formula = y~x) +
#   geom_abline(slope = 1, linetype = 'dashed')
# 
# ggplot(data = tt, aes(x = prop.y, y = prop.x)) +
#   geom_point() +
#   geom_smooth(method = 'lm', formula = y~x) +
#   geom_abline(slope = 1, linetype = 'dashed')
# ggplot(data = tt, aes(x = prop.x, y = prop.y)) +
#   geom_point() +
#   geom_smooth(method = 'lm', formula = y~x) +
#   geom_abline(slope = 1, linetype = 'dashed')
# 
# summary(tt.ols)$coefficients[2, 1]
# summary(tt.ols)$coefficients[2, 2]
# confint(tt.ols)[2,1]
# confint(tt.ols)[2,2]
# tt.white <- coeftest(tt.ols, vcov = vcovHC(tt.ols, "HC1"))   # HC1 gives us the White standard errors
# tt.white[2,1]
# tt.white[2,2]
# confint(tt.white)[2,1]
# confint(tt.white)[2,2]
# 
# 
# library(caret)
# # Define training control
# train.control <- trainControl(method = "repeatedcv", 
#                               number = 10, repeats = 5)
# # Train the model
# model <- train(prop.x ~ prop.y, data = tt, method = "lm",
#                trControl = train.control)
# # model <- train(prop.y ~ prop.x, data = tt, method = "lm",
# #                trControl = train.control)
# # Summarize the results
# print(model)
# summary(model)
# str(model)
# summary(model$finalModel)
# identical(fitted(model), fitted(model$finalModel))
# model$call
# str(model$finalModel)
# tt$FIT_REPCV <- fitted(model)
# g <- ggplot(data = tt, aes(x = prop.y, y = prop.x)) + geom_point(col = 'black')
# # g <- ggplot(data = tt, aes(x = prop.x, y = prop.y)) + geom_point(col = 'black')
# g + geom_abline(slope = model$finalModel$coefficients[2], intercept = model$finalModel$coefficients[1], col = 'red') +
#   geom_abline(slope = 1, linetype = 'dashed') +
#   geom_smooth(method = 'lm', formula = y~x)
# confint(model$finalModel)
# confint(object = lm(formula = prop.x~prop.y, data = tt))
# summary(lm(formula = prop.x~prop.y, data = tt))
# summary(model)
# 
# # install.packages('lmtest')
# library(lmtest)
# 
# tt.ols <- lm(prop.y ~ prop.x, data = tt)
# 
# # test heteroskedasticity
# bptest(tt.ols)
# 
# # Generalized Least Squares With Unknown Form of Variance
# tt$resi <- tt.ols$residuals
# tt$prop.x <- ifelse(test = tt$prop.x==0, yes = 0.00000001, no = tt$prop.x)
# varfunc.ols <- lm(log(resi^2) ~ log(prop.x), data = tt)
# tt$varfunc <- exp(varfunc.ols$fitted.values)
# tt.gls <- lm(prop.y ~ prop.x, weights = 1/sqrt(varfunc), data = tt)
# 
# tt.wlm <- lm(prop.y ~ prop.x, weights = Length, data = tt)
# 
# summary(tt.ols)
# summary(tt.gls)
# summary(tt.wlm)
# 
# g <- ggplot(data = tt, aes(y = prop.y, x = prop.x)) + geom_point(col = 'black')
# g + geom_abline(slope = tt.ols$coefficients[2], intercept = tt.ols$coefficients[1], col = 'red') + 
#   geom_abline(slope = tt.gls$coefficients[2], intercept = tt.gls$coefficients[1], col = 'green') +
#   geom_abline(slope = tt.wlm$coefficients[2], intercept = tt.wlm$coefficients[1], col = 'blue')
# 
# head(mer)
mer$class <- paste(mer$ISL, mer$ECOT, sep = ":")
# invisible(lapply(unique(mer$class), function(x) {
#   dt1 <- mer[mer$class==x, ]
#   # val <- lm(formula = prop.y ~ prop.x, data = dt1)
#   # val <- sum(dt1$prop.x, dt1$prop.y)
#   # cat(x, val, '\n')
#   order(dt1)
# }))

# 
# 
# 
## CORRELATION INDEL-SNP PROPORTIONS
# head(mer)
invisible(lapply(unique(mer$class), function(x) {
  dt1 <- mer[mer$class==x, ]
  val <- cor(dt1$prop.x, dt1$prop.y)
  # val <- sum(dt1$prop.x, dt1$prop.y)
  cat(x, round(val,3), '\n')
}))


