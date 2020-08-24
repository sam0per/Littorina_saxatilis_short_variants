rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "reshape2", "tidyr", "tools", "data.table", "RColorBrewer", "dplyr", "textshape", "plotly",
              "devtools", "BiocManager", "limma")
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

rm_inv = FALSE
dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv")
head(dt)
if (rm_inv) {
  dt <- dt[dt$invRui==FALSE, ]
}

fxd <- dt[dt$DAF==0 | dt$DAF==1, ]
fxd <- unique(fxd[, c('cp', 'VTYPE')])
table(fxd$VTYPE)
dtp <- dt[dt$DAF!=0 & dt$DAF!=1, ]
# nrow(dtp)
dtp <- unique(dtp[, c('cp', 'VTYPE')])
table(dtp$VTYPE)
(tb <- data.frame(table(dtp$VTYPE), table(fxd$VTYPE)))
tb$Freq/tb$Freq.1

dt <- dt[!is.na(dt$av) & dt$DAF!=0 & dt$DAF!=1, ]

dtu <- unique(dt[, c('cp', 'VTYPE', 'LG', 'av', 'invRui')])
table(dtu$VTYPE, dtu$invRui)
table(dtu$VTYPE)
table(dtu$VTYPE, dtu$LG)
vLG <- data.frame(table(dtu$VTYPE, dtu$LG))
sum(vLG[vLG$Var1=='INDEL', 'Freq'])
sum(vLG[vLG$Var1=='SNP', 'Freq'])
table(dtu[is.na(dtu$LG), 'VTYPE'])
# nrow(dtu) - sum(data.frame(table(dtu$VTYPE, dtu$LG))$Freq)


# str(table(dt$ZONE, dt$VTYPE, dt$ECOT))
# table(dt$ZONE, dt$VTYPE, dt$ECOT)
# tot_v <- data.frame(rbind(table(dt$ZONE, dt$VTYPE, dt$ECOT)[1,,],
#                           table(dt$ZONE, dt$VTYPE, dt$ECOT)[2,,],
#                           table(dt$ZONE, dt$VTYPE, dt$ECOT)[3,,]),
#                     VTYPE = rep(c('INDEL', 'SNP'), 3),
#                     ZONE = c('CZA', 'CZA', 'CZB', 'CZB', 'CZD', 'CZD'))

dt$class <- paste(dt$ZONE, dt$ECOT, dt$VTYPE, sep = ":")
(tot_v <- data.frame(table(dt$class)))

pr <- as.data.frame(rbindlist(lapply(tot_v$Var1, function(x) {
  dt1 <- dt[dt$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(Len = dt1$Length), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
head(pr)
# pr[pr$Len==329,]

vt <- split(pr, f = pr$VTYPE)
lapply(vt, head)
identical(vt$INDEL$Len, vt$SNP$Len)
setdiff(vt$INDEL$Len, vt$SNP$Len)
setdiff(vt$SNP$Len, vt$INDEL$Len)

mer <- merge(vt$INDEL, vt$SNP, by = c('ISL', 'ECOT', 'Len'), all = TRUE)
hist(mer$Len)
range(mer$Len)
mer$pal <- cut(mer$Len, breaks = 10, labels = FALSE)
table(mer$pal)
head(mer)
mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
table(mer$ECOT)
mer$ECOT <- factor(mer$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))

len_pal <- colorRampPalette(c("grey", "black"))
pprop <- ggplot(data = mer, aes(x = prop.x, y = prop.y, col = as.factor(pal))) +
  facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
  geom_abline(slope = 1, linetype = 'dashed') +
  geom_point() +
  scale_color_manual(values = len_pal(10)) +
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
pprop
# mer[which.max(mer$prop.x), ]
if (rm_inv) {
  ggsave(filename = "figures/len_snp_vs_indel_prop_noinv.pdf", plot = pprop, width = 10, height = 10, dpi = "print")
} else {
  ggsave(filename = "figures/len_snp_vs_indel_prop.pdf", plot = pprop, width = 10, height = 10, dpi = "print")
}

mer$class <- paste(mer$ISL, mer$ECOT, sep = ":")
head(mer)
invisible(lapply(unique(mer$class), function(x) {
  dt1 <- mer[mer$class==x, ]
  val <- cor(dt1$prop.x, dt1$prop.y)
  # val <- sum(dt1$prop.x, dt1$prop.y)
  cat(x, val, '\n')
}))


