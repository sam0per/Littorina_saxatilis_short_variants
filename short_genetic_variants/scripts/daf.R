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

rm_inv = TRUE
franc_dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv.csv")
head(franc_dt)
if (rm_inv) {
  franc_dt <- franc_dt[franc_dt$invRui==FALSE, ]
}
table(franc_dt$invRui)

# franc_dt <- read.csv(file = "results/Lsax_short_var_czs_daf.csv")
# head(franc_dt)

range(franc_dt$DAF)
franc_dt$DAF_bin <- cut(franc_dt$DAF, breaks = seq(from = 0, to = 1, length.out = 21), include.lowest = TRUE)
table(franc_dt$DAF_bin)
head(franc_dt)
franc_dt$class <- paste(franc_dt$ZONE, franc_dt$ECOT, franc_dt$VTYPE, sep = ":")
tot_v <- data.frame(table(franc_dt$class))

aggr_dt <- aggregate(x = franc_dt$cp, by=list(DAF_bin=franc_dt$DAF_bin, CLASS=franc_dt$class), length)

pr <- as.data.frame(rbindlist(lapply(levels(tot_v$Var1), function(x) {
  dt1 <- aggr_dt[aggr_dt$CLASS==x, ]
  
  dt1$prop <- dt1$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt1 <- separate(data = dt1, col = CLASS, into = c('ISL', 'ECOT', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt1)
})))
pr
str(pr)
pr$ECOT <- factor(pr$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))

vtype_pal <- c("#1B9E77", "#666666")
pdaf <- ggplot(data = pr, aes(x = DAF_bin, y = prop, fill = VTYPE)) +
  facet_grid(rows = vars(ISL), cols = vars(ECOT)) +
  geom_col(col = 'black', position = 'dodge') +
  scale_fill_manual(values = vtype_pal) +
  labs(x = 'DAF', y = 'relative proportion', fill = '') +
  theme(legend.text = element_text(size = 12), legend.position = 'top',
        axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 7, angle = 320, hjust = 0),
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
pdaf
if (rm_inv) {
  ggsave(filename = "figures/daf_hist_prop_noinv.pdf", plot = pdaf, width = 10, height = 9, dpi = "print")
} else {
  ggsave(filename = "figures/daf_hist_prop.pdf", plot = pdaf, width = 10, height = 9, dpi = "print")
}


head(pr)
pr$class <- paste(pr$ISL, pr$ECOT, sep = ":")
lapply(unique(pr$class), function(x) {
  dt1 <- pr[pr$class==x, ]
  dt1 <- split(x = dt1, f = dt1$VTYPE)
  cat(x, '\n')
  # head(dt1)
  # c(nrow(dt1), length(unique(dt1$prop)))
  
  ks.test(x = unique(dt1$INDEL$prop), y = unique(dt1$SNP$prop))
  # ks.test(x = sqrt(unique(dt1$INDEL$x)), y = sqrt(unique(dt1$SNP$x)))
})
