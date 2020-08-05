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
# BiocManager::install("limma")

franc_dt <- read.csv(file = "results/Lsax_short_var_czs_daf.csv")
head(franc_dt)

range(franc_dt$DAF)
franc_dt$DAF_bin <- cut(franc_dt$DAF, breaks = seq(from = 0, to = 1, length.out = 21), include.lowest = TRUE)
table(franc_dt$DAF_bin)
head(franc_dt)

dt <- split(x = franc_dt[, c('cp', 'ZONE', 'VTYPE', 'ECOT', 'DAF', 'DAF_bin')], f = franc_dt$ECOT)
lapply(dt, head)

popx = 'CRAB'
popy = 'WAVE_LEFT'
# head(dt[[wside]])
mer <- merge(x = dt[[popx]], dt[[popy]], by = c('cp', 'ZONE', 'VTYPE'))
head(mer)
table(mer$DAF_bin.x, mer$DAF_bin.y, mer$VTYPE, mer$ZONE)

mer$class <- paste(mer$ZONE, mer$VTYPE, sep = ":")
head(mer)
tot_v <- data.frame(table(mer$class))

pr <- as.data.frame(rbindlist(lapply(levels(tot_v$Var1), function(x) {
  dt1 <- mer[mer$class==x, ]
  
  dt2 <- data.frame(class = x,
                    aggregate(dt1$cp, by = list(DAF_popx = dt1$DAF_bin.x, DAF_popy = dt1$DAF_bin.y), length))
  dt2$prop <- dt2$x / tot_v[tot_v$Var1==x, 'Freq']
  
  dt2 <- separate(data = dt2, col = class, into = c('ISL', 'VTYPE'), sep = ":")
  
  # head(dt2)
  return(dt2)
})))
head(pr)

sqc_c <- ggplot(pr, aes(DAF_popx, DAF_popy)) +
  facet_grid(rows = vars(ISL), cols = vars(VTYPE)) +
  geom_tile(aes(fill = sqrt(prop))) +
  scale_fill_viridis_c() +
  labs(x = paste(popx, 'population uAFS'), y = paste(popy, 'population uAFS'), fill = 'relative sqrt proportion') +
  theme(legend.position = 'top',
        axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 320, hjust = 0),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        # legend.position = "top",
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
  # guides(fill = guide_legend(override.aes = list(size=3), nrow = 1))
sqc_c
ggsave(filename = paste("figures/jafs", popx, popy, "sqrt_prop.pdf", sep = "_"), plot = sqc_c, dpi = "print")

