rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "reshape2", "tidyr", "tools", "data.table", "RColorBrewer", "dplyr", "textshape", "plotly",
              "devtools", "limma")
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

################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################
# Get cline fits
(cl_fl <- list.files(path = "CZCLI006_comp", full.names = TRUE))
(cl_fl <- cl_fl[grep(pattern = "ANG|NoInv", x = cl_fl, invert = TRUE)])
cl_ls <- lapply(cl_fl, read.table, header = TRUE)
# lapply(cl_ls, head)

cl_dt <- as.data.frame(rbindlist(lapply(seq_along(cl_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][2]
  side <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][3]
  zone <- paste(island, side, sep = "_")
  vtype <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][4]
  # paste(zone, vtype)
  
  odt <- mutate(cl_ls[[x]], ZONE = zone, VTYPE = vtype)
  sorted <- odt[order(odt$Var.Ex, decreasing = TRUE), ]
  # head(sorted)
  n_outl <- round((5 * nrow(sorted)) / 100)
  # cat('5% of', nrow(sorted), 'is', n_outl, '\n')
  sorted$sel[1:n_outl] <- TRUE
  
  # out_dt <- data.frame(table(sorted[sorted$sel==TRUE, 'invRui']),
  #                      ZONE = zone, VTYPE = vtype)
  out_dt <- data.frame(sorted[sorted$sel==TRUE, c('cp', 'LG', 'av', 'invRui')],
                       ZONE = zone, VTYPE = vtype)
  out_dt$av <- round(out_dt$av)
  return(out_dt)
})))
head(cl_dt)

table(cl_dt$ZONE, cl_dt$VTYPE)
# table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)[1,,,]
# table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)['CZA_left',,,]
# table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)[1:2,,,]
# table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)

# ggplot(data = cl_dt, aes(x = Var1, y = sqrt(Freq), fill = VTYPE)) +
#   facet_wrap(~ZONE) +
#   geom_col(position = position_dodge(preserve = "single")) +
#   labs(x = '', y = 'square root count', fill = '')

vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")
# lapply(X = levels(cl_dt$VTYPE), FUN = function(x) {
#   ggplot(data = cl_dt[cl_dt$VTYPE==x, ], aes(x = av, fill = VTYPE)) +
#     facet_grid(rows = vars(LG), cols = vars(ZONE)) +
#     geom_histogram(binwidth = 1) +
#     scale_fill_manual(values = as.character(vtype_pal[, x])) +
#     theme(axis.text = element_text(size = 9),
#           axis.title = element_text(size = 16),
#           strip.text = element_text(size = 12),
#           # legend.position = "top",
#           panel.background = element_blank(),
#           strip.background = element_rect(fill = "#91bfdb", color = "black"),
#           panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#           axis.line = element_line(size = 0.2, linetype = "solid",
#                                    colour = "black"),
#           panel.grid = element_line(colour = "gray70", size = 0.2))
# })

cl_agg <- aggregate(x = cl_dt$cp, by = list(LG = cl_dt$LG, av = cl_dt$av, ZONE = cl_dt$ZONE, VTYPE = cl_dt$VTYPE), length)
table(cl_dt$VTYPE)
cl_agg$prop <- round(ifelse(test = cl_agg$VTYPE=='INDEL', yes = cl_agg$x/528, no = cl_agg$x/3366) * 1000)
head(cl_agg)
tail(cl_agg)

vtype_pal <- c("#1B9E77", "#666666")
outl_scat <- ggplot(data = cl_agg, aes(y = ZONE, x = av, col = VTYPE)) +
  facet_wrap(~LG) +
  geom_point(aes(size = prop), position = position_jitter(height = 0.2), alpha = 0.7) +
  scale_color_manual(values = vtype_pal) +
  labs(x = 'map position', y = '', col = '', size = '') +
  theme(legend.text = element_text(size = 12),
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
  guides(col = guide_legend(override.aes = list(size=3)))
outl_scat
ggsave(filename = "figures/outlier_mappos.pdf", plot = outl_scat, width = 10, height = 8)


cl_agg <- aggregate(x = cl_dt$cp, by = list(LG = cl_dt$LG, av = cl_dt$av, VTYPE = cl_dt$VTYPE), length)
head(cl_agg)
sum(cl_agg[cl_agg$VTYPE=='INDEL', 'x'])
sum(cl_agg[cl_agg$VTYPE=='SNP', 'x'])
cl_agg$prop <- ifelse(test = cl_agg$VTYPE=='INDEL', yes = cl_agg$x/528, no = cl_agg$x/3366)

vt <- split(cl_agg, f = cl_agg$VTYPE)
# lapply(vt, head)

mer <- merge(vt$INDEL, vt$SNP, by = c('LG', 'av'), all = TRUE)
mer$prop.x <- ifelse(test = is.na(mer$prop.x), yes = 0, no = mer$prop.x)
mer$prop.y <- ifelse(test = is.na(mer$prop.y), yes = 0, no = mer$prop.y)
mer$x.x <- ifelse(test = is.na(mer$x.x), yes = 0, no = mer$x.x)
mer$x.y <- ifelse(test = is.na(mer$x.y), yes = 0, no = mer$x.y)
# unique(mer$VTYPE.x)
mer$VTYPE.x <- 'INDEL'
# unique(mer$VTYPE.y)
mer$VTYPE.y <- 'SNP'

outl_vt <- ggplot(data = mer, aes(x = prop.y, y = prop.x)) +
  facet_wrap(~LG) +
  geom_abline(slope = 1, linetype = 'dashed') +
  geom_point() +
  labs(x = 'Proportion of SNPs', y = 'Proportion of INDELs')
outl_vt
mer[mer$LG==6,]
cor(x = mer$x.x, y = mer$x.y)
cor(x = mer$prop.x, y = mer$prop.y)
