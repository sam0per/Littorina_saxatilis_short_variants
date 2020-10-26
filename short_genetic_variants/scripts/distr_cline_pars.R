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
# fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
# fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# head(fai)
# Get inversion info
# invRui = read.table("/Users/samuelperini/Documents/research/projects/3.indels/data/20200123/Sweden_inversions_coordinates_2nd_august_2019.csv",
#                     sep=",", header=T, stringsAsFactors=F)
# invRui$LG = gsub("LG", "", invRui$LG)

n <- 1

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
  return(odt)
})))
head(cl_dt)

vtype_pal <- c("#1B9E77", "#666666")


cpars <- c("Centre", "Width", "slope", "p_diff", "Var.Ex")

library(dgof)
np <- 5
cpars[np]
vt <- split(x = cl_dt, f = cl_dt$VTYPE)
vt$INDEL <- vt$INDEL[!is.na(vt$INDEL[, cpars[np]]), ]
vt$SNP <- vt$SNP[!is.na(vt$SNP[, cpars[np]]), ]
ks.test(x = unique(vt$INDEL[, cpars[np]]), unique(vt$SNP[, cpars[np]]))

if (n==1) {
  lapply(X = cpars, function(x) {
    hm <- ggplot(data = cl_dt, aes_string(x = x, fill = "VTYPE")) +
      geom_histogram(bins = 40, position = "dodge") +
      scale_fill_manual(values = vtype_pal) +
      labs(fill = '') +
      theme(legend.position = 'none',
            # legend.text = element_text(size = 12),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.line = element_line(size = 0.2, linetype = "solid",
                                     colour = "black"),
            panel.grid = element_line(colour = "gray70", size = 0.2)) +
      guides(col = guide_legend(override.aes = list(size=3)))
    ggsave(filename = paste('figures/HB', x, 'indel_snp.pdf', sep = '_'), plot = hm, scale = 2/3, dpi = 'screen')
  })
}

(cpar <- cpars[n])
cl_dt_na <- cl_dt[!is.na(cl_dt[, cpar]), ]

range(cl_dt_na[, cpar])
nbins <- 31
cl_dt_na$cpar_bin <- cut(x = cl_dt_na[, cpar], breaks = seq(from = 0, to = round(max(cl_dt_na[, cpar])),
                                                            length.out = nbins),
                         include.lowest = TRUE)
table(cl_dt_na$cpar_bin)
clipar_bin <- merge(expand.grid(VTYPE = unique(cl_dt_na$VTYPE),
                                CLIPAR = levels(cl_dt_na$cpar_bin)),
                    aggregate(x = cl_dt_na[, cpar], list(VTYPE=cl_dt_na$VTYPE,
                                                         CLIPAR=cl_dt_na$cpar_bin),
                              function(x) c(Count = as.integer(length(x)), Mean_cpar = round(mean(x), 3))), all.x = TRUE)
clipar_bin[1:20,]

clipar_bin <- data.frame(clipar_bin[, 1:2], 
                         apply(X = clipar_bin$x, MARGIN = 2, FUN = function(x) ifelse(test = is.na(x), yes = 0, no = x)))
# contig_bin <- merge(x = contig_bin, y = fai, by = "Contig")
# contig_bin$Count_Len <- round(contig_bin$Count/contig_bin$Length, 5)

clipar_bin$Proportion <- NA
table(cl_dt_na$VTYPE)
clipar_bin[clipar_bin$VTYPE=="INDEL", "Proportion"] <- round(clipar_bin[clipar_bin$VTYPE=="INDEL", "Count"]/
                                                               table(cl_dt_na$VTYPE)["INDEL"], 9)
clipar_bin[clipar_bin$VTYPE=="SNP", "Proportion"] <- round(clipar_bin[clipar_bin$VTYPE=="SNP", "Count"]/
                                                             table(cl_dt_na$VTYPE)["SNP"], 9)
clipar_bin[1:30,]
clipar_bin <- clipar_bin[order(clipar_bin$CLIPAR), ]

cpar_split <- split(clipar_bin, f = clipar_bin$VTYPE)
identical(cpar_split$SNP$CLIPAR,cpar_split$INDEL$CLIPAR)
lapply(cpar_split, head)

freq_wide <- data.frame(CLIPAR=cpar_split$SNP[, "CLIPAR"],
                        INDEL=cpar_split$INDEL$Proportion,
                        INDEL_COUNT=cpar_split$INDEL$Count,
                        INDEL_MEAN_CPAR=cpar_split$INDEL$Mean_cpar,
                        SNP=cpar_split$SNP$Proportion,
                        SNP_COUNT=cpar_split$SNP$Count,
                        SNP_MEAN_CPAR=cpar_split$SNP$Mean_cpar)

freq_wide[1:20,]

# str(freq_wide)
# order(levels(freq_wide$FREQ))

cparlev <- c("[0,6.5]", "[0,6.23]", "[0,4.37]", "[0,0.0333]", "[0,2.8]")
freq_wide$CLIPAR <- relevel(freq_wide$CLIPAR, cparlev[n])

# cpar_prop <- ggplot(freq_wide, aes(x = INDEL, y = SNP)) +
#   facet_wrap(~CLIPAR) +
#   geom_point() +
#   geom_abline(slope = 1, linetype = "dashed") +
#   labs(x = paste0("Cline ", cpar, " INDEL proportion"), y = paste0("Cline ", cpar, " SNP proportion")) +
#   theme(axis.text.x = element_text(angle = 320, hjust = 0, size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 12),
#         legend.position = "top",
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "#91bfdb", color = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))

len_pal <- colorRampPalette(c("blueviolet", "seagreen1", "black", "red"))
# len_pal <- wes_palette(name = "Zissou1", n = nrow(freq_wide), type = "continuous")
# length(levels(freq_wide$CLIPAR))

cpar_prop <- ggplot(freq_wide, aes(x = INDEL, y = SNP, col = CLIPAR)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_color_manual(values = len_pal(nrow(freq_wide))) +
  # scale_color_manual(values = topo.colors(n = nrow(freq_wide))) +
  # scale_color_viridis_d() +
  # scale_color_manual(values = wes_palette(name = "Zissou1", n = nrow(freq_wide), type = "continuous")) +
  labs(x = paste0("Cline ", cpar, " INDEL proportion"), y = paste0("Cline ", cpar, " SNP proportion"), col = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  guides(colour = guide_legend(ncol = 2))

cpar_prop
ggsave(filename = paste0("figures/SNP_INDEL_prop_cline_", cpar, ".pdf"), plot = cpar_prop, width = 10, height = 7)

# length(unique(freq_wide$INDEL))
# length(unique(freq_wide$SNP))
# freq_wide$INDEL[duplicated(freq_wide$INDEL)]
# freq_wide$SNP[duplicated(freq_wide$SNP)]
ks.test(x = unique(freq_wide$INDEL), unique(freq_wide$SNP))
# ks.test(x = freq_wide$INDEL, freq_wide$SNP)
# ks.test(x = freq_wide$INDEL, freq_wide$SNP, exact = FALSE)