rm(list = ls())

# islands <- c("CZA", "CZB", "CZD")
# vtype <- c("INDEL", "SNP")

## AFTER FILTERING

(frq_fl <- list.files(path = "allele_freq", full.names = TRUE))
frqs <- lapply(frq_fl, read.table, header = TRUE)
lapply(frqs, head)

# gen_fl <- list.files(path = "genotypes", full.names = TRUE)
# genos <- lapply(gen_fl, read.table, header = TRUE)
# lapply(genos, head)

library(tidyr)
library(dplyr)
library(tools)
library(data.table)
library(ggplot2)
library(RColorBrewer)


frqs_dt <- rbindlist(lapply(seq_along(frq_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][2]
  vtype <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][3]
  
  one_frq <- separate(frqs[[x]], col = REF_FREQ, into = c("REF", "RFREQ"), sep = ":")
  one_frq <- separate(one_frq, col = ALT_FREQ, into = c("ALT", "AFREQ"), sep = ":")
  # head(one_frq)
  maf <- as.numeric(apply(one_frq[, c("RFREQ", "AFREQ")], MARGIN = 1, FUN = min))
  one_frq <- data.frame(ISL = island, VTYPE = vtype, MAF = maf, cp = paste(one_frq[, "CHROM"], one_frq[, "POS"], sep = "_"))
  
  return(one_frq)
  
  # hist(maf, breaks = 20, main = paste(vtype, "maf"))
  # abline(v = 0.1, col = "red")
  # hist(maf[maf>0.1], breaks = 20, main = paste(vtype, "maf"))
}))
head(frqs_dt)
# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(n = 8, name = 'Dark2')
# brewer.pal(n = 8, name = 'Dark2')
display.brewer.pal(n = 9, name = 'Greys')
brewer.pal(n = 9, name = 'Greys')
display.brewer.pal(n = 9, name = 'Greens')
brewer.pal(n = 9, name = 'Greens')

(hist_p <- ggplot(data = frqs_dt, aes(x = MAF)) +
  facet_grid(ISL ~ .) +
  geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
  geom_histogram(aes(fill = VTYPE), position = "dodge", binwidth = 0.05) +
  scale_fill_manual(values = c("#1B9E77", "#666666")) +
  labs(fill = "") +
  theme(legend.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        strip.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(filename = "figures/maf_hist_filt2.pdf", plot = hist_p)

(den_p <- ggplot(data = frqs_dt, aes(x = MAF)) +
    facet_grid(ISL ~ .) +
    geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_density(aes(col = VTYPE)) +
    scale_color_manual(values = c("#1B9E77", "#666666")) +
    labs(col = "") +
    theme(legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(filename = "figures/maf_den_filt2.pdf", plot = den_p)

(sca_p <- ggplot(data = frqs_dt, aes(x = MAF)) +
    facet_grid(ISL ~ .) +
    geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_density(aes(y = ..scaled.., col = VTYPE)) +
    scale_color_manual(values = c("#1B9E77", "#666666")) +
    labs(col = "") +
    theme(legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2)))
ggsave(filename = "figures/maf_scaled_den_filt2.pdf", plot = sca_p)

## AFTER FILTERING AND CLINE ANALYSIS

(comp_fl <- list.files(path = "CZCLI006_comp", full.names = TRUE))
comp_fl <- comp_fl[grep(pattern = "NoInv", x = comp_fl, invert = TRUE)][-1:-2]
comps <- lapply(comp_fl, read.table, header = TRUE)
lapply(comps, head)
lapply(comps, nrow)
lapply(frqs, head)

comp_dt <- rbindlist(lapply(seq_along(comp_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][2]
  side <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][3]
  vtype <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][4]
  # c(island, side, vtype)
  # colnames(comps[[x]])[2:3] <- c("CHROM", "POS")
  one_comp <- mutate(.data = comps[[x]], ISL = island, SIDE = side, VTYPE = vtype)
  return(one_comp)
}))
head(comp_dt)
head(frqs_dt)

comp_freq <- merge(comp_dt, frqs_dt, by = c("cp", "ISL", "VTYPE"))
head(comp_freq)
round(min(comp_freq$MAF), 2)

cline_pal = list(INDEL=c("#00441B", "#74C476"),
                 SNP=c("#000000", "#969696"))

lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  # one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$MAF>0.15, ]
  one_vtype <- comp_freq[comp_freq$VTYPE==x, ]
  # min(one_vtype[, "MAF"])
  hist_one <- ggplot(data = one_vtype, aes(x = MAF)) +
    facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_histogram(aes(fill = Type), position = "dodge", binwidth = 0.05) +
    scale_fill_manual(values = cline_pal[[x]]) +
    labs(fill = "", x = paste(x, "MAF")) +
    theme(legend.position = "top", legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("figures/maf_hist_cline_", x, ".pdf"), plot = hist_one)
})

lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  # one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$MAF>0.15, ]
  one_vtype <- comp_freq[comp_freq$VTYPE==x, ]
  # min(one_vtype[, "MAF"])
  dens_one <- ggplot(data = one_vtype, aes(x = MAF)) +
    facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_density(aes(col = Type)) +
    scale_color_manual(values = cline_pal[[x]]) +
    labs(col = "", x = paste(x, "MAF")) +
    theme(legend.position = "top", legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("figures/maf_den_cline_", x, ".pdf"), plot = dens_one)
})

head(comp_freq)
table(comp_freq$Type)

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
head(fai)
colnames(fai) <- c("Contig", "Length")

comp_fai <- merge(comp_freq, fai, by = "Contig")
lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  # one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$MAF>0.15, ]
  one_vtype <- comp_fai[comp_fai$VTYPE==x, ]
  # min(one_vtype[, "MAF"])
  dens_one <- ggplot(data = one_vtype, aes(x = Length)) +
    facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_density(aes(col = Type)) +
    scale_color_manual(values = cline_pal[[x]]) +
    labs(col = "", x = "Contig length", y = paste(x, "density")) +
    theme(legend.position = "top", legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  # ggsave(filename = paste0("figures/maf_den_cline_", x, ".pdf"), plot = dens_one)
})

lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  # one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$MAF>0.15, ]
  one_vtype <- comp_fai[comp_fai$VTYPE==x, ]
  # min(one_vtype[, "MAF"])
  dens_one <- ggplot(data = one_vtype, aes(x = Length)) +
    facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_histogram(aes(y = ..density.., fill = Type), position = "dodge", binwidth = 15000) +
    scale_fill_manual(values = cline_pal[[x]]) +
    labs(fill = "", x = "Contig length", y = paste(x, "density")) +
    theme(legend.position = "top", legend.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  # ggsave(filename = paste0("figures/maf_den_cline_", x, ".pdf"), plot = dens_one)
})

(dens_hist <- ggplot(data = comp_fai, aes(x = Length)) +
  facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
  # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
  geom_histogram(aes(y = ..density.., fill = VTYPE), position = "dodge", binwidth = 15000) +
  # scale_fill_manual(values = cline_pal[[x]]) +
  labs(fill = "", x = "Contig length") +
  theme(legend.position = "top", legend.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        strip.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)))
# ggsave(filename = paste0("figures/maf_den_cline_", x, ".pdf"), plot = dens_one)

