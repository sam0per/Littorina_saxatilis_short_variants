rm(list = ls())
library(tidyr)
library(dplyr)
library(tools)
library(data.table)
library(ggplot2)
library(RColorBrewer)
# islands <- c("CZA", "CZB", "CZD")
# vtype <- c("INDEL", "SNP")

## AFTER FILTERING

# DERIVED ALLELE FREQUENCY

# gen_fl <- list.files(path = "genotypes", full.names = TRUE)
# genos <- lapply(gen_fl, read.table, header = TRUE)
# lapply(genos, head)

# MINOR ALLELE FREQUENCY

getwd()
(frq_fl <- list.files(path = "allele_freq", full.names = TRUE))
frqs <- lapply(frq_fl, read.table, header = TRUE)
lapply(frqs, head)

frqs_dt <- as.data.frame(rbindlist(lapply(seq_along(frq_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][2]
  vtype <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][3]
  
  one_frq <- separate(frqs[[x]], col = REF_FREQ, into = c("REF", "RFREQ"), sep = ":")
  one_frq <- separate(one_frq, col = ALT_FREQ, into = c("ALT", "AFREQ"), sep = ":")
  # head(one_frq)
  maf <- as.numeric(apply(one_frq[, c("RFREQ", "AFREQ")], MARGIN = 1, FUN = min))
  one_frq <- data.frame(ISL = island, VTYPE = vtype, MAF = maf,
                        cp = paste(one_frq[, "CHROM"], one_frq[, "POS"], sep = "_"))
  one_frq <- mutate(one_frq, ISL = as.character(ISL), VTYPE = as.character(VTYPE), cp = as.character(cp))
  return(one_frq)
  
  # hist(maf, breaks = 20, main = paste(vtype, "maf"))
  # abline(v = 0.1, col = "red")
  # hist(maf[maf>0.1], breaks = 20, main = paste(vtype, "maf"))
})))
head(frqs_dt)
str(frqs_dt)
# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(n = 8, name = 'Dark2')
# brewer.pal(n = 8, name = 'Dark2')
# display.brewer.pal(n = 9, name = 'Greys')
# brewer.pal(n = 9, name = 'Greys')
# display.brewer.pal(n = 9, name = 'Greens')
# brewer.pal(n = 9, name = 'Greens')

vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")

lapply(X = c("INDEL", "SNP"), FUN = function(y) {
  hist_p <- ggplot(data = frqs_dt[frqs_dt$VTYPE==y, ], aes(x = MAF)) +
    facet_grid(rows = vars(VTYPE)) +
    geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed", size = 1.5) +
    geom_histogram(aes(y = ..count../sum(..count..), fill = VTYPE), position = "dodge", binwidth = 0.01, col = "black") +
    scale_fill_manual(values = as.character(vtype_pal[, y])) +
    labs(fill = "") +
    scale_y_continuous(limits = c(0, 0.13)) +
    theme(legend.position = "top", legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("figures/maf_hist_prop_filt2_", y, ".pdf"), plot = hist_p, width = 10, height = 5)
})

# MARKER DENSITY

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

## AFTER FILTERING AND CLINE ANALYSIS

getwd()
# setwd("Anja/Anja_results/20200115")
# setwd("../../../")
(comp_fl <- list.files(path = "CZCLI006_comp", pattern = ".txt", full.names = TRUE))
(comp_fl <- comp_fl[grep(pattern = "NoInv", x = comp_fl, invert = TRUE)])
comps <- lapply(comp_fl, read.table, header = TRUE)
# lapply(comps, head)
lapply(comps, nrow)
# lapply(frqs, head)

comp_dt <- as.data.frame(rbindlist(lapply(seq_along(comp_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][2]
  side <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][3]
  vtype <- strsplit(file_path_sans_ext(basename(comp_fl[[x]])), split = "_")[[1]][4]
  # c(island, side, vtype)
  # colnames(comps[[x]])[2:3] <- c("CHROM", "POS")
  one_comp <- mutate(.data = comps[[x]], ISL = island, SIDE = side, VTYPE = vtype)
  one_comp$cp <- as.character(one_comp$cp)
  one_comp$Contig <- as.character(one_comp$Contig)
  return(one_comp)
})))
head(comp_dt)

cline_pal = list(INDEL=c("#00441B", "#74C476"),
                 SNP=c("#000000", "#969696"))
vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")

# MARKER DENSITY AGAINST CONTIG LENGTH

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# head(fai)
colnames(fai) <- c("Contig", "Length")
fai$Contig <- as.character(fai$Contig)

# comp_fai <- merge(comp_freq, fai, by = "Contig")
comp_fai <- as.data.frame(merge(comp_dt, fai, by = "Contig"))
head(comp_fai)
# str(comp_fai)

if (grepl(pattern = "indels", x = basename(getwd()))) {
  tool <- "gatk"
} else {
  tool <- "sam_gatk"
}

mean_len <- lapply(X = c("INDEL", "SNP"), FUN = function(y) {
  agg_one <- aggregate(x = comp_fai[comp_fai$VTYPE==y, "Length"], list(Type=comp_fai[comp_fai$VTYPE==y, "Type"]), mean)
  colnames(agg_one)[2] <- "Length"
  agg_one <- mutate(.data = agg_one, VTYPE = y)
  return(agg_one)
})
names(mean_len) <- c("INDEL", "SNP")

lapply(X = c("INDEL", "SNP"), FUN = function(y) {
  prop_hist <- ggplot(data = comp_fai[comp_fai$VTYPE==y, ], aes(x = Length)) +
    facet_grid(rows = vars(VTYPE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_histogram(aes(y = ..count../sum(..count..), fill = VTYPE, alpha = Type), binwidth = 15000, position = "dodge") +
    geom_vline(data = mean_len[[y]], aes(xintercept = Length, col = VTYPE, alpha = Type)) +
    scale_fill_manual(values = as.character(vtype_pal[, y])) +
    scale_color_manual(values = as.character(vtype_pal[, y])) +
    scale_alpha_discrete(range = c(1, 0.4)) +
    labs(fill = "", col = "", alpha = "", x = "contig length", y = paste(y, "proportion")) +
    theme(legend.position = "top", legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("/Users/samuelperini/Documents/research/projects/3.indels/figures/prop_cline_", y ,"_", tool, "_contigs.pdf"),
         plot = prop_hist, width = 10, height = 5)
})


# MINOR ALLELE FREQUENCY

head(frqs_dt)
str(comp_dt)
str(frqs_dt)
comp_freq <- merge(comp_dt, frqs_dt, by = c("cp", "ISL", "VTYPE"))
head(comp_freq)
# round(min(comp_freq$MAF), 2)

lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  # one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$MAF>0.15, ]
  one_vtype <- comp_freq[comp_freq$VTYPE==x, ]
  # min(one_vtype[, "MAF"])
  hist_one <- ggplot(data = one_vtype, aes(x = MAF)) +
    facet_grid(rows = vars(VTYPE)) +
    # facet_grid(rows = vars(ISL), cols = vars(SIDE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed") +
    geom_histogram(aes(y = ..count../sum(..count..),fill = VTYPE, alpha = Type),
                   position = "dodge", binwidth = 0.05, col = "black") +
    scale_fill_manual(values = as.character(vtype_pal[, x])) +
    scale_alpha_discrete(range = c(1, 0.4)) +
    scale_y_continuous(limits = c(0, 0.1)) +
    labs(fill = "", alpha = "", x = paste(x, "MAF")) +
    theme(legend.position = "top", legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 11),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 14),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("figures/maf_hist_prop_cline_", x, ".pdf"), plot = hist_one, width = 10, height = 5)
})

cline_type <- "Cline"
isla <- "CZB"
side <- "right"
split_vtype <- split(comp_freq[comp_freq$Type==cline_type, ], f = as.factor(comp_freq[comp_freq$Type==cline_type, "VTYPE"]))
lapply(split_vtype, head)
split_vtype <- lapply(split_vtype, function(x) {
  x[, "MAF"] <- round(x[, "MAF"], 2)
  # x <- x[x[, "ISL"]==isla & x[, "SIDE"]==side, ]
  return(x)
})
range(split_vtype$INDEL$MAF)
range(split_vtype$SNP$MAF)
grid_dt <- expand.grid(INDEL = seq(from = min(split_vtype$INDEL$MAF), to = max(split_vtype$INDEL$MAF), by = 0.05),
                       SNP = seq(from = min(split_vtype$SNP$MAF), to = max(split_vtype$SNP$MAF), by = 0.05))
head(grid_dt)
sum(split_vtype$INDEL$MAF==grid_dt[1,1])/nrow(split_vtype$INDEL) +
  sum(split_vtype$SNP$MAF==grid_dt[1,2])/nrow(split_vtype$SNP)
sum(split_vtype$INDEL$MAF==grid_dt[4,1])/nrow(split_vtype$INDEL) +
  sum(split_vtype$SNP$MAF==grid_dt[4,2])/nrow(split_vtype$SNP)
sum(split_vtype$INDEL$MAF==grid_dt[2,1], split_vtype$SNP$MAF==grid_dt[2,2])

grid_dt$count <- apply(X = grid_dt, MARGIN = 1, FUN = function(x) {
  sum(split_vtype$INDEL$MAF==x[1])/nrow(split_vtype$INDEL) + sum(split_vtype$SNP$MAF==x[2])/nrow(split_vtype$SNP)
})
head(grid_dt)
grid_dt[1:15,]
grid_dt$INDEL <- as.factor(grid_dt$INDEL)
grid_dt$SNP <- as.factor(grid_dt$SNP)

ggp <- ggplot(grid_dt, aes(INDEL, SNP)) +
  geom_tile(aes(fill = count))
ggp

# wide_dt <- reshape(grid_dt, idvar = "INDEL", timevar = "SNP", direction = "wide")
wide_dt <- reshape2::dcast(grid_dt, SNP ~ INDEL)
head(wide_dt)
# install.packages('textshape')
library(textshape)
wide_dt <- as.matrix(column_to_rownames(wide_dt, 'SNP'))
# rownames(wide_dt) <- paste0('row', 1:nrow(wide_dt))
# colnames(wide_dt) <- paste0('col', 1:ncol(wide_dt))
plot_ly(x = rownames(wide_dt), y = colnames(wide_dt), z = wide_dt, colors = "Greys", type = "heatmap") %>%
  add_segments(x = 0.1, y = 0.1, xend = 0.5, yend = 0.5) %>%
  layout(xaxis = list(title = "INDEL MAF"),
         yaxis = list(title = "SNP MAF"))

write.table(x = wide_dt, file = 'results/joint_mafs_CZs.txt', append = FALSE, quote = FALSE, sep = '\t', row.names = TRUE,
            col.names = TRUE)
# class(wide_dt)
# wide_dt <- as.matrix(wide_dt)

# split_vtype$INDEL$MAF[split_vtype$INDEL$MAF==0.15]
# split_vtype$SNP$MAF[split_vtype$SNP$MAF==0.15]
# table(split_vtype$INDEL$MAF)
# table(split_vtype$SNP$MAF)
# re_freq <- data.frame()

head(comp_freq)

prop_vtype <- lapply(X = c("INDEL", "SNP"), FUN = function(x) {
  one_vtype <- comp_freq[comp_freq$VTYPE==x & comp_freq$Type==cline_type, ]
  bin_seq <- seq(from = round(min(one_vtype$MAF), 1), to = round(max(one_vtype$MAF), 1), by = 0.05)
  one_vtype$bins <- cut(one_vtype$MAF, breaks = bin_seq, labels = FALSE, include.lowest = TRUE)
  # head(one_vtype)
  maf_prop <- aggregate(x = one_vtype$MAF, list(MAF_bin = one_vtype$bins), function(y) length(y)/nrow(one_vtype))
  maf_prop$bin_seq <- bin_seq[-length(bin_seq)]
  maf_prop$VTYPE <- x
  
  return(maf_prop)
  # table(one_vtype$bins)
  # nrow(one_vtype)
})
sum(prop_vtype[[1]]$x)
sum(prop_vtype[[2]]$x)

prop_maf <- data.frame(INDEL_PROP = prop_vtype[[1]]$x, SNP_PROP = prop_vtype[[2]]$x)

ggplot(data = prop_maf, aes(x = INDEL_PROP, y = SNP_PROP)) +
  geom_abline(slope = 1) +
  geom_point()

# hist_one <- ggplot(data = one_vtype, aes(x = MAF)) +
#   facet_grid(rows = vars(VTYPE)) +
#   geom_histogram(aes(y = ..count../sum(..count..),fill = VTYPE, alpha = Type),
#                  position = "dodge", binwidth = 0.05, col = "black") +
#   scale_fill_manual(values = as.character(vtype_pal[, x])) +
#   scale_alpha_discrete(range = c(1, 0.4)) +
#   scale_y_continuous(limits = c(0, 0.1)) +
#   labs(fill = "", alpha = "", x = paste(x, "MAF")) +
#   theme(legend.position = "top", legend.text = element_text(size = 14),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 11),
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "#91bfdb", color = "black"),
#         strip.text = element_text(size = 14),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))
# ggsave(filename = paste0("figures/maf_hist_prop_cline_", x, ".pdf"), plot = hist_one, width = 10, height = 5)



