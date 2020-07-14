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

## AFTER FILTERING - BEFORE CLINE ANALYSIS

# DERIVED ALLELE FREQUENCY

(anc_fl <- list.files(path = "allele_freq", pattern = "ANC", full.names = TRUE))
ancs <- as.data.frame(rbindlist(lapply(anc_fl, read.csv)))
head(ancs)
ancs$NE_F1_141_Lc <- ifelse(test = ancs$NE_F1_141_Lc=="-het", yes = NA, no = as.character(ancs$NE_F1_141_Lc))
ancs$W_com_01_Lc <- ifelse(test = ancs$W_com_01_Lc=="-het", yes = NA, no = as.character(ancs$W_com_01_Lc))

(frq_fl <- list.files(path = "summary/allele_freq", pattern = ".txt", full.names = TRUE))
frqs <- lapply(frq_fl, read.table, header = TRUE)
lapply(frqs, head)

frqs_dt <- as.data.frame(rbindlist(lapply(seq_along(frq_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][2]
  ecot <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][3]
  vtype <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][4]
  
  one_frq <- separate(frqs[[x]], col = REF_FREQ, into = c("REF", "RFREQ"), sep = ":")
  one_frq <- separate(one_frq, col = ALT_FREQ, into = c("ALT", "AFREQ"), sep = ":")
  
  # one_frq <- data.frame(ISL = island, VTYPE = vtype, MAF = maf,
  #                       cp = paste(one_frq[, "CHROM"], one_frq[, "POS"], sep = "_"))
  one_frq <- mutate(one_frq, ZONE = as.character(island), ECOT = as.character(ecot),
                    VTYPE = as.character(strsplit(x = vtype, split = ".", fixed = TRUE)[[1]][1]),
                    cp = paste(one_frq[, "CHROM"], one_frq[, "POS"], sep = "_"))
  # one_frq <- separate(data = one_frq, col = cp, into = c("Contig", "Position"), sep = "_", remove = FALSE)
  return(one_frq)
})))
head(frqs_dt)

franc_dt <- merge(x = ancs, y = frqs_dt, by = c("cp", "ZONE", "VTYPE"))
franc_dt <- franc_dt[rowSums(apply(X = franc_dt[, c("NE_F1_141_Lc", "W_com_01_Lc")], MARGIN = 2,
                                   FUN = function(x) is.na(x))) != 2, ]

head(franc_dt)

rm(list = setdiff(x = ls(), y = "franc_dt"))

table(franc_dt$NE_F1_141_Lc)

outg <- which(colnames(franc_dt) %in% c("NE_F1_141_Lc", "W_com_01_Lc"))
# combn(x = unique(franc_dt[, outg[1]]), m = 2)[,4]
franc_dt$NE_W_Lcomp <- paste(franc_dt$NE_F1_141_Lc, franc_dt$W_com_01_Lc, sep = ":")
anc_tb <- data.frame(table(franc_dt$VTYPE, franc_dt$NE_W_Lcomp))
anc_tb <- separate(data = anc_tb, col = Var2, into = c("NE", "W"), sep = ":")
anc_tb <- split(x = anc_tb, f = anc_tb$Var1)
anc_tb
# View(anc_tb$INDEL)
# View(anc_tb$SNP)

# franc_dt <- franc_dt[franc_dt$NE_W_Lcomp!="alt_anc:ref_anc" & franc_dt$NE_W_Lcomp!="ref_anc:alt_anc" &
#                        franc_dt$NE_W_Lcomp!="het:het" & franc_dt$NE_W_Lcomp!="het:NA" &
#                        franc_dt$NE_W_Lcomp!="NA:het", ]
franc_dt <- franc_dt[franc_dt$NE_W_Lcomp=="alt_anc:alt_anc" | franc_dt$NE_W_Lcomp=="ref_anc:ref_anc", ]
data.frame(table(franc_dt$NE_W_Lcomp))
franc_dt$DAF <- NA
franc_dt$DAF <- ifelse(test = grepl(pattern = "alt", x = franc_dt$NE_W_Lcomp), yes = franc_dt$RFREQ, no = franc_dt$AFREQ)
sample_n(tbl = franc_dt, size = 10)
franc_dt$DAF <- as.numeric(franc_dt$DAF)
str(franc_dt)

range(franc_dt$DAF)
franc_dt$DAF_bin <- cut(franc_dt$DAF, breaks = seq(from = 0, to = 1, length.out = 21), include.lowest = TRUE)
table(franc_dt$DAF_bin)
grid_dt <- expand.grid(CHROM=levels(franc_dt$CHROM), DAF_bin=levels(franc_dt$DAF_bin),
                       ZONE=levels(franc_dt$ZONE), ECOT=unique(franc_dt$ECOT), VTYPE=levels(franc_dt$VTYPE))
sample_n(tbl = grid_dt, size = 10)

aggr_dt <- aggregate(x = franc_dt$DAF, by=list(CHROM=franc_dt$CHROM, DAF_bin=franc_dt$DAF_bin, ZONE=franc_dt$ZONE,
                                               ECOT=franc_dt$ECOT, VTYPE=franc_dt$VTYPE), function(x) {
                                                 c(count=length(x), mean_daf=mean(x))
                                               })
aggr_dt[aggr_dt$CHROM=="Contig0", ]
table(franc_dt[franc_dt$CHROM=="Contig0", "VTYPE"])
# franc_dt$DAF <- round(franc_dt$DAF, 1)
sample_n(tbl = franc_dt, size = 10)
sample_n(tbl = aggr_dt, size = 10)

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# head(fai)
colnames(fai) <- c("CHROM", "Length")

contig_bin <- merge(grid_dt, aggr_dt, all.x = TRUE)
contig_bin[1:20,]
# contig_bin[contig_bin$Contig=="Contig0", ]
# frqs_fai[frqs_fai$Contig=="Contig0", ]
contig_bin <- data.frame(contig_bin[, 1:5], 
                         apply(X = contig_bin$x, MARGIN = 2, FUN = function(x) ifelse(test = is.na(x), yes = 0, no = x)))
contig_bin <- merge(x = contig_bin, y = fai, by = "CHROM")

str(contig_bin)
table(contig_bin$ZONE, contig_bin$ECOT, contig_bin$VTYPE)

corfun <- function(dd, pop, eco, type) {
  idx <- rowSums(data.frame(pop=dd[,1]==pop, eco=dd[,2]==eco, type=dd[,3]==type))==3
  dt <- aggregate(x = dd[idx, "count"], by = list(dd[idx, "Length"]), sum)
  # return(dt)
  corrv <- cor(x = dt[,1], y = dt[,2])
  # cat("For", pop, eco, type, ", correlation between count and contig length is\n", corrv, "\n")
  corplot <- ggplot(data = dt, aes(x = Group.1, y = x)) +
    geom_point() +
    labs(x = "contig length", y = paste(type, "count"),
         title = paste(pop, eco, type, "correlation =", round(corrv, 2))) +
    theme(title = element_text(size = 16), axis.text = element_text(size = 12))
  # geom_text(x=200, y=15, label=corrv)
  ggsave(filename = paste("figures/corr", pop, eco, type, "count_length.pdf", sep = "_"), plot = corplot)
}
# corfun(dd = contig_bin[, -1:-2], pop = "CZA", eco = "CRAB", type = "INDEL")
# corfun(dd = contig_bin[, -1:-2], pop = "CZA", eco = "CRAB", type = "SNP")
# 
# contig_bin[contig_bin$ZONE=="CZA" & contig_bin$ECOT=="CRAB" & contig_bin$VTYPE=="INDEL", ]
# contig_bin[contig_bin$ZONE=="CZA" & contig_bin$ECOT=="CRAB" & contig_bin$VTYPE=="INDEL" & contig_bin$Length==12830, ]
# cor(x = contig_bin[contig_bin$ZONE=="CZA" & contig_bin$ECOT=="CRAB" & contig_bin$VTYPE=="INDEL", "count"],
#     y = contig_bin[contig_bin$ZONE=="CZA" & contig_bin$ECOT=="CRAB" & contig_bin$VTYPE=="INDEL", "Length"])
# 
# contig_bin[579:581,]
comb_dt <- expand.grid(levels(contig_bin$ZONE), levels(contig_bin$ECOT), levels(contig_bin$VTYPE))
apply(X = comb_dt, MARGIN = 1, FUN = function(x) {
  corfun(dd = contig_bin[, -1:-2], pop = x[1], eco = x[2], type = x[3])
})

contig_bin[1:20,]
contig_bin <- contig_bin[order(contig_bin$DAF_bin), ]
daf_sum <- aggregate(contig_bin$count, by = list(ZONE=contig_bin$ZONE, ECOT=contig_bin$ECOT,
                                                 VTYPE=contig_bin$VTYPE, DAF_bin=contig_bin$DAF_bin), sum)
head(daf_sum)

vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")
lapply(X = c("INDEL", "SNP"), FUN = function(y) {
  hist_p <- ggplot(data = daf_sum[daf_sum$VTYPE==y, ], aes(x = as.integer(DAF_bin)-1)) +
    facet_grid(rows = vars(ECOT), cols = vars(ZONE)) +
    # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed", size = 1.5) +
    geom_histogram(aes(y = sqrt(x), fill = VTYPE), stat = "identity", col = "black") +
    scale_fill_manual(values = as.character(vtype_pal[, y])) +
    labs(fill = "", x = "DAF", y = "squared count") +
    scale_y_continuous(limits = c(0, 200)) +
    # scale_x_discrete(labels = paste(data.frame(levels(daf_sum$DAF_bin), seq(0.05,1,0.05))[,1],
    #                                 data.frame(levels(daf_sum$DAF_bin), seq(0.05,1,0.05))[,2], sep = '"="')) +
    # scale_x_discrete(breaks = seq(0, to = 1, length.out = 11), labels = seq(0, to = 1, length.out = 11)) +
    theme(legend.position = "top", legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 11),
          # axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
          panel.background = element_blank(),
          strip.background = element_rect(fill = "#91bfdb", color = "black"),
          strip.text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  ggsave(filename = paste0("figures/daf_hist_sqrt_count_filt2_", y, ".pdf"), plot = hist_p)
})
sample_n(tbl = daf_sum, size = 10)
sum(daf_sum$x)

zev_sum <- aggregate(daf_sum$x, by = list(ZONE=daf_sum$ZONE, ECOT=daf_sum$ECOT, VTYPE=daf_sum$VTYPE), sum)
str(zev_sum)
str(daf_sum)
library(prodlim)
row.match(c("CZB", "WAVE", "SNP"), zev_sum[,1:3])
daf_sum$count_prop <- apply(X = daf_sum, MARGIN = 1, FUN = function(x) {
  idx_tot <- row.match(x[1:3], zev_sum[,1:3])
  # idx_tot
  x_prop <- as.numeric(x[5]) / zev_sum[idx_tot, 4]
  return(x_prop)
})
head(daf_sum)
sum(daf_sum[daf_sum$ZONE=='CZA' & daf_sum$ECOT=='WAVE' & daf_sum$VTYPE=='INDEL', 'count_prop'])
sum(daf_sum[daf_sum$ZONE=='CZA' & daf_sum$ECOT=='WAVE', 'count_prop'])

hist_p <- ggplot(data = daf_sum, aes(x = as.integer(DAF_bin)-1)) +
  facet_grid(rows = vars(ECOT), cols = vars(ZONE)) +
  # geom_vline(xintercept = 0.1, col = "#D95F02", linetype = "dashed", size = 1.5) +
  geom_histogram(aes(y = count_prop, fill = VTYPE), stat = "identity", position = "dodge", col = "black") +
  scale_fill_manual(values = c("#1B9E77", "#666666")) +
  labs(fill = "", x = "DAF", y = "relative proportion") +
  # scale_y_continuous(limits = c(0, 200)) +
  # scale_x_discrete(labels = paste(data.frame(levels(daf_sum$DAF_bin), seq(0.05,1,0.05))[,1],
  #                                 data.frame(levels(daf_sum$DAF_bin), seq(0.05,1,0.05))[,2], sep = '"="')) +
  # scale_x_discrete(breaks = seq(0, to = 1, length.out = 11), labels = seq(0, to = 1, length.out = 11)) +
  theme(legend.position = "top", legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 11),
        # axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        strip.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
ggsave(filename = "figures/daf_hist_rel_prop_filt2.pdf", plot = hist_p)

# table(daf_sum$VTYPE, round(sqrt(daf_sum$x)))
# as.integer(daf_sum$DAF_bin)

colnames(daf_sum) <- c("grp1", "ECOT", "grp2", "DAF_bin", "count", "count_prop")
daf_ecot <- split(daf_sum, f = daf_sum$ECOT)
lapply(daf_ecot, head)
colnames(daf_ecot$WAVE)[5:6] <- c("count_pop1", "count_prop_pop1")
colnames(daf_ecot$CRAB)[5:6] <- c("count_pop2", "count_prop_pop2")
daf_ecot$WAVE$ECOT <- NULL
daf_ecot$CRAB$ECOT <- NULL

join2_daf <- function(dd1, dd2, colnm, grps) {
  cn1 <- paste0(colnames(dd1)[names(dd1) == colnm], "_pop1")
  cn2 <- paste0(colnames(dd2)[names(dd2) == colnm], "_pop2")
  # c(cn1, cn2)
  dt <- expand.grid(levels(dd1[, names(dd1) == colnm]),
                    levels(dd2[, names(dd2) == colnm]))
  colnames(dt) <- c(cn1, cn2)
  # head(dt)
  for (i in seq_along(grps)) {
    dt[, paste0("grp", i)] <- grps[i]
  }
  colnames(dd1)[names(dd1) == colnm] <- cn1
  colnames(dd2)[names(dd2) == colnm] <- cn2
  dtm <- merge(dt, dd1)
  dtm <- merge(x = dtm, y = dd2)
  # head(dd2)
  # dtm2 <- merge(dt, dd2)
  return(dtm)
}
# zone <- "CZA"
# vtype <- "SNP"
zv_grid <- expand.grid(zone=c("CZA", "CZB", "CZD"), vtype=c("INDEL", "SNP"))
apply(X = zv_grid, MARGIN = 1, FUN = function(x) {
  zv_dt <- join2_daf(dd1 = daf_ecot$WAVE, dd2 = daf_ecot$CRAB, colnm = "DAF_bin", grps = c(x[1], x[2]))
  
  zv_dt$DAF_bin_pop1 <- relevel(zv_dt$DAF_bin_pop1, "[0,0.05]")
  zv_dt$DAF_bin_pop2 <- relevel(zv_dt$DAF_bin_pop2, "[0,0.05]")
  
  zv_dt$sqrt_pop1 <- sqrt(zv_dt$count_pop1)
  zv_dt$sqrt_pop2 <- sqrt(zv_dt$count_pop2)
  zv_dt$sqrt_sum <- rowSums(x = zv_dt[, 9:10])
  
  sqc_b <- ggplot(zv_dt, aes(DAF_bin_pop1, DAF_bin_pop2)) +
    geom_tile(aes(fill = round(sqrt_sum))) +
    scale_fill_viridis_b() +
    labs(x = 'WAVE population uAFS', y = 'CRAB population uAFS', fill = paste(x[1], x[2], 'square root count')) +
    theme(axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
          legend.position = 'top')
  ggsave(filename = paste("figures/jafs", x[1], x[2], "sqrt_count_virb.pdf", sep = "_"), plot = sqc_b, width = 8, height = 8)
  
  sqc_c <- ggplot(zv_dt, aes(DAF_bin_pop1, DAF_bin_pop2)) +
    geom_tile(aes(fill = round(sqrt_sum))) +
    scale_fill_viridis_c() +
    labs(x = 'WAVE population uAFS', y = 'CRAB population uAFS', fill = paste(x[1], x[2], 'square root count')) +
    theme(axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
          legend.position = 'top')
  ggsave(filename = paste("figures/jafs", x[1], x[2], "sqrt_count_virc.pdf", sep = "_"), plot = sqc_c, width = 8, height = 8)
  
  p_colu <- grepl(pattern = 'prop', colnames(zv_dt))
  zv_dt$diff_count_prop <- apply(X = zv_dt[, p_colu], MARGIN = 1, FUN = function(y) y[1] - y[2])
  
  pdiff_b <- ggplot(zv_dt, aes(DAF_bin_pop1, DAF_bin_pop2)) +
    geom_tile(aes(fill = round(diff_count_prop, 2))) +
    scale_fill_viridis_b() +
    labs(x = 'WAVE population uAFS', y = 'CRAB population uAFS', fill = paste(x[1], x[2], 'prop. difference')) +
    theme(axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
          legend.position = 'top')
  ggsave(filename = paste("figures/jafs", x[1], x[2], "prop_diff_virb.pdf", sep = "_"), plot = pdiff_b, width = 8, height = 8)
  
  pdiff_c <- ggplot(zv_dt, aes(DAF_bin_pop1, DAF_bin_pop2)) +
    geom_tile(aes(fill = round(diff_count_prop, 2))) +
    scale_fill_viridis_c() +
    labs(x = 'WAVE population uAFS', y = 'CRAB population uAFS', fill = paste(x[1], x[2], 'prop. difference')) +
    theme(axis.text.x = element_text(angle = 320, size = 8, hjust = 0),
          legend.position = 'top')
  ggsave(filename = paste("figures/jafs", x[1], x[2], "prop_diff_virc.pdf", sep = "_"), plot = pdiff_c, width = 8, height = 8)
})

# CZA_INDEL <- join2_daf(dd1 = daf_ecot$WAVE, dd2 = daf_ecot$CRAB, colnm = "DAF_bin", grps = c("CZA", "INDEL"))
# CZA_INDEL <- zv_dt
# sum(unique(CZA_INDEL$count_prop_pop1))
# sum(unique(CZA_INDEL$count_pop1))
# sum(unique(daf_sum[daf_sum$grp1==zone & daf_sum$ECOT=='WAVE' & daf_sum$grp2==vtype, 'count']))
# sum(unique(CZA_INDEL$count_pop2))
# sum(unique(daf_sum[daf_sum$grp1==zone & daf_sum$ECOT=='CRAB' & daf_sum$grp2==vtype, 'count']))
# sum(unique(CZA_INDEL$count_prop_pop2))
# rm(CZA_INDEL)





