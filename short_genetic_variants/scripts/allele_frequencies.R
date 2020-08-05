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
ancs <- unique(ancs)
head(ancs)
# ancs[ancs$cp=='Contig23485_2700', ]

ancs$NE_F1_141_Lc <- ifelse(test = ancs$NE_F1_141_Lc=="-het", yes = NA, no = as.character(ancs$NE_F1_141_Lc))
ancs$W_com_01_Lc <- ifelse(test = ancs$W_com_01_Lc=="-het", yes = NA, no = as.character(ancs$W_com_01_Lc))

(frq_fl <- list.files(path = "summary/allele_freq", pattern = ".txt", full.names = TRUE))
frqs <- lapply(frq_fl, read.table, header = TRUE)
lapply(frqs, head)

frqs_dt <- as.data.frame(rbindlist(lapply(seq_along(frq_fl), function(x) {
  
  island <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][2]
  ecot <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][3]
  if (ecot == 'WAVE') {
    side <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][4]
    ecot <- paste(ecot, side, sep = "_")
    vtype <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][5]
  } else {
    vtype <- strsplit(file_path_sans_ext(basename(frq_fl[[x]])), split = "_")[[1]][4]
  }
  
  one_frq <- separate(frqs[[x]], col = REF_FREQ, into = c("REF", "RFREQ"), sep = ":")
  one_frq <- separate(one_frq, col = ALT_FREQ, into = c("ALT", "AFREQ"), sep = ":")
  
  one_frq <- mutate(one_frq, ZONE = as.character(island), ECOT = as.character(ecot),
                    VTYPE = as.character(strsplit(x = vtype, split = ".", fixed = TRUE)[[1]][1]),
                    cp = paste(one_frq[, "CHROM"], one_frq[, "POS"], sep = "_"))
  one_frq <- unique(one_frq)
  return(one_frq)
})))
head(frqs_dt)
table(frqs_dt$ECOT)

franc_dt <- merge(x = ancs, y = frqs_dt, by = c("cp", "ZONE", "VTYPE"))
franc_dt <- franc_dt[rowSums(apply(X = franc_dt[, c("NE_F1_141_Lc", "W_com_01_Lc")], MARGIN = 2,
                                   FUN = function(x) is.na(x))) != 2, ]

# franc_dt[duplicated(franc_dt), ]
head(franc_dt)

rm(list = setdiff(x = ls(), y = "franc_dt"))

table(franc_dt$NE_F1_141_Lc)

# len_pal <- colorRampPalette("Green")
# 
# getVenn <- function(dt, grp1, grp2, vtype, cpal) {
#   set1 <- as.character(dt[dt$NE_F1_141_Lc==grp1 & dt$VTYPE==vtype, "cp"])
#   set2 <- as.character(dt[dt$W_com_01_Lc==grp2 & dt$VTYPE==vtype, "cp"])
#   universe <- sort(unique(c(set1, set2)))
#   Counts <- matrix(0, nrow=length(universe), ncol=2)
#   for (i in 1:length(universe)) {
#     Counts[i,1] <- universe[i] %in% set1
#     Counts[i,2] <- universe[i] %in% set2
#   }
#   colnames(Counts) <- c(paste(vtype, grp1, "NE", sep = "_"), paste(vtype, grp2, "W", sep = "_"))
#   
#   # Specify the colors for the sets
#   # cols <- c("Red", "Green")
#   cols <- brewer.pal(n = 9, name = cpal)[c(4,7)]
#   
#   if (grp1!=grp2) {
#     pdf(file = paste("figures/Venn", vtype, grp1, "NE", grp2, "W", "Lcomp.pdf", sep = "_"), width = 10)
#     vennDiagram(vennCounts(Counts), circle.col=cols,  cex=c(1,2,1))
#     dev.off()
#   } else {
#     pdf(file = paste("figures/Venn", vtype, grp1, "Lcomp.pdf", sep = "_"), width = 10)
#     vennDiagram(vennCounts(Counts), circle.col=cols,  cex=c(1,2,1))
#     dev.off()
#   }
# }
# getVenn(dt = franc_dt, grp1 = "ref_anc", grp2 = "ref_anc", vtype = "INDEL", cpal = "Greens")
# getVenn(dt = franc_dt, grp1 = "ref_anc", grp2 = "alt_anc", vtype = "INDEL", cpal = "Greens")
# getVenn(dt = franc_dt, grp1 = "alt_anc", grp2 = "ref_anc", vtype = "INDEL", cpal = "Greens")
# getVenn(dt = franc_dt, grp1 = "alt_anc", grp2 = "alt_anc", vtype = "INDEL", cpal = "Greens")
# getVenn(dt = franc_dt, grp = "het", vtype = "INDEL", cpal = "Greens")
# 
# getVenn(dt = franc_dt, grp = "ref_anc", vtype = "SNP", cpal = "Greys")
# getVenn(dt = franc_dt, grp = "alt_anc", vtype = "SNP", cpal = "Greys")
# getVenn(dt = franc_dt, grp = "het", vtype = "SNP", cpal = "Greys")

# head(franc_dt)

outg <- which(colnames(franc_dt) %in% c("NE_F1_141_Lc", "W_com_01_Lc"))
# combn(x = unique(franc_dt[, outg[1]]), m = 2)[,4]
franc_dt$NE_W_Lcomp <- paste(franc_dt$NE_F1_141_Lc, franc_dt$W_com_01_Lc, sep = ":")
anc_tb <- data.frame(table(franc_dt$VTYPE, franc_dt$NE_W_Lcomp))
anc_tb <- separate(data = anc_tb, col = Var2, into = c("NE", "W"), sep = ":")
anc_tb <- split(x = anc_tb, f = anc_tb$Var1)
# View(anc_tb$INDEL)
# View(anc_tb$SNP)

# franc_dt <- franc_dt[franc_dt$NE_W_Lcomp!="alt_anc:ref_anc" & franc_dt$NE_W_Lcomp!="ref_anc:alt_anc" &
#                        franc_dt$NE_W_Lcomp!="het:het" & franc_dt$NE_W_Lcomp!="het:NA" &
#                        franc_dt$NE_W_Lcomp!="NA:het", ]
franc_dt <- franc_dt[franc_dt$NE_W_Lcomp=="alt_anc:alt_anc" | franc_dt$NE_W_Lcomp=="ref_anc:ref_anc", ]
data.frame(table(franc_dt$NE_W_Lcomp))
franc_dt$DAF <- NA
franc_dt$DAF <- ifelse(test = grepl(pattern = "alt", x = franc_dt$NE_W_Lcomp), yes = franc_dt$RFREQ, no = franc_dt$AFREQ)
sample_n(tbl = franc_dt, size = 30)

franc_dt$DAF <- as.numeric(franc_dt$DAF)
franc_dt$cp <- as.character(franc_dt$cp)
str(franc_dt)

which.max(aggregate(franc_dt$DAF, by = list(cp = franc_dt$cp), sum)$x)
aggregate(franc_dt$DAF, by = list(cp = franc_dt$cp), sum)[81329, ]
cpdaf_sum <- aggregate(franc_dt$DAF, by = list(cp = franc_dt$cp), sum)
head(cpdaf_sum)
fixed_v <- as.character(unlist(lapply(X = 0:9, function(x) cpdaf_sum[cpdaf_sum$x==x, 'cp'])))

# franc_dt[franc_dt$cp=='Contig3492_13320', ]
# franc_dt[franc_dt$cp=='Contig2844_2649', ]

table(unique(franc_dt[, c('cp', 'VTYPE')])$VTYPE)
table(unique(franc_dt[franc_dt$cp %in% fixed_v, c('cp', 'VTYPE')])$VTYPE)
franc_dt <- franc_dt[!franc_dt$cp %in% fixed_v, ]
head(franc_dt)

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
head(fai)
colnames(fai) <- c("CHROM", "Length")

franc_dt <- merge(dt, fai, by = 'CHROM')
write.csv(x = franc_dt, file = "results/Lsax_short_var_czs_daf.csv", quote = FALSE, row.names = FALSE)

mapdata = read.table("data/map_v11.txt", header = TRUE, stringsAsFactors = FALSE)
head(mapdata)
mapdata$cp = paste(mapdata$contig, mapdata$pos, sep="_")

# Assign SNPs not on map to closest map position, if within 1000bp
means <- franc_dt
rm(franc_dt)
for (number1 in seq(1,length(means$cp))){
  contig = means[number1, "CHROM"]
  pos = means[number1, "POS"]
  focal = mapdata[mapdata$contig==contig, ]
  closest = focal[abs(focal$pos-pos) == min(abs(focal$pos-pos)), ][1,]
  
  if((abs(closest$pos-pos))<=1000 & (is.na(closest$pos)==F)){
    means[number1, "LG"] = closest$LG
    means[number1, "av"] = closest$av
  }
}
head(means)
# Round map positions
means$av = round(means$av, 1)
table(means$LG)

invRui = read.table("data/20200123/Sweden_inversions_coordinates_2nd_august_2019.csv", sep=",", header=TRUE,
                    stringsAsFactors=FALSE)
invRui$LG = gsub("LG", "", invRui$LG)

getInvRui = function(x){
  x$invRui = F
  for (number1 in 1:length(x$av)){
    av = x[number1, "av"]
    LG = x[number1, "LG"]
    if (LG %in% invRui$LG){
      focal = invRui[invRui$LG==LG,]
      for (number2 in seq(1, length(focal$Cluster))){
        if (av>=(focal[number2, "Start"]-0.1) & av<=(focal[number2, "End"]+0.1)){
          x$invRui[number1] = focal$Cluster[number2]
        }
      }
    }
  }
  return(x)
}
means <- getInvRui(means)
head(means)
table(means$invRui)
write.csv(x = means, file = "results/Lsax_short_var_czs_daf_inv.csv", quote = FALSE, row.names = FALSE)






# contig_bin <- merge(grid_dt, aggr_dt, all.x = TRUE)

vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")
lapply(X = c("INDEL", "SNP"), FUN = function(y) {
  hist_p <- ggplot(data = franc_dt[franc_dt$VTYPE==y, ], aes(x = DAF)) +
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
  ggsave(filename = paste0("figures/daf_hist_prop_filt2_", y, ".pdf"), plot = hist_p, width = 10, height = 5)
})
# franc_dt[franc_dt$DAF>0.75,]

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
  one_frq <- separate(data = one_frq, col = cp, into = c("Contig", "Position"), sep = "_", remove = FALSE)
  return(one_frq)
  
  # hist(maf, breaks = 20, main = paste(vtype, "maf"))
  # abline(v = 0.1, col = "red")
  # hist(maf[maf>0.1], breaks = 20, main = paste(vtype, "maf"))
})))

if (exists("franc_dt")) {
  frqs_dt <- franc_dt
  names(frqs_dt)[names(frqs_dt) == 'CHROM'] <- 'Contig'
  names(frqs_dt)[names(frqs_dt) == 'POS'] <- 'Position'
  frqs_dt$MAF <- frqs_dt$DAF
  frqs_dt$cp <- as.character(frqs_dt$cp)
  frqs_dt$Contig <- as.character(frqs_dt$Contig)
}
head(frqs_dt)
str(frqs_dt)

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# head(fai)
colnames(fai) <- c("CHROM", "Length")
# fai$Contig <- as.character(fai$Contig)

# comp_fai <- merge(comp_freq, fai, by = "Contig")
# frqs_fai <- as.data.frame(merge(frqs_dt, fai, by = "Contig"))
# head(frqs_fai)
# frqs_fai <- frqs_fai[!is.na(frqs_fai$MAF), ]
# frqs_fai$frq_bin <- cut(x = frqs_fai$MAF, breaks = seq(from = 0, to = round(max(frqs_fai$MAF)), by = 0.05),
#                         include.lowest = TRUE)
# table(frqs_fai$frq_bin)
# contig_bin <- merge(expand.grid(VTYPE = unique(frqs_fai$VTYPE),
#                                 Contig = unique(frqs_fai$Contig),
#                                 FREQ = levels(frqs_fai$frq_bin)),
#                     aggregate(x = frqs_fai$MAF, list(VTYPE=frqs_fai$VTYPE,
#                                                      Contig=frqs_fai$Contig,
#                                                      FREQ=frqs_fai$frq_bin),
#                               function(x) c(Count = as.integer(length(x)), Mean_frq = round(mean(x), 3))), all.x = TRUE)

contig_bin <- merge(grid_dt, aggr_dt, all.x = TRUE)
contig_bin[1:20,]
# contig_bin[contig_bin$Contig=="Contig0", ]
# frqs_fai[frqs_fai$Contig=="Contig0", ]
contig_bin <- data.frame(contig_bin[, 1:5], 
                         apply(X = contig_bin$x, MARGIN = 2, FUN = function(x) ifelse(test = is.na(x), yes = 0, no = x)))
contig_bin <- merge(x = contig_bin, y = fai, by = "CHROM")
# contig_bin$Count_Len <- round(contig_bin$Count/contig_bin$Length, 5)

# contig_bin$Proportion <- NA
# table(frqs_dt$VTYPE)
# contig_bin[contig_bin$VTYPE=="INDEL", "Proportion"] <- round(contig_bin[contig_bin$VTYPE=="INDEL", "Count"]/
#                                                                table(frqs_dt$VTYPE)["INDEL"], 5)
# contig_bin[contig_bin$VTYPE=="SNP", "Proportion"] <- round(contig_bin[contig_bin$VTYPE=="SNP", "Count"]/
#                                                              table(frqs_dt$VTYPE)["SNP"], 5)
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

contig_bin[1:30,]
contig_bin <- contig_bin[order(contig_bin$DAF_bin), ]
daf_sum <- aggregate(contig_bin$count, by = list(ZONE=contig_bin$ZONE, ECOT=contig_bin$ECOT,
                                                 VTYPE=contig_bin$VTYPE, DAF_bin=contig_bin$DAF_bin), sum)
head(daf_sum)
colnames(daf_sum) <- c("grp1", "ECOT", "grp2", "DAF_bin", "count")
daf_ecot <- split(daf_sum, f = daf_sum$ECOT)
lapply(daf_ecot, head)
colnames(daf_ecot$WAVE)[5] <- "count_pop1"
colnames(daf_ecot$CRAB)[5] <- "count_pop2"
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
CZA_INDEL <- join2_daf(dd1 = daf_ecot$WAVE, dd2 = daf_ecot$CRAB, colnm = "DAF_bin", grps = c("CZA", "INDEL"))
CZA_INDEL$DAF_bin_pop1 <- relevel(CZA_INDEL$DAF_bin_pop1, "[0,0.05]")
CZA_INDEL$DAF_bin_pop2 <- relevel(CZA_INDEL$DAF_bin_pop2, "[0,0.05]")
str(CZA_INDEL)
head(CZA_INDEL)
CZA_INDEL$sqrt_pop1 <- sqrt(CZA_INDEL$count_pop1)
CZA_INDEL$sqrt_pop2 <- sqrt(CZA_INDEL$count_pop2)
CZA_INDEL$sqrt_sum <- rowSums(x = CZA_INDEL[, 7:8])

ggp <- ggplot(CZA_INDEL, aes(DAF_bin_pop1, DAF_bin_pop2)) +
  geom_tile(aes(fill = sqrt_sum))
ggp

# order(levels(contig_bin$FREQ))
# contig_bin[which.max(contig_bin$Count_Tot), ]

# ggplot(contig_bin, aes(x = Length, y = Proportion, col = VTYPE)) +
#   facet_wrap(~FREQ) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   scale_color_manual(values = c("#666666", "#1B9E77")) +
#   labs(x = "contig length", y = "variant proportion", col = "") +
#   theme(axis.text.x = element_text(angle = 320, hjust = 0))

# freq_split <- split(contig_bin, f = contig_bin$VTYPE)
freq_split <- split(contig_bin, f = contig_bin$ECOT)
# identical(freq_split$SNP$Contig,freq_split$INDEL$Contig)
# identical(freq_split$SNP$FREQ,freq_split$INDEL$FREQ)
# identical(freq_split$SNP$Length,freq_split$INDEL$Length)
identical(freq_split$WAVE$CHROM,freq_split$CRAB$CHROM)
identical(freq_split$WAVE$DAF_bin,freq_split$CRAB$DAF_bin)
identical(freq_split$WAVE$Length,freq_split$CRAB$Length)
identical(freq_split$WAVE$ZONE,freq_split$CRAB$ZONE)
identical(freq_split$WAVE$VTYPE,freq_split$CRAB$VTYPE)
lapply(freq_split, head)

# freq_wide <- data.frame(freq_split$SNP[, c("Contig", "Length", "FREQ")],
#                         INDEL=freq_split$INDEL$Proportion,
#                         INDEL_COUNT=freq_split$INDEL$Count,
#                         SNP=freq_split$SNP$Proportion,
#                         SNP_COUNT=freq_split$SNP$Count)
freq_wide <- data.frame(freq_split$WAVE[, c("CHROM", "Length", "DAF_bin", "ZONE", "VTYPE")],
                        CRAB_COUNT=freq_split$CRAB$count,
                        WAVE_COUNT=freq_split$WAVE$count)


freq_wide[1:20,]

# str(freq_wide)
# order(levels(freq_wide$FREQ))
freq_wide$DAF_bin <- relevel(freq_wide$DAF_bin, "[0,0.05]")

# freq_wide$Len_bin <- cut(freq_wide$Length,
#                          breaks = as.integer(seq(from = 0, to = 500000, length.out = 11)),
#                          include.lowest = TRUE, labels = FALSE)
# freq_wide$Len_bin <- factor(freq_wide$Len_bin, levels = as.character(1:10))
# table(freq_wide$Len_bin)


if (exists("franc_dt")) {
  # write.csv(x = freq_wide, file = "results/variant_prop_DAF_bin.csv", row.names = FALSE)
  write.csv(x = freq_wide, file = "results/ecotype_count_DAF_bin.csv", row.names = FALSE)
} else {
  write.csv(x = freq_wide, file = "results/variant_prop_MAF_bin.csv", row.names = FALSE)
}

# len_pal <- colorRampPalette(c("grey", "black"))
# 
# maf_prop <- ggplot(freq_wide, aes(x = INDEL, y = SNP, col = Len_bin)) +
#   facet_wrap(~FREQ) +
#   geom_point() +
#   geom_abline(slope = 1, linetype = "dashed") +
#   # geom_smooth(method='lm') +
#   scale_color_manual(values = len_pal(10)) +
#   labs(x = "INDEL proportion", y = "SNP proportion", col = "") +
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
# maf_prop

# vtype="INDEL"
ecot_count <- ggplot(freq_wide, aes(x = sqrt(CRAB_COUNT), y = sqrt(WAVE_COUNT))) +
  facet_grid(rows = vars(VTYPE), cols = vars(ZONE)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(from = 0, to = 8, by = 2)) +
  # geom_smooth(method='lm') +
  # scale_color_manual(values = len_pal(10)) +
  # labs(x = "INDEL proportion", y = "SNP proportion", col = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        # legend.position = "top",
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#91bfdb", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
ecot_count

if (exists("franc_dt")) {
  ggsave(filename = "figures/daf_snp_vs_indel_prop.pdf", plot = maf_prop, width = 10, height = 8, dpi = "print")
} else {
  ggsave(filename = "figures/maf_snp_vs_indel_prop.pdf", plot = maf_prop, width = 10, height = 8, dpi = "print")
}
franc_dt[franc_dt$DAF>0.95,]
ggplotly(maf_prop)

freq_wide[freq_wide$FREQ=="(0.1,0.15]", ][which.max(freq_wide[freq_wide$FREQ=="(0.1,0.15]", "INDEL"]), ]
contig_bin[contig_bin$Proportion==0.00053,]

contig_bin[which.max(contig_bin$Proportion), ]
contig_bin[which.max(contig_bin$Length), ]
contig_bin[contig_bin$Length==237503, ]
# ggplot(contig_bin, aes(x = Length, y = Count_Tot, fill = VTYPE)) +
#   facet_wrap(~FREQ) +
#   geom_bar(stat = "identity", position = "dodge")

# contig_bin$Count_Tot <- contig_bin$Count_Tot/100
# summary(lm(formula = Count_Tot~Length, data = contig_bin[contig_bin$VTYPE=="SNP" & contig_bin$FREQ=="[0,0.05]", ]))

# with(data = contig_bin, plot(x = Length, Mean_frq))

freq_bin <- merge(expand.grid(VTYPE = levels(contig_bin$VTYPE),
                              FREQ = levels(contig_bin$FREQ)),
                    aggregate(x = contig_bin$Count_Len, list(VTYPE=contig_bin$VTYPE,
                                                             FREQ=contig_bin$FREQ), sum), all.x = TRUE)
freq_bin
levels(freq_bin$FREQ)
sort(levels(freq_bin$FREQ))
freq_bin$SORT <- rep(order(levels(freq_bin$FREQ)), 2)
freq_bin <- freq_bin[order(freq_bin$SORT), ]

freq_wide <- rbind(INDEL=freq_bin[freq_bin$VTYPE=="INDEL", "x"],
                   SNP=freq_bin[freq_bin$VTYPE=="SNP", "x"])
freq_wide <- cbind(INDEL=freq_bin[freq_bin$VTYPE=="INDEL", "x"],
                   SNP=freq_bin[freq_bin$VTYPE=="SNP", "x"])
cor(freq_wide)

vtype_pal <- data.frame(INDEL="#1B9E77", SNP="#666666")

ggplot(data = freq_bin, aes(x = FREQ, y = x, fill = VTYPE)) +
  geom_histogram(stat = "identity", position = "dodge", col = "black") +
  scale_fill_manual(values = levels(unlist(vtype_pal)))

# display.brewer.all(colorblindFriendly = TRUE)
# display.brewer.pal(n = 8, name = 'Dark2')
# brewer.pal(n = 11, name = 'YlOrRd')
# display.brewer.pal(n = 9, name = 'Greys')
# brewer.pal(n = 9, name = 'Greys')
# display.brewer.pal(n = 9, name = 'Greens')
# brewer.pal(n = 9, name = 'Greens')



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

# cline_type <- "Cline"
# isla <- "CZB"
# side <- "right"
split_vtype <- split(comp_freq[comp_freq$Type==cline_type, ], f = as.factor(comp_freq[comp_freq$Type==cline_type, "VTYPE"]))
split_vtype <- split(comp_freq, f = as.factor(comp_freq[, "VTYPE"]))
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
# sum(split_vtype$INDEL$MAF==grid_dt[1,1])/nrow(split_vtype$INDEL) +
#   sum(split_vtype$SNP$MAF==grid_dt[1,2])/nrow(split_vtype$SNP)
# sum(split_vtype$INDEL$MAF==grid_dt[4,1])/nrow(split_vtype$INDEL) +
#   sum(split_vtype$SNP$MAF==grid_dt[4,2])/nrow(split_vtype$SNP)
# sum(split_vtype$INDEL$MAF==grid_dt[2,1], split_vtype$SNP$MAF==grid_dt[2,2])

grid_dt$count <- apply(X = grid_dt, MARGIN = 1, FUN = function(x) {
  sum(split_vtype$INDEL$MAF==x[1])/nrow(split_vtype$INDEL) + sum(split_vtype$SNP$MAF==x[2])/nrow(split_vtype$SNP)
})
head(grid_dt)
grid_dt[1:15,]
grid_dt$INDEL <- as.factor(grid_dt$INDEL)
grid_dt$SNP <- as.factor(grid_dt$SNP)

# ggp <- ggplot(grid_dt, aes(INDEL, SNP)) +
#   geom_tile(aes(fill = count))
# ggp

# wide_dt <- reshape(grid_dt, idvar = "INDEL", timevar = "SNP", direction = "wide")
wide_dt <- reshape2::dcast(grid_dt, SNP ~ INDEL)
head(wide_dt)
# install.packages('textshape')
library(textshape)
wide_dt <- as.matrix(column_to_rownames(wide_dt, 'SNP'))
wide_dt <- round(wide_dt, 2)
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



