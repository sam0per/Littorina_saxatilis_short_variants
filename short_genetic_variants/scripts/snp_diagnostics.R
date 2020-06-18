# rm (list=ls())
setwd("Anja/Anja_results/20200115/")
setwd("../../../")
getwd()
################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################

liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")
vtype <- "SNP"
tool <- "gatk"
fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
# head(fai)
# Get inversion info
invRui = read.table("/Users/samuelperini/Documents/research/projects/3.indels/data/20200123/Sweden_inversions_coordinates_2nd_august_2019.csv",
                    sep=",", header=T, stringsAsFactors=F)
invRui$LG = gsub("LG", "", invRui$LG)

# Get cline fits
ANG_right = read.table(paste0("CZCLI006_comp/CZCLI006_ANG_right_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZA_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZA_left_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZA_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZA_right_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZB_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZB_left_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZB_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZB_right_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZD_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZD_left_", vtype, ".txt"), header=T, stringsAsFactors=F)
CZD_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZD_right_", vtype, ".txt"), header=T, stringsAsFactors=F)
table(ANG_right$Type)
table(CZB_right$Type)

################################################################################################################
##### MAKE OUTLIER TABLES ######################################################################################
################################################################################################################

outliers = ANG_right[, c("cp", "LG", "av", "invRui")]
outliers$ANG_right = outliers$cp %in% ANG_right$cp[ANG_right$sel==T]
outliers$CZA_left = outliers$cp %in% CZA_left$cp[CZA_left$sel==T]
outliers$CZA_right = outliers$cp %in% CZA_right$cp[CZA_right$sel==T]
outliers$CZB_left = outliers$cp %in% CZB_left$cp[CZB_left$sel==T]
outliers$CZB_right = outliers$cp %in% CZB_right$cp[CZB_right$sel==T]
outliers$CZD_left = outliers$cp %in% CZD_left$cp[CZD_left$sel==T]
outliers$CZD_right = outliers$cp %in% CZD_right$cp[CZD_right$sel==T]

head(outliers)

outliers$all = rowSums(outliers[, 6:11]) == 6
outliers$any = rowSums(outliers[, 6:11]) >= 1
table(outliers$all)
# outliers[outliers$all==TRUE, ]
assign(paste("outliers", tool, vtype, sep="_"), outliers)

# sh_outl <- merge(outliers_samtools_SNP, outliers_gatk_SNP, by = "cp")
# head(sh_outl)
# ggplot(data = sh_outl[sh_outl$all.x==TRUE, ], aes(x = av.x, y = av.y, col = all.y)) +
#   facet_wrap(~LG.x) +
#   geom_abline(slope = 1) +
#   geom_point()
# ggplot(data = sh_outl[sh_outl$all.y==TRUE, ], aes(x = av.x, y = av.y, col = all.x)) +
#   facet_wrap(~LG.x) +
#   geom_abline(slope = 1) +
#   geom_point()
# sum(sh_outl$all.x)
# sum(sh_outl$all.y)
# nrow(sh_outl[sh_outl$all.x==TRUE & sh_outl$all.y==TRUE, ])
# nrow(sh_outl[sh_outl$all.x==FALSE & sh_outl$all.y==FALSE, ])
# sh_outl[sh_outl$all.x==TRUE & sh_outl$all.y==FALSE, ]
# sh_outl[sh_outl$all.x==FALSE & sh_outl$all.y==TRUE, ]
# sh_outl[sh_outl$all.y==TRUE, ]

# Load the library
library(limma)

# Generate example data
outl = "any"
# set1<-letters[1:5]
set1<-outliers_samtools_SNP$cp
# set2<-letters[4:8]
set2<-outliers_gatk_SNP$cp
# set3<-letters[5:9]
set3<-outliers_gatk_SNP[outliers_gatk_SNP[,outl]==TRUE, "cp"]
set4<-outliers_samtools_SNP[outliers_samtools_SNP[,outl]==TRUE, "cp"]
# intersect(set3,set4)
# What are the possible letters in the universe?
universe <- sort(unique(c(set1, set2, set3, set4)))

# Generate a matrix, with the sets in columns and possible letters on rows
Counts <- matrix(0, nrow=length(universe), ncol=4)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% set1
  Counts[i,2] <- universe[i] %in% set2
  Counts[i,3] <- universe[i] %in% set3
  Counts[i,4] <- universe[i] %in% set4
}

# Name the columns with the sample names
# colnames(Counts) <- c("set1","set2")
colnames(Counts) <- c("samtools","gatk",paste0(outl,"_outliers_gatk"),paste0(outl,"_outliers_samtools"))

# Specify the colors for the sets
# cols<-c("Red", "Green")
cols<-c("Red", "Green", "Blue", "Gray")
pdf(file = paste0("figures/snp_diagnostics_", outl, "_outliers.pdf"), width = 10)
vennDiagram(vennCounts(Counts), circle.col=cols)
dev.off()
###################################################################################
##### MARKER DENSITIES ############################################################
###################################################################################
library(tidyr)
library(ggplot2)
outliers <- separate(data = outliers, col = cp, into = c("chr", "pos"), sep = "_")
head(outliers)
colnames(fai) <- c("chr", "len")
head(fai)
outliers <- merge(outliers, fai, by = "chr")
length(unique(outliers$chr))
table(outliers[outliers$any==TRUE, "LG"])

contig_count <- aggregate(x = outliers$len, by = list(outliers$chr), function(x) c(count = length(x), len = unique(x)))
head(contig_count)
# m_count$nN <- round(m_count$x[, 1] / nrow(outliers), 3)
m_count <- data.frame(chr = contig_count$Group.1, contig_count$x, nN = round(contig_count$x[, 1] / nrow(outliers), 3))
# m_count$density <- round(m_count$nN * m_count$x[, 2], 3)
head(m_count)
tail(m_count)

(p_nN <- ggplot(data = m_count, aes(x = len, y = nN)) +
    geom_point(alpha = 0.6, size = 3) +
    labs(x = "contig length", y = paste(vtype, "n / N")) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)))
m_count[order(m_count$nN, decreasing = TRUE), ][1:10,]
ggsave(filename = paste0("figures/prop_nN_all_", vtype, ".pdf"), plot = p_nN)

(p_count <- ggplot(data = m_count, aes(x = len, y = count)) +
    geom_point(alpha = 0.6, size = 3) +
    labs(x = "contig length", y = paste(vtype, "count")) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)))
m_count[order(m_count$count, decreasing = TRUE), ][1:10,]
ggsave(filename = paste0("figures/count_all_", vtype, ".pdf"), plot = p_count)

assign(paste("nN_samtools", vtype, sep="_"), m_count)

rm(list = setdiff(ls(), ls()[grepl(pattern = "samtools", x = ls())]))
# rm(list = setdiff(ls(), "INDEL_nN"))

shared_chr <- data.frame(chr = intersect(SNP_nN$chr, INDEL_nN$chr))
head(shared_chr)
shared <- merge(shared_chr, INDEL_nN, by = "chr")
head(shared)
shared <- merge(shared, SNP_nN, by = "chr")
shared$sum <- shared$count.x + shared$count.y

(p2_nN <- ggplot(data = shared, aes(x = nN.x, y = nN.y, col = len.x)) +
    geom_abline(slope = 1) +
    # facet_wrap(~shared) +
    geom_point(aes(size = sum)) +
    labs(x = "INDEL n / N", y = "SNP n / N", col = "contig length") +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)))
ggsave(filename = paste0("figures/prop_nN_all_INDEL_SNP.pdf"), plot = p2_nN)

(p2_count <- ggplot(data = shared, aes(x = count.x, y = count.y, col = len.x)) +
    geom_abline(slope = 1) +
    # facet_wrap(~shared) +
    geom_point(aes(size = sum)) +
    labs(x = "INDEL count", y = "SNP count", col = "contig length") +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 10)))
ggsave(filename = paste0("figures/count_all_INDEL_SNP.pdf"), plot = p2_count)