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
# franc_dt <- franc_dt[!franc_dt$cp %in% fixed_v, ]
head(franc_dt)

fai_path <- "/Users/samuelperini/Documents/research/projects/3.indels/data/reference/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta.fai"
fai <- read.table(file = fai_path, header = FALSE, sep = "\t")[, 1:2]
head(fai)
colnames(fai) <- c("CHROM", "Length")

dt <- franc_dt
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