rm (list=ls())

################################################################################################################
##### INPUT ####################################################################################################
args = commandArgs(trailingOnly=TRUE)

# zone = "CZA"
zone = args[2]

# vartype = "indels"
vartype = args[3]

# side = "left"
side = args[4]

# means_out = paste0("clines/CZCLI005_", vartype, "_", zone, "_", side, "_means.txt")
means_out = args[5]

################################################################################################################
##### INPUT & OVERVIEW #########################################################################################

# results = read.table(paste0("clines/CZCLI003_", vartype, "_", zone, ".txt"), header=T, stringsAsFactors=F)
results = read.table(args[1], header=T, stringsAsFactors=F)
results = results[(is.na(results$Side))==F, ]

# get results only for the focal contact
results_side = results[results$Side == side, ]
results_side$cp = paste(results_side$Contig, results_side$Position, sep="_")


results_side$Wave[results_side$Wave=="a1"] = 1
results_side$Wave[results_side$Wave=="a2"] = 2
results_side$Wave = as.numeric(results_side$Wave)

# total number of SNPs:
length(unique(results_side$cp))



################################################################################################################
##### FILTERING OF RESULTS #####################################################################################

# keep only Cline and NC SNPs
results_side = results_side[(results_side$Type %in% c("Dup>HWE", "SL", "C_Peak", "Stuck"))==F, ]
unique(results_side$Type)

dim(results_side)
length(unique(results_side$cp)) # should be ~1/3rd of the previous


# make table which for each SNP says how many of the replicates are of each type
tab = as.data.frame.matrix(table(results_side[, c("cp", "Type")]))
tab[tab$Cline >= 2, "cat"] = "Cline"
tab[tab$Cline < 2, "cat"] = "NoCline"



################################################################################################################
##### CLINAL SNPS: AVERAGE ACROSS REPLICATES ETC. ##############################################################

# only clinal SNPs, and for those only the clinal replicates
Clines_cp = row.names(tab[tab$cat == "Cline" & (is.na(tab$cat)==F), ])
Clines = results_side[results_side$cp %in% Clines_cp, ]
Clines = Clines[Clines$Type == "Cline", ]

# get everything on normal scale and calculate some additional parameters
Clines$p_crab = exp(Clines$p_crab) / (1+exp(Clines$p_crab))
Clines$p_wave_left = exp(Clines$p_wave_left) / (1+exp(Clines$p_wave_left))
Clines$p_wave_right = exp(Clines$p_wave_right) / (1+exp(Clines$p_wave_right))
Clines$p_diff_left = Clines$p_wave_left - Clines$p_crab
Clines$p_diff_right = Clines$p_wave_right - Clines$p_crab
Clines$Width_left = exp(Clines$Width_left)
Clines$Width_right = exp(Clines$Width_right)
Clines$slope_left = Clines$p_diff_left/Clines$Width_left
Clines$slope_right = Clines$p_diff_right/Clines$Width_right

Clines$Hs_left =
  (2*Clines$p_crab*(1-Clines$p_crab) + 2*Clines$p_wave_left*(1-Clines$p_wave_left))/2
Clines$Ht_left = 2 * ((Clines$p_crab + Clines$p_wave_left)/2) * (((1-Clines$p_crab) + (1-Clines$p_wave_left))/2)
Clines$Fst_left = (Clines$Ht_left - Clines$Hs_left) / Clines$Ht_left

Clines$Hs_right =
  (2*Clines$p_crab*(1-Clines$p_crab) + 2*Clines$p_wave_right*(1-Clines$p_wave_right))/2
Clines$Ht_right = 2 * ((Clines$p_crab + Clines$p_wave_right)/2) * (((1-Clines$p_crab) + (1-Clines$p_wave_right))/2)
Clines$Fst_right = (Clines$Ht_right - Clines$Hs_right) / Clines$Ht_right


# get average over replicate cline fittings
Clines_means = aggregate(Clines, by=list(Clines$cp), FUN = function(x) mean(x, na.rm=T))
Clines_means$cp=NULL
names(Clines_means)[1] = "cp"
spl = data.frame(do.call('rbind', strsplit(as.character(Clines_means$cp), "_", fixed=T)), stringsAsFactors=F)
Clines_means$Contig = spl$X1
Clines_means$Position = as.numeric(spl$X2)
Clines_means$Type = "Cline"



################################################################################################################
##### NON-CLINAL SNPS: AVERAGE ACROSS REPLICATES ETC. ##########################################################

No.Clines_cp = row.names(tab[(tab$cat == "NoCline") & (is.na(tab$cat)==F), ])
No.Clines = results_side[results_side$cp %in% No.Clines_cp, ]
No.Clines[, 4:31]=NA
No.Clines$p_diff_left = NA
No.Clines$p_diff_right = NA
No.Clines$slope_left = NA
No.Clines$slope_right = NA
No.Clines$Hs_left = NA
No.Clines$Ht_left = NA
No.Clines$Fst_left = NA
No.Clines$Hs_right = NA
No.Clines$Ht_right = NA
No.Clines$Fst_right = NA


# get average over replicates - this is not really relevant for No.Clines data
# but generates a single copy of each SNP and correct column order
No.Clines_means = aggregate(No.Clines, by=list(No.Clines$cp), FUN = function(x) mean(x, na.rm=T))
No.Clines_means$cp=NULL
names(No.Clines_means)[1] = "cp"
spl = data.frame(do.call('rbind', strsplit(as.character(No.Clines_means$cp), "_", fixed=T)), stringsAsFactors=F)
No.Clines_means$Contig = spl$X1
No.Clines_means$Position = as.numeric(spl$X2)
No.Clines_means$Type = "NoCline"



################################################################################################################
##### PUT CLINAL AND NON-CLINAL SNPS IN SAME DF ################################################################

means = rbind.data.frame(Clines_means, No.Clines_means)
means$Index = NULL



################################################################################################################
##### ADD OTHER INFO ###########################################################################################

# get map data
mapdata = read.table("data/map_v11.txt", header=T, stringsAsFactors = F)
mapdata$cp = paste(mapdata$contig, mapdata$pos, sep="_")

# assign SNPs not on map to closest map position, if within 1000bp
for (number1 in seq(1,length(means$cp))){
  contig = means[number1, "Contig"]
  pos = means[number1, "Position"]
  focal = mapdata[mapdata$contig==contig, ]
  closest = focal[abs(focal$pos-pos) == min(abs(focal$pos-pos)), ][1,]
  
  if((abs(closest$pos-pos))<=1000 & (is.na(closest$pos)==F)){
    means[number1, "LG"] = closest$LG
    means[number1, "av"] = closest$av
  }
  print(number1)
}


# add inversion info (just the three main inversions from ANG)
means$inv = F
means[  (means$LG==6 & means$av<29.5) & (is.na(means$LG)==F)|
          (means$LG==14 & means$av<12.5) & (is.na(means$LG)==F)|
          (means$LG==17 & means$av>46) & (is.na(means$LG)==F), "inv"] = T



################################################################################################################
##### OUTPUT ###################################################################################################

# remove SNPs where fitting is somehow messed up (should be very few, otherwise there is a problem!)
means = means[is.na(means$Fst_left)==T | means$Fst_left<=1, ]
means = means[is.na(means$Fst_right)==T | means$Fst_right<=1, ]
means = means[is.na(means$Var.Ex)==T | means$Var.Ex>=0, ]
means = means[is.na(means$p_diff_left)==T | means$p_diff_left>=0, ]
means = means[is.na(means$p_diff_right)==T | means$p_diff_right>=0, ]

# add Fst calculated based on raw genotypes
# W_side_file = list.files(pattern = paste("CZ000_", zone, "_W_", side, "_endfreqs.txt", sep=""))
# W_side = read.table(W_side_file, header=T, stringsAsFactors=F)
# C_file = list.files(pattern = paste("CZ000_", zone, "_C_endfreqs.txt", sep=""))
# C = read.table(C_file, header=T, stringsAsFactors=F)
# 
# zoneRaw = merge(W_side, C, by="cp")
# names(zoneRaw) = c("cp", "allele1.W", "allele2.W", "allele1.C", "allele2.C")
# zoneRaw$Hs = ((2*zoneRaw$allele1.W*zoneRaw$allele2.W)+
#                 (2*zoneRaw$allele1.C*zoneRaw$allele2.C))/2
# zoneRaw$Ht = 0.5*(zoneRaw$allele1.W+zoneRaw$allele1.C)*(zoneRaw$allele2.W+zoneRaw$allele2.C)
# zoneRaw$FstRaw = (zoneRaw$Ht-zoneRaw$Hs) / zoneRaw$Ht
# 
# means = merge(means, zoneRaw[, c("cp", "FstRaw", "allele1.W", "allele2.W", "allele1.C", "allele2.C")], by="cp", all.x=T, all.y=F)
# means$p_waveRaw = NA
# means$p_crabRaw = NA
# means$p_waveRaw[means$Wave==1&is.na(means$Wave)==F] = means$allele1.W[means$Wave==1&is.na(means$Wave)==F]
# means$p_crabRaw[means$Wave==1&is.na(means$Wave)==F] = means$allele1.C[means$Wave==1&is.na(means$Wave)==F]
# means$p_waveRaw[means$Wave==2&is.na(means$Wave)==F] = means$allele2.W[means$Wave==2&is.na(means$Wave)==F]
# means$p_crabRaw[means$Wave==2&is.na(means$Wave)==F] = means$allele2.C[means$Wave==2&is.na(means$Wave)==F]
# means$p_waveRaw[is.na(means$Wave)] = means$allele1.W[is.na(means$Wave)]
# means$p_crabRaw[is.na(means$Wave)] = means$allele1.C[is.na(means$Wave)]
# 
# plot(means$p_crab, means$p_crabRaw)
# plot(means$p_wave_right, means$p_waveRaw)

# write out
# write.table(means, paste("CZCLI005_", zone, "_", side, "_means.txt", sep=""),
#             append = F, quote=F, col.names=T, row.names=T)
write.table(means, means_out, append = F, quote=F, col.names=T, row.names=T)


