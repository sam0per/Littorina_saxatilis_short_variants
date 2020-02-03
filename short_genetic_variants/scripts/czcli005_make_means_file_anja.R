rm(list=ls())
# setwd("Anja/Anja_results/20200115/")
.packages = c("data.table", "Rmisc", "optparse", "dplyr")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="DIRECTORY with input data", metavar="character"),
  make_option(c("-z", "--zone"), type="character", default=NULL,
              help="Name of the zone found in the input filenames [example: CZA_right]", metavar="character"),
  make_option(c("-v", "--variant"), type="character", default=NULL,
              help="VARIANT TYPE, either SNP or INDEL", metavar="character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file", metavar = "character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Compute summary statistics on the cline analysis output keeping zones and variant types separated.",
                          epilogue = "Example: Rscript scripts/czcli005_make_means_file_anja.R -d Anja/Anja_results/20200115/ -z CZA_right -v INDEL -o Anja/Anja_results/20200115/CZCLI005_CZA_right_means.txt")
opt = parse_args(opt_parser)

if (is.null(opt$dir) | is.null(opt$zone) | is.null(opt$variant) | is.null(opt$output)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

setwd(opt$dir)
# zone = "CZD_right"
zone = opt$zone

################################################################################################################
##### INPUT & OVERVIEW #########################################################################################

results = read.table(paste("CZCLI004_", zone, "_clines.txt", sep=""), header=T, stringsAsFactors=F,sep=' ')
results$cp = paste(results$Contig, results$Position, sep="_")

results$Wave[results$Wave=="a1"] = 1
results$Wave[results$Wave=="a2"] = 2
results$Wave = as.numeric(results$Wave)

# total number of SNPs:
pre_res = length(unique(results$cp))
cat("Total number of variants:", pre_res, "\n")
# head(results)


################################################################################################################
##### FILTERING OF RESULTS #####################################################################################

# keep only Cline and NC SNPs
results = results[(results$Type %in% c("Dup>HWE", "SL", "C_Peak", "Stuck"))==F, ]

# dim(results)
cat("Number of variants after filtering:", length(unique(results$cp)),
    "\nProportion is:", length(unique(results$cp))/pre_res, "and it should be ~1/3rd of the previous\n")
# should be ~1/3rd of the previous
# length(unique(results$cp))/pre_res

# vartype = "SNP"
vartype = opt$variant
id_csv = read.csv(file = paste0("CZCLI01_", substr(x = zone, start = 1, stop = 3), ".filt2.vcf.csv"))
# head(id_csv)
colnames(id_csv) = c("Contig", "Position", "REF", "ALT")
len_var = function(n_ref, n_alt, o_dt) {
  c_ref = as.character(n_ref)
  c_alt = as.character(n_alt)
  l_var = nchar(c_ref) - nchar(c_alt)
  if (abs(l_var)==0) {
    vtype = "SNP"
  } else {
    vtype = "INDEL"
  }
  return(vtype)
}
# len_var(n_ref = "AC", n_alt = "CCAA")
id_csv[, 'VTYPE'] = apply(X = id_csv[, c('REF', 'ALT')], MARGIN = 1,
                          FUN = function(x) len_var(n_ref = x[1], n_alt = x[2]))
# table(id_csv$VTYPE)
# head(results)
res_vtype = merge(x = id_csv, y = results, by = c('Contig', 'Position'))
# nrow(res_vtype)
# nrow(results)
# table(res_vtype$VTYPE)
# table(results$Type)
cat("\nTotal number of variants by type:\n", table(res_vtype$VTYPE),
    "\nDownstream analysis is going to be run on", vartype, "\n")

results = res_vtype[res_vtype$VTYPE==vartype, ]
# table(results$VTYPE)
# make table which for each SNP says how many of the replicates are of each type
tab = as.data.frame.matrix(table(results[, c("cp", "Type")]))
# head(tab)
tab[tab$Cline >= 2, "cat"] = "Cline"
tab[tab$Cline < 2, "cat"] = "NoCline"



################################################################################################################
##### CLINAL SNPS: AVERAGE ACROSS REPLICATES ETC. ##############################################################

# only clinal SNPs, and for those only the clinal replicates
Clines_cp = row.names(tab[tab$cat == "Cline" & (is.na(tab$cat)==F), ])
Clines = results[results$cp %in% Clines_cp, ]
Clines = Clines[Clines$Type == "Cline", ]

# get everything on normal scale and calculate some additional parameters
Clines$p_crab = exp(Clines$p_crab) / (1+exp(Clines$p_crab))
Clines$p_wave_right = exp(Clines$p_wave_right) / (1+exp(Clines$p_wave_right))
Clines$p_diff_right = Clines$p_wave_right - Clines$p_crab
Clines$Width_right = exp(Clines$Width_right)
Clines$slope_right = Clines$p_diff_right/Clines$Width_right

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
No.Clines = results[results$cp %in% No.Clines_cp, ]
#No.Clines$p_crab = exp(No.Clines$p_crab) / (1+exp(No.Clines$p_crab))
#No.Clines$p_wave_right = exp(No.Clines$p_wave_right) / (1+exp(No.Clines$p_wave_right))
No.Clines[, 5:24] = NA
No.Clines$p_diff_right = NA
No.Clines$slope_right = NA
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

# Use only columns of interest
use_names = c("cp","Contig","Position","Type","Wave","Centre_right","Width_right","slope_right",
              "p_crab","p_wave_right","p_diff_right","Var.Ex","Fst_right")   
means = means[, use_names]
names(means) = c("cp","Contig","Position","Type","Wave","Centre","Width","slope",
                 "p_crab","p_wave","p_diff","Var.Ex","Fst") 
# head(means)

################################################################################################################
##### ADD OTHER INFO ###########################################################################################

# Add map data #################################################################################################
# Get map
mapdata = read.table("../../../data/map_v11.txt", header=T, stringsAsFactors = F)
# head(mapdata)
mapdata$cp = paste(mapdata$contig, mapdata$pos, sep="_")

# Assign SNPs not on map to closest map position, if within 1000bp
for (number1 in seq(1,length(means$cp))){
  contig = means[number1, "Contig"]
  pos = means[number1, "Position"]
  focal = mapdata[mapdata$contig==contig, ]
  closest = focal[abs(focal$pos-pos) == min(abs(focal$pos-pos)), ][1,]
  
  if((abs(closest$pos-pos))<=1000 & (is.na(closest$pos)==F)){
    means[number1, "LG"] = closest$LG
    means[number1, "av"] = closest$av
  }
  # print(number1)
}

# Round map positions
means$av = round(means$av, 1)


# Add unswitched allele frequencies. The frequencies are from the fitting, but the frequency ###################
# is not of the "Wave" allele but of the reference allele from the assembly.
means$ref_p_wave[means$Wave==1&is.na(means$Wave)==F] = means$p_wave[means$Wave==1&is.na(means$Wave)==F]
means$ref_p_wave[means$Wave==2&is.na(means$Wave)==F] = 1-means$p_wave[means$Wave==2&is.na(means$Wave)==F]
means$ref_p_crab[means$Wave==1&is.na(means$Wave)==F] = means$p_crab[means$Wave==1&is.na(means$Wave)==F]
means$ref_p_crab[means$Wave==2&is.na(means$Wave)==F] = 1-means$p_crab[means$Wave==2&is.na(means$Wave)==F]


# Add end frequencies based on raw genotypes ###################################################################
# W_file = list.files(pattern = paste("CZ000_", zone, "_W_endfreqs.txt", sep=""))
# W = read.table(W_file, header=T, stringsAsFactors=F)
# C_file = list.files(pattern = paste("CZ000_", zone, "_C_endfreqs.txt", sep=""))
# C = read.table(C_file, header=T, stringsAsFactors=F)
# 
# zoneRaw = merge(W, C, by="cp")
# names(zoneRaw) = c("cp", "allele1.W", "allele2.W", "allele1.C", "allele2.C")
# 
# means = merge(means, zoneRaw[, c("cp", "allele1.W", "allele1.C")], by="cp", all.x=T, all.y=F)
# means$p_waveRaw = means$allele1.W
# means$p_crabRaw = means$allele1.C
# means$allele1.C = means$allele1.W = NULL


# Add info about Rui's inversions ##############################################################################
invRui = read.table("../../../data/20200123/Sweden_inversions_coordinates_2nd_august_2019.csv", sep=",", header=T,
                    stringsAsFactors=F)
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

means = getInvRui(means)



################################################################################################################
##### FILTERING ################################################################################################

# ~~~~ EXCLUDE UNMAPPED SNPs!
means = means[is.na(means$LG)==F, ]

# Remove SNPs where fitting is somehow messed up (should be very few, otherwise there is a problem!)
means = means[is.na(means$Fst)==T | means$Fst<=1, ]
means = means[is.na(means$Var.Ex)==T | means$Var.Ex>=0, ]
means = means[is.na(means$p_diff)==T | means$p_diff>=0, ]


# plot(means$p_crabRaw, means$p_crab)
# plot(means$p_waveRaw, means$p_wave)



################################################################################################################
##### OUTPUT ###################################################################################################

# write out
# write.table(means, paste("CZCLI005_", zone, "_means.txt", sep=""),
#             append = F, quote=F, col.names=T, row.names=T)
write.table(means, opt$output,
            append = F, quote=F, col.names=T, row.names=T)
cat("MISSION COMPLETE!\nOutput file is:", opt$output, "\n")
