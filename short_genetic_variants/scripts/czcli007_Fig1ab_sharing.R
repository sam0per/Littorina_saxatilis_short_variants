rm (list=ls())
source(file = "/Users/samuelperini/Documents/research/projects/3.indels/Littorina_saxatilis/short_genetic_variants/scripts/colour_text_hadley.R")
# setwd("/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115/")

.packages = c("optparse", "dplyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-v", "--variant"), type="character", default=NULL,
              help="SNP or INDEL", metavar="character"),
  make_option(c("-i", "--inversion"), type="logical", default=FALSE,
              help="TRUE or FALSE [default: %default]", metavar="logical"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$variant)){
  print_help(opt_parser)
  stop("The type of the variant must be selected, choose between SNP or INDEL.\n", call.=FALSE)
}

################################################################################################################
##### INPUT - READ EITHER DATA WITH OR WITHOUT INV #############################################################
################################################################################################################
# setwd("Anja/Anja_results/20200115/CZCLI006_comp/")
liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")
# vartype = "SNP"
vartype = opt$variant
# YNinv = "NoInv"
if (opt$inversion) {
  YNinv = ""
} else {
  YNinv = "NoInv"
}

# Get cline fits
# ANG_right = CZA_right
ANG_right = read.table(paste0("CZCLI006_comp/CZCLI006_ANG_right", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZA_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZA_left", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZA_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZA_right", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZB_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZB_left", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZB_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZB_right", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZD_left = read.table(paste0("CZCLI006_comp/CZCLI006_CZD_left", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)
CZD_right = read.table(paste0("CZCLI006_comp/CZCLI006_CZD_right", YNinv, "_", vartype, ".txt"), header=T, stringsAsFactors=F)



################################################################################################################
##### TOTAL NUMBER OF SNPs ANALYSED ############################################################################
################################################################################################################

# dim(ANG_right)[1]
cat(colourise(paste("Total number of", vartype, YNinv, sep = " "), "blue"), "\n")
dim(CZA_left)[1]



################################################################################################################
##### PROPORTION OF SNPs WITH SIGNIFICANT CLINES ###############################################################
################################################################################################################
cat("\n")
cat(colourise(paste("Proportion of", vartype, YNinv, "with significant clines.", sep = " "), "blue"), "\n")
# length(ANG_right$cp[ANG_right$Type=="Cline"]) / length(ANG_right$cp)
cat(colourise("CZA left:", "light blue"), "\n")
length(CZA_left$cp[CZA_left$Type=="Cline"]) / length(CZA_left$cp)
cat(colourise("CZA right:", "light blue"), "\n")
length(CZA_right$cp[CZA_right$Type=="Cline"]) / length(CZA_right$cp)
cat(colourise("CZB left:", "light blue"), "\n")
length(CZB_left$cp[CZB_left$Type=="Cline"]) / length(CZB_left$cp)
cat(colourise("CZB right:", "light blue"), "\n")
length(CZB_right$cp[CZB_right$Type=="Cline"]) / length(CZB_right$cp)
cat(colourise("CZD left:", "light blue"), "\n")
length(CZD_left$cp[CZD_left$Type=="Cline"]) / length(CZD_left$cp)
cat(colourise("CZD right:", "light blue"), "\n")
length(CZD_right$cp[CZD_right$Type=="Cline"]) / length(CZD_right$cp)
cat("\n")


################################################################################################################
##### PROP OF OUTLIERS THAT ARE SHARED #########################################################################
################################################################################################################
cat(colourise(paste("Proportions of", vartype, YNinv, "outliers that are shared.", sep = " "), "blue"), "\n")
cat(colourise("CZA left and right:", "light blue"), "\n")
length(CZA_left$cp[CZA_left$sel==T & CZA_right$sel==T]) / length(CZA_left$cp[CZA_left$sel==T])
cat(colourise("CZB left and right:", "light blue"), "\n")
length(CZB_left$cp[CZB_left$sel==T & CZB_right$sel==T]) / length(CZB_left$cp[CZB_left$sel==T])
cat(colourise("CZD left and right:", "light blue"), "\n")
length(CZD_left$cp[CZD_left$sel==T & CZD_right$sel==T]) / length(CZD_left$cp[CZD_left$sel==T])

# mean(c(
#   length(ANG_right$cp[ANG_right$sel==T & CZA_left$sel==T]) / length(ANG_right$cp[ANG_right$sel==T]),
#   length(ANG_right$cp[ANG_right$sel==T & CZA_right$sel==T]) / length(ANG_right$cp[ANG_right$sel==T])))
# 
# mean(c(
#   length(ANG_right$cp[ANG_right$sel==T & CZB_left$sel==T]) / length(ANG_right$cp[ANG_right$sel==T]),
#   length(ANG_right$cp[ANG_right$sel==T & CZB_right$sel==T]) / length(ANG_right$cp[ANG_right$sel==T])))
# 
# mean(c(
#   length(ANG_right$cp[ANG_right$sel==T & CZD_left$sel==T]) / length(ANG_right$cp[ANG_right$sel==T]),
#   length(ANG_right$cp[ANG_right$sel==T & CZD_right$sel==T]) / length(ANG_right$cp[ANG_right$sel==T])))

cat(colourise("CZA and CZB:", "light blue"), "\n")
mean(c(
  length(CZA_left$cp[CZA_left$sel==T & CZB_left$sel==T]) / length(CZA_left$cp[CZA_left$sel==T]),
  length(CZA_left$cp[CZA_left$sel==T & CZB_right$sel==T]) / length(CZA_left$cp[CZA_left$sel==T]),
  length(CZA_left$cp[CZA_right$sel==T & CZB_left$sel==T]) / length(CZA_right$cp[CZA_right$sel==T]),
  length(CZA_left$cp[CZA_right$sel==T & CZB_right$sel==T]) / length(CZA_right$cp[CZA_right$sel==T])))
cat(colourise("CZA and CZD:", "light blue"), "\n")
mean(c(
  length(CZA_left$cp[CZA_left$sel==T & CZD_left$sel==T]) / length(CZA_left$cp[CZA_left$sel==T]),
  length(CZA_left$cp[CZA_left$sel==T & CZD_right$sel==T]) / length(CZA_left$cp[CZA_left$sel==T]),
  length(CZA_left$cp[CZA_right$sel==T & CZD_left$sel==T]) / length(CZA_right$cp[CZA_right$sel==T]),
  length(CZA_left$cp[CZA_right$sel==T & CZD_right$sel==T]) / length(CZA_right$cp[CZA_right$sel==T])))
cat(colourise("CZB and CZD:", "light blue"), "\n")
mean(c(
  length(CZB_left$cp[CZB_left$sel==T & CZD_left$sel==T]) / length(CZB_left$cp[CZB_left$sel==T]),
  length(CZB_left$cp[CZB_left$sel==T & CZD_right$sel==T]) / length(CZB_left$cp[CZB_left$sel==T]),
  length(CZB_left$cp[CZB_right$sel==T & CZD_left$sel==T]) / length(CZB_right$cp[CZB_right$sel==T]),
  length(CZB_left$cp[CZB_right$sel==T & CZD_right$sel==T]) / length(CZB_right$cp[CZB_right$sel==T])))
cat("\n")


################################################################################################################
##### SHARING WITH ANY ZONE ####################################################################################
################################################################################################################
???????
# zone1 <- "CZA_left"
cat(colourise(paste("Proportions of", vartype, YNinv, "outliers that are shared with any zone.", sep = " "), "blue"), "\n")
for (zone1 in liste){
  out = get(zone1)
  out = out$cp[out$sel==T]
  others = liste[liste!=zone1]
  out_others = c()
  for (zone2 in others){
    out2 = get(zone2)
    out2 = out2$cp[out2$sel==T]
    out_others = append(out_others, out2)
  }
  print(paste(zone1, length(out[out %in% out_others == F]) / length(out)))
}
cat("\n")


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
outliers$ANG_right = NULL

outliers$all = rowSums(outliers[, 5:10]) == 6
outliers$any = rowSums(outliers[, 5:10]) >= 1
outliers$count = rowSums(outliers[, 5:10])
# sample_n(outliers[outliers$count>1, ], 30)

invisible(lapply(1:6, function(x) {
  len_out = length(outliers$cp[outliers$count==x])
  cat(colourise(paste("Number of", vartype, YNinv, "outliers found in", x, "hybrid zone(s):", len_out, sep = " "),
                "light blue"), "\n")
}))
# outliers[outliers$count==6, ]
# length(outliers$cp[outliers$count==1])
# length(outliers$cp[outliers$count==2])
# length(outliers$cp[outliers$count==3])
# length(outliers$cp[outliers$count==4])
# length(outliers$cp[outliers$count==5])
# length(outliers$cp[outliers$count==6])
# length(outliers$cp[outliers$count==7])
# 
# outliers[outliers$count==1 & outliers$invRui!=F, ]
cat("\n")
invisible(lapply(1:6, function(x) {
  prop_out = length(outliers$cp[outliers$count==x & outliers$invRui!=F]) / length(outliers$cp[outliers$count==x])
  cat(colourise(paste("Proportion of", vartype, YNinv, "outliers in inversions found in", x, "hybrid zone(s):", prop_out, sep = " "),
                "light blue"), "\n")
}))
# length(outliers$cp[outliers$count==1 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==1])
# length(outliers$cp[outliers$count==2 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==2])
# length(outliers$cp[outliers$count==3 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==3])
# length(outliers$cp[outliers$count==4 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==4])
# length(outliers$cp[outliers$count==5 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==5])
# length(outliers$cp[outliers$count==6 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==6])
# length(outliers$cp[outliers$count==7 & outliers$invRui!=F]) / length(outliers$cp[outliers$count==7])

# if (vartype == "indel") {
#   outliers[outliers$count==6, ]
# }
