rm (list=ls())
setwd("/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115")

.packages = c("optparse")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-v", "--variant"), type="character", default=NULL,
              help="snp or indel", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$variant)){
  print_help(opt_parser)
  stop("The type of the variant must be selected, choose between snp or indel.\n", call.=FALSE)
}
################################################################################################################
##### INPUT ####################################################################################################
# vartype = "snp"
vartype = opt$variant
# Get cline fits
# ANG_right = read.table("CZCLI005_ANG_right_means.txt", header=T, stringsAsFactors=F)
CZA_left = read.table(paste0("CZCLI005_means/CZCLI005_CZA_left_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
CZA_right = read.table(paste0("CZCLI005_means/CZCLI005_CZA_right_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
CZB_left = read.table(paste0("CZCLI005_means/CZCLI005_CZB_left_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
CZB_right = read.table(paste0("CZCLI005_means/CZCLI005_CZB_right_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
CZD_left = read.table(paste0("CZCLI005_means/CZCLI005_CZD_left_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
CZD_right = read.table(paste0("CZCLI005_means/CZCLI005_CZD_right_", vartype, "_means.txt"), header=T, stringsAsFactors=F)
ANG_right = CZA_right

# Keep only SNPs with data in all zones
use = intersect(intersect(intersect(intersect(intersect(intersect(
  ANG_right$cp,
  CZA_left$cp),
  CZA_right$cp),
  CZB_left$cp),
  CZB_right$cp),
  CZD_left$cp),
  CZD_right$cp)

length(use)

ANG_right = ANG_right[ANG_right$cp %in% use, ]
ANG_right = ANG_right[order(ANG_right$cp),]
CZA_left = CZA_left[CZA_left$cp %in% use, ]
CZA_left = CZA_left[order(CZA_left$cp),]
CZA_right = CZA_right[CZA_right$cp %in% use, ]
CZA_right = CZA_right[order(CZA_right$cp),]
CZB_left = CZB_left[CZB_left$cp %in% use, ]
CZB_left = CZB_left[order(CZB_left$cp),]
CZB_right = CZB_right[CZB_right$cp %in% use, ]
CZB_right = CZB_right[order(CZB_right$cp),]
CZD_left = CZD_left[CZD_left$cp %in% use, ]
CZD_left = CZD_left[order(CZD_left$cp),]
CZD_right = CZD_right[CZD_right$cp %in% use, ]
CZD_right = CZD_right[order(CZD_right$cp),]


# Function to find nnSNPs
getOutl = function(x, y){
  thresh = sort(x$Var.Ex, decreasing=T, na.last=T)[length(x$cp)/100]
  x$sel = F
  x$sel[x$Var.Ex>=thresh & is.na(x$Var.Ex)==F] = T
  return(x)
}


# Identify outliers
ANG_right = getOutl(ANG_right, "ANG")
CZA_left = getOutl(CZA_left, "CZA_left")
CZA_right = getOutl(CZA_right, "CZA_right")
CZB_left = getOutl(CZB_left, "CZB_left")
CZB_right = getOutl(CZB_right, "CZB_right")
CZD_left = getOutl(CZD_left, "CZD_left")
CZD_right = getOutl(CZD_right, "CZD_right")


# Write out
liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")

for (zone in liste){
  dat = get(zone)
  write.table(dat, paste("CZCLI006_comp/CZCLI006_", zone, "_", vartype, ".txt", sep=""), col.names=T, quote=F, row.names=F, append=F)
}



################################################################################################################
##### MAKE FILE WITHOUT RUI'S INVERSIONS #######################################################################

# Remove inversions and LG12
ANG_rightNoInv = ANG_right[ANG_right$LG!=12 & ANG_right$invRui==F, ]
CZA_leftNoInv = CZA_left[CZA_left$LG!=12 & CZA_left$invRui==F, ]
CZA_rightNoInv = CZA_right[CZA_right$LG!=12 & CZA_right$invRui==F, ]
CZB_leftNoInv = CZB_left[CZB_left$LG!=12 & CZB_left$invRui==F, ]
CZB_rightNoInv = CZB_right[CZB_right$LG!=12 & CZB_right$invRui==F, ]
CZD_leftNoInv = CZD_left[CZD_left$LG!=12 & CZD_left$invRui==F, ]
CZD_rightNoInv = CZD_right[CZD_right$LG!=12 & CZD_right$invRui==F, ]


# Find outliers when inversions are excluded
ANG_rightNoInv = getOutl(ANG_rightNoInv, "ANG")
CZA_leftNoInv = getOutl(CZA_leftNoInv, "CZA_left")
CZA_rightNoInv = getOutl(CZA_rightNoInv, "CZA_right")
CZB_leftNoInv = getOutl(CZB_leftNoInv, "CZB_left")
CZB_rightNoInv = getOutl(CZB_rightNoInv, "CZB_right")
CZD_leftNoInv = getOutl(CZD_leftNoInv, "CZD_left")
CZD_rightNoInv = getOutl(CZD_rightNoInv, "CZD_right")


# Write out
liste = c("ANG_rightNoInv", "CZA_leftNoInv", "CZA_rightNoInv", "CZB_leftNoInv", "CZB_rightNoInv",
          "CZD_leftNoInv", "CZD_rightNoInv")

for (zone in liste){
  dat = get(zone)
  write.table(dat, paste("CZCLI006_comp/CZCLI006_", zone, "_", vartype, ".txt", sep=""), col.names=T, quote=F, row.names=F, append=F)
}
