# rm (list=ls())
setwd("Anja/Anja_results/20200115/")
setwd("../../../")
getwd()
################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################

liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")
vtype <- "INDEL"
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

# Or cline fits without inversions
ANG_right = read.table("CZCLI006_ANG_rightNoInv.txt", header=T, stringsAsFactors=F)
CZA_left = read.table("CZCLI006_CZA_leftNoInv.txt", header=T, stringsAsFactors=F)
CZA_right = read.table("CZCLI006_CZA_rightNoInv.txt", header=T, stringsAsFactors=F)
CZB_left = read.table("CZCLI006_CZB_leftNoInv.txt", header=T, stringsAsFactors=F)
CZB_right = read.table("CZCLI006_CZB_rightNoInv.txt", header=T, stringsAsFactors=F)
CZD_left = read.table("CZCLI006_CZD_leftNoInv.txt", header=T, stringsAsFactors=F)
CZD_right = read.table("CZCLI006_CZD_rightNoInv.txt", header=T, stringsAsFactors=F)



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
outliers[outliers$all==TRUE, ]

################################################################################################################
##### TEST FOR ENRICHMENT OF OUTLIERS IN INVERSIONS ############################################################
################################################################################################################

invs = unique(outliers$invRui)
invs = sort(invs[invs!=FALSE])

# Make table of all inversions x all zones
tab = data.frame(matrix(nrow=length(invs), ncol=7))
rownames(tab) = invs
colnames(tab) = liste
for (inv in invs){
  for (zone in liste){
    num_out_inv = length(outliers$cp[outliers$invRui==inv & outliers[,zone]==T])
    num_in_inv =length(outliers$cp[outliers$invRui==inv & outliers[,zone]==F])
    
    # Do permutations
    # Sample as many loci as located in the inversion from collinear regions and ask if they are outliers
    # (in focal zone). This tells us how many outliers are expected if there is no enrichment in invs.
    num_out_perm = c()
    for (num in 1:1000){
      perm = sample(outliers[outliers$invRui==F, zone], num_out_inv+num_in_inv, replace=T)
      num_out_perm = append(num_out_perm, length(perm[perm==T]))
    }
    
    # Plot
    hist(num_out_perm, xlim=c(0,200), main=paste(inv, zone, (num_out_inv > sort(num_out_perm)[950])))
    lines(c(num_out_inv, num_out_inv), c(-10,1000), col="red")
    
    # Add to table whether significant or not
    tab[inv, zone] = num_out_inv > sort(num_out_perm)[950]
    
  }
}



################################################################################################################
##### PLOT GENOMIC LOCATION OF OUTLIERS ########################################################################
################################################################################################################

# Get number of SNPs and of outliers per 1cM interval
outliers_cM = outliers
outliers_cM$av = round(outliers_cM$av)
outliers_cM$lgav = paste(outliers_cM$LG, outliers_cM$av, sep="_")
outliers_cM$n = TRUE # To maintain total count of SNPs per cM in the next step

outliers_cM = aggregate(outliers_cM[c(liste, "n")],
                        list(outliers_cM$lgav), function(x) length(x[x==T]))
names(outliers_cM)[1] = "lgav"

outliers_cM$LG = as.numeric(data.frame(do.call('rbind', strsplit(as.character(outliers_cM$lgav), "_", fixed=T)), stringsAsFactors=F)$X1)
outliers_cM$av = as.numeric(data.frame(do.call('rbind', strsplit(as.character(outliers_cM$lgav), "_", fixed=T)), stringsAsFactors=F)$X2)


# Get proportion of SNPs that are outliers per 1cM interval
outliers_cM[,liste] = outliers_cM[,liste]/outliers_cM$n

# Keep only positions with a reasonable number of SNPs
outliers_cM10 = outliers_cM[outliers_cM$n>=10, ]


pdf(file = "CZCLI007_Fig2_along_genome.pdf", width = 10, height = 10)
par(mar=c(3.1,3.5,1.1,0.1), mfrow=c(4,5), oma=c(1,3.5,1,1)) # plotting parameters

# All LGs with any outliers
outl_LGs = sort(unique(outliers_cM10$LG[apply(outliers_cM10[,liste], 1, sum)>0]))

for(LG in sort(outl_LGs)){
  # Get results for one LG
  focal = outliers_cM10[outliers_cM10$LG==LG, ]
  focal = focal[order(focal$av), ]
  
  # Get highest map position for this LG
  max_av = read.table(paste0("CZCLI006_comp/CZCLI006_ANG_right_", vtype, ".txt"), header=T, stringsAsFactors=F)
  max_av = max(max_av[max_av$LG==LG, "av"], na.rm=T)
  
  # Plot box for this LG
  plot(focal$av, rep(5,length(focal$av)), col=NA, main=LG, ylim=c(0,7), yaxt='n', pch=16,
       ylab="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, xaxt='n', yaxs='i', xlim=c(0, max_av))
  axis(1, at=seq(0,80,20), seq(0,80,20))
  
  # Make lines separating islands
  lines(c(-10,100),c(6,6), col="grey")
  lines(c(-10,100),c(4,4), col="grey")
  lines(c(-10,100),c(2,2), col="grey")

  # Label islands
  if (LG %in% c(1,6,11,16)){
    mtext("ANG", 2, las=2, at=6.5, line=0.5)
    mtext("CZA", 2, las=2, at=5, line=0.5)
    mtext("CZB", 2, las=2, at=3, line=0.5)
    mtext("CZD", 2, las=2, at=1, line=0.5)
  }
  
  # Plot Rui's inversions as rectangle as background
  focalInvs = invRui[invRui$LG==LG, ]
  if (length(focalInvs$Cluster)>0){
  
    for (number1 in seq(1, length(focalInvs$Cluster))){ # For each inversion in this LG...
      focalInv = focalInvs[number1, ]
      
      for(number2 in 1:7){ # For each zone...
        
        # If significant enrichment of outliers, make rectangle orange; otherwise grey
        col="grey"
        if(length(tab$ANG_right)>0){
          if(tab[focalInv$Cluster, number2]==T){col=rgb(0.9,0.8,0.1)}
        }
        
        rect(focalInv$Start, -number2+7, focalInv$End, -number2+8, col=col, border=NA)
      }
    }
  }
  # Add line between 2 directly adjacent inversions on LG10
  if (LG==10){
    lines(c(invRui[invRui$LG==10,][1,"End"], (invRui[invRui$LG==10,][1,"End"])),
          c(0,100), col="black")
  }
  if (LG==14){
    lines(c(invRui[invRui$LG==14,][1,"End"], (invRui[invRui$LG==14,][1,"End"])),
          c(0,100), col="black")
  }
  
  # Function to determine size of points
  # If less that 2% of SNPs are outliers, set to minimum size
  # If more than 50% of SNPs are outliers, set to maximum size
  cexf = function(x){
    if(is.na(x)){return(0)}
    x = 20*x
    if(x>0 & x<0.4){x=0.4}
    if(x>10){x=10}
    return(x)
  }

  # Plot outlier proportions
  if(sum(focal$ANG, na.rm=T)>0){
    points(focal$av, rep(6.5, length(focal$av)), cex=sapply(focal$ANG_right, cexf))
  }
  if(sum(focal$CZA_left, na.rm=T)>0){
    points(focal$av, rep(5.5, length(focal$av)), cex=sapply(focal$CZA_left, cexf))
  }
  if(sum(focal$CZA_right, na.rm=T)>0){
    points(focal$av, rep(4.5, length(focal$av)), cex=sapply(focal$CZA_right, cexf))
  }
  if(sum(focal$CZB_left, na.rm=T)>0){
    points(focal$av, rep(3.5, length(focal$av)), cex=sapply(focal$CZB_left, cexf))
  }
  if(sum(focal$CZB_right, na.rm=T)>0){
    points(focal$av, rep(2.5, length(focal$av)), cex=sapply(focal$CZB_right, cexf))
  }
  if(sum(focal$CZD_left, na.rm=T)>0){
    points(focal$av, rep(1.5, length(focal$av)), cex=sapply(focal$CZD_left, cexf))
  }
  if(sum(focal$CZD_right, na.rm=T)>0){
    points(focal$av, rep(0.5, length(focal$av)), cex=sapply(focal$CZD_right, cexf))
  }
}
mtext("Map position (cM)", 1, outer=T, at = 0.5)

dev.off()

