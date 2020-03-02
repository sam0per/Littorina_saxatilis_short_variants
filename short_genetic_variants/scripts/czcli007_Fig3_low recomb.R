rm (list=ls())

################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################

liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")

# Get inversion info
invRui = read.table("./data/Sweden_inversions_coordinates_2nd_august_2019.csv", sep=",", header=T,
                    stringsAsFactors=F)
invRui$LG = gsub("LG", "", invRui$LG)

# Get cline fits without inversions
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

outliers$all = rowSums(outliers[, 5:11]) == 7
outliers$any = rowSums(outliers[, 5:11]) >= 1


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



################################################################################################################
### ASSOCIATION WITH LOW-RECOMBINATION REGIONS #################################################################
################################################################################################################

pdf(file = "CZCLI007_Fig3_low_recomb.pdf", width = 12, height = 2)
par(mar=c(2.1,3.5,1.1,0.1), mfrow=c(1,7), oma=c(1,3.5,1,1)) # plotting parameters

for (zone in liste){
  print(zone)
  
  numOut = outliers_cM[,zone]*outliers_cM$n # Number of outliers per 1cM interval
  numNon = outliers_cM$n - numOut # Number of non-outliers per 1cM interval
  
  slopes_perm = c()
  intercepts_perm = c()
  
  # Generate 1000 permutations of outliers across the genome and calculate slope # outliers vs. # non-outliers
  for (number1 in 1:1000){
    
    # Permute genomic positions randomly across SNPs
    perm = sample(rep(outliers_cM$lgav, outliers_cM$n))
    status = c(rep(F, sum(numNon)), rep(T, sum(numOut)))
    
    t1 = table(perm, status)
    
    slopes_perm = append(slopes_perm, lm(t1[,2] ~ t1[,1])$coefficients[[2]])
    intercepts_perm = append(intercepts_perm, lm(t1[,2] ~ t1[,1])$coefficients[[1]])
  }
  
  slope = lm(numOut ~ numNon)$coefficients[[2]]
  intercept = lm(numOut ~ numNon)$coefficients[[1]]

  plot(numNon, numOut, col=rgb(0.5,0.5,0.5,0.4), yaxt='n', xlab="", ylab="",main=zone)
  axis(2, seq(0,100,10), seq(0,100,10), las=1)
  if(zone=="ANG_right"){mtext("Number of outliers", 2, 3)}
  abline(intercept, slope, col=rgb(0.5,0.5,0.5,0.5))
  abline(mean(intercepts_perm), mean(slopes_perm), lty="dashed")
  
  print(slope > sort(slopes_perm,decreasing=F)[950])
}

mtext("Number of non-outliers", side=1, outer=T)
dev.off()
