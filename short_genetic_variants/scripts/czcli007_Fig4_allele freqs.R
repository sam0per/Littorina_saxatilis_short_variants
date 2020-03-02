rm (list=ls())

################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################

liste = c("ANG_right", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")

# Get inversion info
invRui = read.table("./data/Sweden_inversions_coordinates_2nd_August_2019.csv", sep=",", header=T,
                    stringsAsFactors=F)
invRui$LG = gsub("LG", "", invRui$LG)

# Get cline fits
ANG_right = read.table("CZCLI006_ANG_right.txt", header=T, stringsAsFactors=F)
CZA_left = read.table("CZCLI006_CZA_left.txt", header=T, stringsAsFactors=F)
CZA_right = read.table("CZCLI006_CZA_right.txt", header=T, stringsAsFactors=F)
CZB_left = read.table("CZCLI006_CZB_left.txt", header=T, stringsAsFactors=F)
CZB_right = read.table("CZCLI006_CZB_right.txt", header=T, stringsAsFactors=F)
CZD_left = read.table("CZCLI006_CZD_left.txt", header=T, stringsAsFactors=F)
CZD_right = read.table("CZCLI006_CZD_right.txt", header=T, stringsAsFactors=F)



################################################################################################################
##### REPEATABILITY OF ALLELE FREQ PATTERNS ####################################################################
################################################################################################################

invClinesANG_right = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "ANG", "_inv_clines_free_20190717.txt", sep=""),
                                header=T, stringsAsFactors=F)
invClinesCZA_left = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZA", "_inv_clines_free_20190717.txt", sep=""),
                               header=T, stringsAsFactors=F)
invClinesCZA_right = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZA", "_inv_clines_free_20190717.txt", sep=""),
                                header=T, stringsAsFactors=F)
invClinesCZB_left = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZB", "_inv_clines_free_20190717.txt", sep=""),
                               header=T, stringsAsFactors=F)
invClinesCZB_right = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZB", "_inv_clines_free_20190717.txt", sep=""),
                                header=T, stringsAsFactors=F)
invClinesCZD_left = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZD", "_inv_clines_free_20190717.txt", sep=""),
                               header=T, stringsAsFactors=F)
invClinesCZD_right = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/", "CZD", "_inv_clines_free_20190717.txt", sep=""),
                                header=T, stringsAsFactors=F)


# Check which of the haplotypes of complex arrangements changes most between C and W end
# This is the one that should be plotted
# (Do this for the LG6 and LG14 inversion)
for (zone in liste){
  
  # Get SNPs outside inversions for this zone
  snps = get(zone)
  snps = snps[snps$invRui==F, ]
  
  # Get inversion cline fits for this zone
  invs = get(paste("invClines", zone, sep=""))
  invs = invs[invs$Side==substr(zone, 5, 9), ]
  invs$p_crab = exp(invs$lp_crab) / (1+exp(invs$lp_crab))
  invs$p_wave = exp(invs$lp_wave) / (1+exp(invs$lp_wave))
  
  invs$p_crab[invs$Allele==0] = 1-invs$p_crab[invs$Allele==0]
  invs$p_wave[invs$Allele==0] = 1-invs$p_wave[invs$Allele==0]
  
  invs$p_diff = invs$p_wave - invs$p_crab
  
  print(invs[invs$Inv %in% c("LGC14.1/2a", "LGC14.1/2b", "LGC14.1/2c"),
             c("Inv", "Allele", "p_crab", "p_wave", "p_diff")])
}


pdf("CZCLI007_Fig4_allele_freqs.pdf", 4, 7.5)
par(mfrow=c(4,2), mar=c(1,1,1,1), oma=c(3.5,3.5,0,0))

#Plot legend first
plot(1, 1, axes=F, col="white")
legend(0.85, 1.45,
       legend=c("LGC1.1","LGC1.2","LGC2.1","LGC4.1","LGC6.1/2a","LGC7.1","LGC7.2","LGC9.1","LGC10.1","LGC10.2",
                "LGC11.1","LGC14.1/2c","LGC14.3","LGC17.1"),
       col= c("blue","blue","light blue","dark green","green","black","black","red",
              "purple","purple","turquoise","darkgoldenrod1","darkgoldenrod1","pink"),
       pch = c(1,2,1,1,8,1,2,1,1,2,1,8,2,1), cex=0.8, bty='n')

# Plot allele freqs for SNPs and inversions, C against W 
for (zone in liste){
  snps = get(zone)
  snps = snps[snps$invRui==F, ]
  
  invs = get(paste("invClines", zone, sep=""))
  invs = invs[invs$Side==substr(zone, 5, 9), ]
  invs$p_crab = exp(invs$lp_crab) / (1+exp(invs$lp_crab))
  invs$p_wave = exp(invs$lp_wave) / (1+exp(invs$lp_wave))

  # Plot ref allele freqs for non-neutliers and outliers separately
  plot(snps$ref_p_wave[snps$invRui==F & snps$sel==F], snps$ref_p_crab[snps$invRui==F & snps$sel==F],
       col=rgb(0.2,0.2,0.2,0.1), xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', pch=16, cex=0.5, main=zone)
  points(snps$ref_p_wave[snps$invRui==F & snps$sel==T], snps$ref_p_crab[snps$invRui==F & snps$sel==T],
       col=rgb(1,0,0,0.3), pch=16)
  axis(1, seq(0,1,0.1), labels = NA)
  axis(2, seq(0,1,0.1), labels = NA)
  axis(1, 0:1, 0:1, las=1)
  axis(2, 0:1, 0:1, las=1)

  # Colours & symbols to indicate inversions
  cols = rep(NA, length(invs$Inv))
  cols[invs$Inv=="LGC1.1"] = "blue"
  cols[invs$Inv=="LGC1.2"] = "blue"
  cols[invs$Inv=="LGC2.1"] = "light blue"
  cols[invs$Inv=="LGC4.1"] = "dark green"
  cols[invs$Inv=="LGC6.1/2a"] = "green"
  cols[invs$Inv=="LGC7.1"] = "black"
  cols[invs$Inv=="LGC7.2"] = "black"
  cols[invs$Inv=="LGC9.1"] = "red"
  cols[invs$Inv=="LGC10.1"] = "purple"
  cols[invs$Inv=="LGC10.2"] = "purple"
  cols[invs$Inv=="LGC11.1"] = "turquoise"
  cols[invs$Inv=="LGC14.1/2c"] = "darkgoldenrod1"
  cols[invs$Inv=="LGC14.3"] = "darkgoldenrod1"
  cols[invs$Inv=="LGC17.1"] = "pink"
  
  symbs = rep(NA, length(invs$Inv))
  symbs[invs$Inv=="LGC1.1"] = 1
  symbs[invs$Inv=="LGC1.2"] = 2
  symbs[invs$Inv=="LGC2.1"] = 1
  symbs[invs$Inv=="LGC4.1"] = 1
  symbs[invs$Inv=="LGC6.1/2a"] = 8
  symbs[invs$Inv=="LGC7.1"] = 1
  symbs[invs$Inv=="LGC7.2"] = 2
  symbs[invs$Inv=="LGC9.1"] = 1
  symbs[invs$Inv=="LGC10.1"] = 1
  symbs[invs$Inv=="LGC10.2"] = 2
  symbs[invs$Inv=="LGC11.1"] = 1
  symbs[invs$Inv=="LGC14.1/2c"] = 8
  symbs[invs$Inv=="LGC14.3"] = 2
  symbs[invs$Inv=="LGC17.1"] = 1

  # Plot inversion frequencies, making sure directionality is always the same
  points(1-invs$p_wave[invs$Allele==0], 1-invs$p_crab[invs$Allele==0], col=cols[invs$Allele==0], pch=symbs[invs$Allele==0], cex=1.5)
  points(invs$p_wave[invs$Allele==2], invs$p_crab[invs$Allele==2], col=cols[invs$Allele==2], pch=symbs[invs$Allele==2], cex=1.5)

}
mtext("Allele frequency (Wave)", 1, outer=T, line=1)
mtext("Allele frequency (Crab)", 2, outer=T, line=1)

dev.off()



# for (zone1 in liste){
#   for (zone2 in liste){
#     
#     focal1 = get(zone1)
#     focal2 = get(zone2)
#     plot(focal1$ref_p_wave[focal1$invRui==F & focal1$sel==F & focal2$sel==F], focal2$ref_p_wave[focal1$invRui==F & focal1$sel==F & focal2$sel==F],
#          col=rgb(0.2,0.2,0.2,0.2), xlim=c(0,1), ylim=c(0,1),
#          xaxt='n', yaxt='n', main=paste(zone1, zone2),
#          xlab="p(Wave) zone 1", ylab="p(Wave) zone 2")
#     outs = focal1$cp[focal1$sel==T & focal2$sel==T & focal1$invRui==F]
#     points(focal1$ref_p_wave[focal1$cp %in% outs], focal2$ref_p_wave[focal2$cp %in% outs],
#            col=rgb(1,0,0,0.3))
#     axis(1, seq(-1,2,0.1), seq(-1,2,0.1), las=1)
#     axis(2, seq(-1,2,0.1), seq(-1,2,0.1), las=1)
#     
#     print(paste(zone1, zone2))
#   }
# }
