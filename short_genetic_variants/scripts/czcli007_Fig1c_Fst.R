# Calculate Fst between same-ecotype samples from different islands, same-ecotype samples from the same
# island, and between ecotypes. Need to use the raw allele frequencies - cline-estimated allele freqs
# exist only for clinal SNPs and therefore Fst between ecotypes will be overestimated.

rm (list=ls())

df = data.frame(matrix(ncol=7, nrow=0))
names(df) = c("island1", "side1", "ecotype1", "island2", "side2", "ecotype2", "Fst")
counter = 1
for (island1 in c("ANG", "CZA", "CZB", "CZD")){
  for (side1 in c("left", "right")){
    for (ecotype1 in c("p_crabRaw", "p_waveRaw")){
      for (island2 in c("ANG", "CZA", "CZB", "CZD")){
        for (side2 in c("left", "right")){
          for (ecotype2 in c("p_crabRaw", "p_waveRaw")){
            
            if((island1=="ANG" & side1=="left") | (island2=="ANG" & side2=="left")){break}

            freqs1 = read.table(paste("CZCLI006_", island1, "_", side1, "noInv", "_snp.txt", sep=""), header=T,
                                stringsAsFactors=F)[,c("cp", ecotype1)]
            freqs2 = read.table(paste("CZCLI006_", island2, "_", side2, "noInv", "_snp.txt", sep=""), header=T,
                                stringsAsFactors=F)[,c("cp", ecotype2)]
            freqs1 = na.omit(freqs1)
            freqs2 = na.omit(freqs2)
            
            use = intersect(freqs1$cp, freqs2$cp)
            freqs1 = freqs1[freqs1$cp %in% use, ]
            freqs2 = freqs2[freqs2$cp %in% use, ]
            freqs1 = freqs1[order(freqs1$cp), ]
            freqs2 = freqs2[order(freqs2$cp), ]
            names(freqs1) = c("cp,", "p")
            names(freqs2) = c("cp,", "p")
            
            
            Hs = (2*freqs1$p*(1-freqs1$p) + 2*freqs2$p*(1-freqs2$p))/2
            Ht = 2 * ((freqs1$p + freqs2$p)/2) * (((1-freqs1$p) + (1-freqs2$p))/2)
            Fst = (Ht - Hs) / Ht

            df[counter, ] = c(island1, side1, ecotype1, island2, side2, ecotype2, mean(Fst, na.rm=T))
            counter = counter+1
            print(counter)
            
          }
        }
      }
    }
  }
}

# Indicate which category a sample pair belongs to ("different islands, Crab ecotype" etc.)
df$cat[df$island1!=df$island2 & df$ecotype1=="p_crabRaw" & df$ecotype2=="p_crabRaw"] = "different_C"
df$cat[df$island1==df$island2 & df$ecotype1=="p_waveRaw" & df$ecotype2=="p_waveRaw" & df$side1!=df$side2] = "same_W"
df$cat[df$island1!=df$island2 & df$ecotype1=="p_waveRaw" & df$ecotype2=="p_waveRaw"] = "different_W"
df$cat[df$island1==df$island2 & df$ecotype1=="p_crabRaw" & df$ecotype2=="p_waveRaw" & df$side1==df$side2] = "same_betwEco"
df$cat[df$island1==df$island2 & df$ecotype1=="p_waveRaw" & df$ecotype2=="p_crabRaw" & df$side1==df$side2] = "same_betwEco"

df$Fst = as.numeric(df$Fst)

pdf(file = "CZCLI007_Fig1c_general_Fst.pdf", width = 9, height = 6)

#par(mar=c(0.5,1,0.5,0.5), mfcol=c(4,6), oma=c(5,5,2,0)) # plotting parameters
par(mar=c(9,4.5,0.5,0.5)) # plotting parameters

pl=boxplot(df$Fst[df$cat=="different_C" & is.na(df$cat)==F], df$Fst[df$cat=="same_W" & is.na(df$cat)==F],
        df$Fst[df$cat=="different_W" & is.na(df$cat)==F], df$Fst[df$cat=="same_betwEco" & is.na(df$cat)==F],
        col=c("dark red", "blue", "blue", "grey"),
        xaxt="n", yaxt="n", ylab=expression('Average F'[ST]), ylim=c(0,0.065),
        cex.lab=1)
axis(1, 1:4, c("Within Crab,\nbetween islands",
               "Within Wave,\nwithin island", "Within Wave,\nbetween islands",
               "Between ecotypes,\nwithin islands"), tick=F,
     las=2, cex=0.7)

axis(2, seq(0,0.11,0.01), seq(0,0.11,0.01), las=2, cex.lab=0.1)
dev.off()

