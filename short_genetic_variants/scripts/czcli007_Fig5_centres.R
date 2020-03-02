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

pdf(file = "CZCLI007_Fig5_centres.pdf", width = 10, height = 15)
#par(mfrow=c(4,2)) # plotting parameters

# Empty plot for ANG_left
#a=hist(snps$Centre, breaks=seq(0,max(snps$Centre,na.rm=T)+2.5, 2.5), col="white",
       #axes=F, border=F, main="", xlab="",ylab="")

par(mfrow=c(4,2), oma=c(3.5,3.5,0,0), xpd=F) # plotting parameters
plot(1, 1, axes=F, col="white", xlab="", ylab="")

for (zone in liste){
  par(xpd=F)
  snps = get(zone)
  
  #Histogram of cline centres for SNPs
  a=hist(snps$Centre, breaks=seq(0,max(snps$Centre,na.rm=T)+2.5, 2.5), pl=F)
  hist(snps$Centre, main=zone,
       breaks=seq(0,max(snps$Centre,na.rm=T)+2.5, 2.5), ylim=c(0, max(a$counts)+10),
       xaxs='i', yaxs='i', xaxt='n', yaxt='n', xlab="", ylab="", col="azure2", border=NA)
  rect(0, -max(a$counts)/4, max(a$breaks), max(a$counts)+10)
  lines(c(-10,300), c(0,0))
  axis(1, seq(0,10000,10), seq(0,10000,10), cex=13)
  axis(2, seq(0,10000,200), seq(0,10000,200), las=2)
  
  # Get info about habitat transitions
  habitat = read.table("./data/Habitat transitions/2_Habitat transitions_new centre.csv", sep=",",
                       header=T, stringsAsFactors=F)
  habitat = habitat[habitat$Island==substr(zone, 1, 3) & habitat$Side==substr(zone, 5, 10), ]
  
  # Add arrows for habitat transitions
  arrows(habitat$Step[habitat$Env=="Rock"], max(a$counts), habitat$Step[habitat$Env=="Rock"], max(a$counts)-max(a$counts)/6,
         lwd=2.5, length=0.1, xpd=T, col="red")
  arrows(habitat$Step[habitat$Env=="Barnacle"], max(a$counts), habitat$Step[habitat$Env=="Barnacle"], max(a$counts)-max(a$counts)/6,
         lwd=2.5, length=0.1, xpd=T, col="black")
  
  # Add phenotypic transitions
  if(zone=="ANG_right"){size=86.09; shape=70.87}
  if(zone=="CZA_left"){size=64.23; shape=59.2}
  if(zone=="CZA_right"){size=79.11; shape=58.24}
  if(zone=="CZB_left"){size=25.1; shape=21.8}
  if(zone=="CZB_right"){size=44.67; shape=21.83}
  if(zone=="CZD_left"){size=71.83; shape=52.94}
  if(zone=="CZD_right"){size=48.75; shape=46.91}
  lines(c(size,size), c(0,10000), lty=1, lwd=1, col="grey")
  lines(c(shape,shape), c(0,10000), lty=2, lwd=1, col="grey")
  
  invs = read.table(paste("./data/Inversion clines/Inversion cline fits 20190717_final invs only/",
                          substr(zone, 1, 3), "_inv_clines_free_20190717.txt", sep=""), header=T, stringsAsFactors=F)
  invs = invs[invs$Side == substr(zone, 5, 10), ]
  invs = invs[invs$Inv %in% c("LGC1.1","LGC1.2","LGC2.1", "LGC4.1", "LGC6.1/2a","LGC7.1", "LGC7.2",
                              "LGC9.1","LGC10.1","LGC10.2","LGC11.1","LGC14.1/2c","LGC14.3","LGC17.1"),]
  
  # Colours to indicate inversions
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
  
  # ltys = rep(NA, length(invs$Inv))
  # ltys[invs$Inv=="LGC1.1"] = 1
  # ltys[invs$Inv=="LGC1.2"] = 2
  # ltys[invs$Inv=="LGC2.1"] = 1
  # ltys[invs$Inv=="LGC4.1"] = 1
  # ltys[invs$Inv=="LGC6.1/2a"] = 3
  # ltys[invs$Inv=="LGC7.1"] = 1
  # ltys[invs$Inv=="LGC7.2"] = 2
  # ltys[invs$Inv=="LGC9.1"] = 1
  # ltys[invs$Inv=="LGC10.1"] = 1
  # ltys[invs$Inv=="LGC10.2"] = 2
  # ltys[invs$Inv=="LGC11.1"] = 1
  # ltys[invs$Inv=="LGC14.1/2c"] = 3
  # ltys[invs$Inv=="LGC14.3"] = 2
  # ltys[invs$Inv=="LGC17.1"] = 1
  
  # Add inversion names to indicate centre of each inversion
  par(xpd=T)
  for (num in 1:length(invs$Inv)){
    if (is.na(cols[num])==F){
      if(invs$Centre[num] > 0 | is.na(invs$Centre[num])){
        text(invs$Centre[num], num*(max(a$counts)/23.4), invs$Inv[num], col=cols[num], adj=0)
      }
    }
  }
}
mtext("Position along shore (m)", 1, -1, outer=T)
mtext("Frequency", 2, -1, outer=T)

dev.off()
