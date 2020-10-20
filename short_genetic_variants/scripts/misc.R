A_dir <- "/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115/"
txt_fl <- list.files(path = "CZCLI006_comp", pattern = "_CZ", full.names = TRUE)[c(1,7)]
tbl <- lapply(txt_fl, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
identical(tbl[[1]]$cp, tbl[[2]]$cp)

A_txt_fl <- list.files(path = paste0(A_dir, "CZCLI006_comp"), pattern = "_CZ", full.names = TRUE)[c(2,14)]
A_tbl <- lapply(A_txt_fl, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
identical(A_tbl[[1]]$cp, A_tbl[[2]]$cp)

head(tbl[[1]])
head(A_tbl[[1]])
cp_diff <- data.frame(cp = setdiff(A_tbl[[1]]$cp, tbl[[1]]$cp))
cp_diff <- data.frame(cp = setdiff(A_tbl[[1]][A_tbl[[1]]$sel==FALSE, "cp"], tbl[[1]][tbl[[1]]$sel==FALSE, "cp"]))
sum(A_tbl[[1]]$sel)
sum(tbl[[1]]$sel)
cp_diff <- data.frame(cp = setdiff(tbl[[1]]$cp, A_tbl[[1]]$cp))
head(cp_diff)
library(tidyr)
contig_diff <- separate(data = cp_diff, col = cp, into = c("contig", "pos"), sep = "_")
head(contig_diff)

contig_int <- read.table(file = "Littorina_saxatilis/short_genetic_variants/contig_intervals.list")
head(contig_int)
id <- round(runif(n = 1, min = 1, max = nrow(contig_diff)))
as.character(contig_int[, 1])[grepl(pattern = paste0(contig_diff$contig[id], ":"), x = as.character(contig_int[, 1]))]
contig_diff$pos[id]
contig_diff$contig[id]

# 
# 
# 

# MARKER DENSITY
tt <- mer[mer$ISL=='CZA' & mer$ECOT=='WAVE_LEFT', ]
head(tt)
tt <- tt[tt$x.x!=0 & tt$x.y!=0, ]

# perspective plot
library(MVN)
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "persp")
# contour plot
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "contour")

library(bestNormalize)
(BNobject.x <- bestNormalize(tt$prop.x))
(BNobject.y <- bestNormalize(tt$prop.y))

identical(BNobject.x$x, tt$prop.x)
identical(BNobject.x$x.t, BNobject.x$chosen_transform$x.t)

hist(tt$prop.x)
hist(tt$prop.y)
MASS::truehist(BNobject.x$chosen_transform$x.t)
MASS::truehist(BNobject.y$x.t)
result <- mvn(data.frame(BNobject.x$x.t, BNobject.y$x.t), mvnTest = "hz", multivariatePlot = "persp")
# contour plot
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "contour")

sunflowerplot(tt$prop.x, tt$prop.y)

# install the package 
# install.packages("ggstatsplot")

# Load the package
library(ggstatsplot)

# Load the dataset 
data("warpbreaks")

# Create a boxplot of the dataset, outliers are shown as two distinct points
boxplot(tt$prop.x)$out

#Create a boxplot that labels the outliers  
head(tt)
# ggbetweenstats(tt, x.x, x.y, outlier.tagging = TRUE)
Q <- quantile(tt$prop.x, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(tt$prop.x)
up <- Q[2]+1.5*iqr # Upper Range  
low <- Q[1]-1.5*iqr # Lower Range
eliminated <- subset(tt, tt$prop.x > (Q[1] - 1.5*iqr) & tt$prop.x < (Q[2]+1.5*iqr))

outliers <- boxplot(tt$prop.x, plot=FALSE)$out
tt_notl <- tt[-which(tt$prop.x %in% outliers),]

# boxplot(eliminated$prop.x)$out
boxplot(tt_notl$prop.x)$out

outliers <- boxplot(tt_notl$prop.y, plot=FALSE)$out
tt_notl <- tt_notl[-which(tt_notl$prop.y %in% outliers),]

tt <- tt_notl
hist(tt$prop.x)

# 
# 
# 

library(tidyverse)
library(caret)
# Load the data
data("swiss")
# Inspect the data
sample_n(swiss, 3)
# Define training control
set.seed(123) 
train.control <- trainControl(method = "cv", number = 10)
# Train the model
model <- train(Fertility ~., data = swiss, method = "lm",
               trControl = train.control)
# Summarize the results
print(model)
train.control <- trainControl(method = "repeatedcv", 
                              number = 10, repeats = 3)
# Train the model
model <- train(Fertility ~., data = swiss, method = "lm",
               trControl = train.control)
# Summarize the results
print(model)

model$finalModel
summary(lm(formula = Fertility ~., data = swiss))

confint(model$finalModel)
confint(object = lm(formula = Fertility ~., data = swiss))

# 
# 
# 
# Generate the site frequency spectrum
# For a neutral model with two populations and migration:
install.packages('coala')
library(coala)
model <- coal_model(20, 2000) +
  feat_mutation(2) +
  feat_recombination(1) +
  sumstat_tajimas_d()
stats <- simulate(model, seed = 15)
plot(density(stats$tajimas_d, na.rm = TRUE), 
     main = "Neutral Distribution of Tajiam's D")

model2a <- coal_model(c(10, 10), 100) +
  feat_mutation(10) +
  feat_recombination(5) +
  feat_migration(0.5, symmetric = TRUE) +
  sumstat_sfs(population = "all")
stats <- simulate(model2a, seed = 20)
barplot(stats$sfs / sum(stats$sfs), 
        names.arg = seq_along(stats$sfs), 
        col = 3)



install.packages('shallot')
library(shallot)
pd1 <- ewens(mass(1),50)
decay <- decay.exponential(temperature(1.0),dist(scale(USArrests)))
attraction <- attraction(permutation(n.items=50,fixed=FALSE), decay)
pd2 <- ewens.pitman.attraction(mass(1), discount(0.05), attraction)
pd3 <- ddcrp(mass(1), attraction)
# 
# 
# 
eff <- c('coding_sequence_variant', 'chromosome', 'duplication', 'inversion','inframe_insertion',
         'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion',
         'downstream_gene_variant', 'exon_variant', 'exon_loss_variant', 'frameshift_variant',
         'gene_variant', 'feature_ablation', 'gene_fusion', 'bidirectional_gene_fusion', 'rearranged_at_DNA_level',
         'intergenic_region', 'conserved_intergenic_variant', 'intragenic_variant', 'intron_variant',
         'conserved_intron_variant', 'miRNA', 'missense_variant', 'initiator_codon_variant', 'stop_retained_variant',
         'protein_protein_contact', 'structural_interaction_variant', 'rare_amino_acid_variant', 'splice_acceptor_variant',
         'splice_donor_variant', 'splice_region_variant', 'stop_lost', '5_prime_UTR_premature start_codon_gain_variant',
         'start_lost', 'stop_gained', 'synonymous_variant', 'start_retained', 'stop_retained_variant', 'transcript_variant',
         'regulatory_region_variant', 'upstream_gene_variant', '3_prime_UTR_variant', '3_prime_UTR_truncation + exon_loss',
         '5_prime_UTR_variant', '5_prime_UTR_truncation + exon_loss_variant', 'sequence_feature + exon_loss_variant')
# unique(length(eff))
# data.frame(eff, t = 'FALSE')
coding <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
            TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
            FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
            FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
            FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE)
impact <- c('LOW', 'HIGH', 'HIGH', 'HIGH', 'MODERATE', 'MODERATE', 'MODERATE', 'MODERATE', 'MODIFIER', 'MODIFIER',
            'HIGH', 'HIGH', 'MODIFIER', 'HIGH', 'HIGH', 'HIGH', 'HIGH', 'MODIFIER', 'MODIFIER', 'MODIFIER',
            'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODERATE', 'LOW', 'LOW', 'HIGH', 'HIGH', 'HIGH', 'HIGH',
            'HIGH', 'LOW', 'HIGH', 'LOW', 'HIGH', 'HIGH', 'LOW', 'LOW', 'LOW', 'MODIFIER',
            'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODERATE', 'MODIFIER', 'MODERATE', 'MODIFIER')
nnsyn <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
           TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
           FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
           FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
           FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE)
dt_eff <- data.frame(eff, impact, coding, noncoding=!coding, nnsyn)
dt_eff$syn <- ifelse(grepl(pattern = 'synonymous|_retained', x = dt_eff$eff), yes = TRUE, no = FALSE)
# sum(dt_eff$syn)
dt_eff$ins <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
                TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
                TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
dt_eff$del <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
                TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
                TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
                TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

write.table(x = dt_eff, file = 'data/ANN_snpeff_classes.csv', quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
# 
# 
# 
as.character(dtn$ANN)[sapply(dtn$ANN, FUN = function(x) {
  stsp <- strsplit(x = as.character(x), split = '\\|')[[1]]
  hi <- sum(stsp %in% snpeff2)>0
})]
dtn[sapply(dtn$ANN, FUN = function(x) {
  stsp <- strsplit(x = as.character(x), split = '\\|')[[1]]
  hi <- sum(stsp %in% snpeff2)>0
}),]
# 
# 
# 
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 12, name = "Paired")
brewer.pal(n = 12, name = "Paired")[7:8]

display.brewer.pal(n = 9, name = "Greens")
brewer.pal(n = 9, name = "Greens")

display.brewer.pal(n = 9, name = "YlGnBu")
brewer.pal(n = 9, name = "YlGnBu")

display.brewer.pal(n = 9, name = "Greys")
brewer.pal(n = 9, name = "Greys")

display.brewer.pal(n = 9, name = "Purples")
brewer.pal(n = 9, name = "Purples")

display.brewer.pal(n = 11, name = "Dark2")
brewer.pal(n = 8, name = "Dark2")

display.brewer.pal(n = 9, name = "Reds")
brewer.pal(n = 9, name = "Reds")[c(3,4,5)]
display.brewer.pal(n = 9, name = "Blues")
brewer.pal(n = 9, name = "Blues")[c(2,4,5)]
display.brewer.pal(n = 9, name = "YlOrBr")
brewer.pal(n = 9, name = "YlOrBr")[c(3,4,5)]
# 
# 
# 
## PREPARE TABLE FOR TAJIMA'S D EXPECTATIONS
vart <- c('INDEL', 'INS', 'DEL', 'SNP')
cate <- c('noncoding', 'coding', 'frameshift', 'inframe', 'synonymous', 'nonsynonymous')
vc <- unlist(lapply(X = vart, FUN = function(x) paste(x, cate, sep = '_')))

colClasses = rep(x = "character", length(vc))
col.names = vc
df <- read.table(text = "",
                 colClasses = colClasses,
                 col.names = col.names, nrows = length(vc))
row.names(df) <- vc

df <- data.frame(matrix(data = 0, nrow = length(vc), ncol = length(vc),
                        dimnames = list(vc, vc)),
                 stringsAsFactors=F)

# write.table(x = df, file = 'data/ANN_expected_tajD.csv', quote = FALSE, sep = ',', row.names = TRUE, col.names = TRUE)
# 
# 
# 
## TAJIMA'S D
dtp_or <- dtp
dtp <- dtp_or[dtp_or$ANN==tv[1],]

pi_n <- 2 * dtp$DAC * ((2*dtp$N) - dtp$DAC)
# 2*10
# (2*56)-10
# 20*102
pi_d <- (2*dtp$N) * ((2*dtp$N) - 1)
# 2*56
# 112*111
# pi_d[1]
e_pi <- sum(pi_n / pi_d)
# sum(c(2,4,1)/c(2,2,10))

head(dtp)
e_theta <- nrow(dtp) / sum(1/seq(1, 2*unique(dtp$N)-1))

variance.d <- function(n,S) {
  a1=sum(1/(seq(from=1, to=(n-1), by=1)))
  a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
  b1=(n+1)/(3*(n-1))
  b2=(2*((n**2)+n+3))/((9*n)*(n-1))
  c1=b1 - (1/a1)
  c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
  e1=c1/a1
  e2=c2/((a1**2)+a2)
  var=(e1*S) + (e2*S*(S-1))
  return(var)
}
varD <- variance.d(n = 2*unique(dtp$N), S = nrow(dtp))
(e_pi - e_theta) / sqrt(varD)

# (7.120912 - 8.541293) / sqrt(9.629897)
# 
# 
# 
## GC-biased gene conversion
dt$gBGC <- NULL
dt$ANC <- ifelse(test = dt$NE_W_Lcomp=='ref_anc:ref_anc', yes = as.character(dt$REF), no = as.character(dt$ALT))
dt$DER <- ifelse(test = dt$NE_W_Lcomp=='ref_anc:ref_anc', yes = as.character(dt$ALT), no = as.character(dt$REF))
dt$A_gBGC <- ifelse(test = dt$ANC=='A' | dt$ANC=='T', yes = 'W', no = 'S')
head(dt)
tail(dt)
dt$D_gBGC <- ifelse(test = dt$DER=='A' | dt$DER=='T', yes = 'W', no = 'S')
dt$gBGC <- paste(dt$A_gBGC, dt$D_gBGC, sep = '')
table(dt$gBGC)
write.table(x = dt, file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv', quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
# 
# 
# 
