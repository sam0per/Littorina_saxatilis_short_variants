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
# dtp_or <- dtp
# dtp <- dtp_or[dtp_or$ANN==tv[1],]
dtp <- read.table(file = 'results/marker_density/MD_CZA_CRAB_coding_WWSS_count.txt', header = TRUE)
dtp <- read.table(file = 'results/marker_density/MD_CZA_CRAB_noncoding_WWSS_count.txt', header = TRUE)
head(dtp)

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


dt <- unique(read.csv(file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv'))
table(dt$gBGC)
head(dt[dt$gBGC=='WW', ])

sum(paste0(dt$REF, dt$ALT)=='AA')
sum(paste0(dt$REF, dt$ALT)=='AT')
sum(paste0(dt$REF, dt$ALT)=='TA')
sum(paste0(dt$REF, dt$ALT)=='TT')

sum(paste0(dt$REF, dt$ALT)=='GG')
sum(paste0(dt$REF, dt$ALT)=='GC')
sum(paste0(dt$REF, dt$ALT)=='CG')
sum(paste0(dt$REF, dt$ALT)=='CC')

dt$CLASS <- dt$gBGC
table(dt[dt$gBGC=='WW' | dt$gBGC=='SS', 'gBGC'])
table(dt[dt$gBGC=='WW' | dt$gBGC=='SS', 'CLASS'])
write.table(x = dt[dt$gBGC=='WW' | dt$gBGC=='SS', ], file = 'results/Lsax_short_WWSS_czs_daf_inv_findv.csv', quote = FALSE,
            sep = ',', row.names = FALSE, col.names = TRUE)

dt <- unique(read.csv(file = 'results/Lsax_short_WWSS_czs_daf_inv_findv.csv'))
head(dt)
dt$CLASS <- 'WWSS'
table(dt$CLASS)
table(dt$gBGC)
write.table(x = dt, file = 'results/Lsax_short_WWSS_czs_daf_inv_findv.csv', quote = FALSE,
            sep = ',', row.names = FALSE, col.names = TRUE)

table(dt$gBGC)
table(dt[dt$gBGC=='SW' | dt$gBGC=='WS', 'gBGC'])
dt <- split(x = dt, f = dt$gBGC)
write.table(x = dt$SW, file = 'results/Lsax_short_SW_czs_daf_inv_findv.csv', quote = FALSE,
            sep = ',', row.names = FALSE, col.names = TRUE)
write.table(x = dt$WS, file = 'results/Lsax_short_WS_czs_daf_inv_findv.csv', quote = FALSE,
            sep = ',', row.names = FALSE, col.names = TRUE)

dt <- unique(read.csv(file = 'results/Lsax_short_SW_czs_daf_inv_findv.csv'))
head(dt)
dt$CLASS <- dt$gBGC
table(dt$CLASS)
table(dt$gBGC)
write.table(x = dt, file = 'results/Lsax_short_SW_czs_daf_inv_findv.csv', append = FALSE, quote = FALSE,
            sep = ',', row.names = FALSE, col.names = TRUE)
# 
# 
# 
# 
ACns_WWSS <- read.table('results/marker_density/MD_CZA_CRAB_nonsyn_WWSS_count.txt', header = TRUE)
head(ACns_WWSS)
length(unique(as.character(ACns_WWSS$CHROM)))
p <- ggplot(data = ACns_WWSS, aes(x = DAC)) +
  geom_histogram(binwidth = 1, fill = 'grey', col ='black') +
  labs(title = 'CZA CRAB nonsyn WWSS') +
  theme(title = element_text(size = 20))
p
ggsave(filename = 'test/figures/SFS_CZA_CRAB_nonsyn_WWSS.pdf', plot = p, dpi = 'screen')

Ad <- unique(as.character(ACns_WWSS[duplicated(ACns_WWSS$CHROM), 'CHROM']))
unique(Ad)
nrow(ACns_WWSS[ACns_WWSS$CHROM %in% unique(Ad), ])
pos_dac <- ACns_WWSS[ACns_WWSS$CHROM %in% Ad, c('CHROM', 'POS', 'cp', 'LG', 'av','DAC', 'ZONE')]
data.frame(table(pos_dac$CHROM))
cor(x = pos_dac$POS, y = pos_dac$DAC)
ggplot(data = pos_dac, aes(x = POS, y = DAC)) +
  facet_wrap(~CHROM) +
  geom_point(alpha = 0.3)
ggplot(data = pos_dac, aes(x = av, y = DAC)) +
  facet_wrap(~LG) +
  geom_point(alpha = 0.3)
# plot(x = pos_dac$POS, y = pos_dac$DAC)
# diff(x = pos_dac$POS)
data.frame(table(ACns_WWSS[ACns_WWSS$CHROM %in% Ad, 'DAC']))

# hist(ACns_WWSS$DAC, breaks = unique(ACns_WWSS$N)*2+1, include.lowest = TRUE)
ggplot(data = ACns_WWSS, aes(x = DAC)) +
  geom_histogram(binwidth = 1, fill = 'grey', col ='black') +
  geom_histogram(data = ACns_WWSS[ACns_WWSS$CHROM %in% Ad, ], aes(x = DAC), fill = 'blue', col = 'black', binwidth = 1)
data.frame(table(ACns_WWSS$DAC))
ACns_WWSS[ACns_WWSS$DAC==82, 'cp']
ACns_WWSS[ACns_WWSS$CHROM=='Contig40948', ]
ACns_WWSS[ACns_WWSS$CHROM=='Contig9412', ]

BCns_WWSS <- read.table('results/marker_density/MD_CZB_CRAB_nonsyn_WWSS_count.txt', header = TRUE)
head(BCns_WWSS)
length(unique(as.character(BCns_WWSS$CHROM)))
# hist(BCns_WWSS$DAC, breaks = unique(BCns_WWSS$N)*2)
ggplot(data = BCns_WWSS, aes(x = DAC)) +
  geom_histogram(binwidth = 1)

Bd <- unique(as.character(BCns_WWSS[duplicated(BCns_WWSS$CHROM), 'CHROM']))
unique(Bd)
nrow(BCns_WWSS[BCns_WWSS$CHROM %in% unique(Bd), ])
Bpos_dac <- BCns_WWSS[BCns_WWSS$CHROM %in% Bd, c('CHROM', 'POS', 'cp', 'LG', 'av','DAC', 'ZONE')]




DCns_WWSS <- read.table('results/marker_density/MD_CZD_CRAB_nonsyn_WWSS_count.txt', header = TRUE)
head(DCns_WWSS)
length(unique(as.character(DCns_WWSS$CHROM)))
# hist(DCns_WWSS$DAC, breaks = unique(DCns_WWSS$N)*2)
ggplot(data = DCns_WWSS, aes(x = DAC)) +
  geom_histogram(binwidth = 1)
data.frame(table(DCns_WWSS$DAC))
DCns_WWSS[DCns_WWSS$DAC==1, ]
DCns_WWSS[DCns_WWSS$DAC==1, 'cp']

DCns_WWSS[DCns_WWSS$cp=='Contig40948_77418', 'DAC']
DCns_WWSS[DCns_WWSS$cp=='Contig9412_15107', 'DAC']

BCns_WWSS[BCns_WWSS$cp=='Contig40948_77418', 'DAC']
BCns_WWSS[BCns_WWSS$cp=='Contig9412_15107', 'DAC']

Dd <- unique(as.character(DCns_WWSS[duplicated(DCns_WWSS$CHROM), 'CHROM']))
unique(Dd)
nrow(DCns_WWSS[DCns_WWSS$CHROM %in% unique(Dd), ])
Dpos_dac <- DCns_WWSS[DCns_WWSS$CHROM %in% Dd, c('CHROM', 'POS', 'cp', 'LG', 'av','DAC', 'ZONE')]

pd <- rbind(pos_dac, Bpos_dac, Dpos_dac)

ggplot(data = pd, aes(x = POS, y = DAC, col = ZONE)) +
  facet_wrap(~CHROM) +
  geom_point(alpha = 0.4, size = 3)
ggplot(data = pd, aes(x = av, y = DAC, col = ZONE)) +
  facet_wrap(~LG) +
  geom_point(alpha = 0.3)

rm(list = ls())
(fl <- list.files(path = 'results/marker_density', pattern = 'nonsyn_WWSS', full.names = TRUE))

ll <- list()
for (i in 1:length(fl)) {
  dt <- read.table(file = fl[i], header = TRUE)
  dd <- unique(as.character(dt[duplicated(dt$CHROM), 'CHROM']))
  ll[[i]] <- dd
}
sum(sapply(ll, length))
length(unlist(ll))
cc <- data.frame(table(unlist(ll)))
sum(cc$Freq)
dc <- as.character(cc[cc$Freq==length(fl), 1])

library(data.table)
dt <- as.data.frame(rbindlist(lapply(fl, read.table, header = TRUE)))
# ds <- split(x = dt, f = dt$ZONE)
pos_dac <- dt[dt$CHROM %in% dc, c('CHROM', 'POS', 'cp', 'LG', 'av','DAC', 'ZONE', 'ECOT')]

dd <- unique(as.character(dt[duplicated(dt$CHROM), 'CHROM']))
length(unique(dd))
# nrow(DCns_WWSS[DCns_WWSS$CHROM %in% unique(Dd), ])
pos_dac <- dt[dt$CHROM %in% dd, c('CHROM', 'POS', 'cp', 'LG', 'av','DAC', 'ZONE', 'ECOT')]
pos_dac$ZE <- paste(pos_dac$ZONE, pos_dac$ECOT)

p <- ggplot(data = pos_dac, aes(x = POS, y = DAC, col = ZONE)) +
  facet_wrap(~CHROM) +
  geom_point(alpha = 0.4, size = 3) +
  theme(title = element_text(size = 20), text = element_text(size = 14))
ggsave(filename = 'test/figures/pos_dac_czs_dup.pdf', plot = p, dpi = 'screen')
q <- ggplot(data = pos_dac, aes(x = av, y = DAC, col = ZONE)) +
  facet_wrap(~LG) +
  geom_point(alpha = 0.4, size = 3) +
  labs(x = 'Map position') +
  theme(title = element_text(size = 20), text = element_text(size = 14))
ggsave(filename = 'test/figures/av_dac_czs_dup.pdf', plot = q, dpi = 'screen')

aa <- read.table('results/marker_density/MD_CZD_CRAB_nonsyn_WWSS_count.txt', header = TRUE)
aa[aa$CHROM=='Contig8325',]
aa[aa$DAC==82,]

# ggplot(data = pos_dac, aes(x = av, y = ZE, col = DAC)) +
#   facet_wrap(~LG) +
#   geom_point(alpha = 0.3, size = 2)

ff <- list.files(path = 'CZCLI006_comp', pattern = 'NoInv_SNP', full.names = TRUE)
et <- as.data.frame(rbindlist(lapply(ff, read.table, header = TRUE)))
table(et$sel)
cl <- et[et$Contig %in% dc, ]
table(cl$sel)
table(cl$Type)

plot(x = pos_dac$POS, y = pos_dac$DAC)

pd <- split(x = pos_dac, f = pos_dac$ZONE)
pd <- merge(x = pd$CZA, y = pd$CZB, by = 'CHROM')
plot(x = pd$DAC.x, y = pd$DAC.y)
cor.test(x = pd$DAC.x, y = pd$DAC.y)

pd <- split(x = pos_dac, f = pos_dac$ZONE)
pd <- merge(x = pd$CZA, y = pd$CZD, by = 'CHROM')
plot(x = pd$DAC.x, y = pd$DAC.y)
cor.test(x = pd$DAC.x, y = pd$DAC.y)

pd <- split(x = pos_dac, f = pos_dac$ZONE)
pd <- merge(x = pd$CZB, y = pd$CZD, by = 'CHROM')
plot(x = pd$DAC.x, y = pd$DAC.y)
cor.test(x = pd$DAC.x, y = pd$DAC.y)
summary(lm(formula = DAC.x ~ DAC.y * ZE.y, data = pd))
ggplot(data = pd, aes(x = DAC.x, y = DAC.y)) +
  geom_abline(slope = 1) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x)


aa <- read.table('results/marker_density/MD_CZB_CRAB_nonsyn_WWSS_count.txt', header = TRUE)
head(aa)
p <- ggplot(data = aa, aes(x = DAC)) +
  geom_histogram(binwidth = 1, fill = 'grey', col ='black') +
  labs(title = 'CZB CRAB nonsyn WWSS') +
  theme(title = element_text(size = 20))
p
ggsave(filename = 'test/figures/SFS_CZB_CRAB_nonsyn_WWSS.pdf', plot = p, dpi = 'screen')
# 
# 
# 
# INDEL VS SNP COUNT
rm(list = ls())
(ff <- list.files(path = 'results/marker_density', pattern = 'noncoding', full.names = TRUE))
(ff <- ff[grepl(pattern = "DEL|INS", x = ff)])
# INDEL
tt <- vector(mode = "integer", length = length(ff))
for (i in seq(from = 1, to = 18, by = 2)) {
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  oo <- c(as.character(dd$cp), as.character(ee$cp))
  # print(length(unique(oo)))
  tt[i] <- length(unique(oo))
}
tt
nonc <- tt[tt!=0]
pop <- c('CZA_C', 'CZA_WL', 'CZA_WR', 'CZB_C', 'CZB_WL', 'CZB_WR', 'CZD_C', 'CZD_WL', 'CZD_WR')
# 
# 
# 
# NON-CODING VS CODING
(ff <- list.files(path = 'results/marker_density', pattern = '_coding', full.names = TRUE))
(ff <- ff[grepl(pattern = "DEL|INS", x = ff)])

tt <- vector(mode = "integer", length = length(ff))
tt
for (i in seq(from = 1, to = 18, by = 2)) {
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  oo <- c(as.character(dd$cp), as.character(ee$cp))
  # print(length(unique(oo)))
  tt[i] <- length(unique(oo))
}
cod <- tt[tt!=0]
da <- data.frame(pop, nonc, cod)
da$rat <- da$nonc/da$cod
da[which.min(da$rat),]
da[which.max(da$rat),]
da$vt <- 'INDEL'
ind <- da[, c(-1,-4)]
# SNP
rm(list = setdiff(ls(), 'ind'))
(ff <- list.files(path = 'results/marker_density', pattern = 'noncoding', full.names = TRUE))
(ff <- ff[!grepl(pattern = "DEL|INS", x = ff)])
tt <- vector(mode = "integer", length = length(ff))
for (i in seq(from = 1, to = 27, by = 3)) {
  # i <- 4
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  gg <- read.table(file = ff[i+2], header = TRUE)
  oo <- c(as.character(dd$cp), as.character(ee$cp), as.character(gg$cp))
  # print(length(unique(oo)))
  tt[i] <- length(unique(oo))
}
nonc <- tt[tt!=0]
pop <- c('CZA_C', 'CZA_WL', 'CZA_WR', 'CZB_C', 'CZB_WL', 'CZB_WR', 'CZD_C', 'CZD_WL', 'CZD_WR')
# 
# 
# 
# NON-CODING VS CODING
(ff <- list.files(path = 'results/marker_density', pattern = '_coding', full.names = TRUE))
(ff <- ff[!grepl(pattern = "DEL|INS", x = ff)])

tt <- vector(mode = "integer", length = length(ff))
tt
for (i in seq(from = 1, to = 27, by = 3)) {
  # i <- 4
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  gg <- read.table(file = ff[i+2], header = TRUE)
  oo <- c(as.character(dd$cp), as.character(ee$cp), as.character(gg$cp))
  # print(length(unique(oo)))
  tt[i] <- length(unique(oo))
}
cod <- tt[tt!=0]
da <- data.frame(pop, nonc, cod)
da$rat <- da$nonc/da$cod
mean(da$rat)
da[which.min(da$rat),]
da[which.max(da$rat),]
da$vt <- 'SNP'
vv <- rbind(ind, da[, c(-1,-4)])
cq <- chisq.test(vv[,1:2])
cq$observed
cq$expected
# 
# 
# 
# INDEL SIZE
(ff <- list.files(path = 'results/marker_density', pattern = '_noncoding_DEL', full.names = TRUE))
for (i in 1:9) {
  dd <- read.table(file = ff[i], header = TRUE)
  dd <- dd[dd$SIZE %in% c(1, 2, 3, 6), ]
  # dd <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_DEL_count.txt', header = TRUE)
  # ggplot(data = rbind(dd,ii), aes(x = SIZE, fill = VTYPE)) +
  #   geom_histogram(binwidth = 1, col = 'black', position = 'dodge')
  slices <- data.frame(table(dd$SIZE))[, 2]
  # lbls <- c("US", "UK", "Australia", "Germany", "France")
  # lbls <- c('SW','WS', 'WWSS')
  lbls <- data.frame(table(dd$SIZE))[, 1]
  size <- as.character(lbls) %in% c("1", "2", "3", "6")
  pct <- round(slices/sum(slices)*100, digits = 1)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main='DEL')
  # print(sum(pct[1:9]))
  # sum(pct[1:6])
  # print(sum(pct[1:3]))
  # print(sum(pct[1:2]))
  # print(sum(pct[size]))
  print(chisq.test(slices))
}
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("cardiomoon/moonBook")
# devtools::install_github("cardiomoon/webr")
# require(ggplot2)
# require(moonBook)
# install.packages('webr')
# library(webr)
# PieDonut(acs,aes(pies=Dx,donuts=smoking))
# 
# 
# 
rm(list = ls())
(ff <- list.files(path = 'results/marker_density', pattern = '_coding_INS', full.names = TRUE))
for (i in 1:9) {
  dd <- read.table(file = ff[i], header = TRUE)
  dd <- dd[dd$SIZE %in% c(1, 2, 3, 6), ]
  # dd <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_DEL_count.txt', header = TRUE)
  # ggplot(data = rbind(dd,ii), aes(x = SIZE, fill = VTYPE)) +
  #   geom_histogram(binwidth = 1, col = 'black', position = 'dodge')
  slices <- data.frame(table(dd$SIZE))[, 2]
  # lbls <- c("US", "UK", "Australia", "Germany", "France")
  # lbls <- c('SW','WS', 'WWSS')
  lbls <- data.frame(table(dd$SIZE))[, 1]
  size <- as.character(lbls) %in% c("1", "2", "3", "6")
  pct <- round(slices/sum(slices)*100, digits = 1)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main='INS')
  # print(sum(pct[1:9]))
  # sum(pct[1:6])
  # print(sum(pct[1:3]))
  # print(sum(pct[size]))
  # print(sum(pct[1:2]))
  # print(sum(pct[3:5]))
  if (length(slices)==3) {
    slices <- c(slices,0)
  }
  print(chisq.test(slices))
}

# 
# 
# 
ee <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_syn_DEL_count.txt', header = TRUE)
er <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_DEL_count.txt', header = TRUE)
et <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_syn_INS_count.txt', header = TRUE)
ey <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_INS_count.txt', header = TRUE)
ee <- rbind(ee, er, et, ey)
we <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_syn_DEL_count.txt', header = TRUE)
wr <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_nonsyn_DEL_count.txt', header = TRUE)
wt <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_syn_INS_count.txt', header = TRUE)
wy <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_nonsyn_INS_count.txt', header = TRUE)
we <- rbind(we, wr, wt, wy)
length(intersect(ee$cp, we$cp))
length(ee$cp)
length(we$cp)

ee <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_syn_SW_count.txt', header = TRUE)
er <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_SW_count.txt', header = TRUE)
et <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_syn_WS_count.txt', header = TRUE)
ey <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_WS_count.txt', header = TRUE)
eu <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_syn_WWSS_count.txt', header = TRUE)
ei <- read.table(file = 'results/marker_density/MD_CZA_WAVE_LEFT_nonsyn_WWSS_count.txt', header = TRUE)
ee <- rbind(ee, er, et, ey, eu, ei)
we <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_syn_SW_count.txt', header = TRUE)
wr <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_nonsyn_SW_count.txt', header = TRUE)
wt <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_syn_WS_count.txt', header = TRUE)
wy <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_nonsyn_WS_count.txt', header = TRUE)
wu <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_syn_WWSS_count.txt', header = TRUE)
wi <- read.table(file = 'results/marker_density/MD_CZB_WAVE_LEFT_nonsyn_WWSS_count.txt', header = TRUE)
we <- rbind(we, wr, wt, wy, wu, wi)
length(intersect(ee$cp, we$cp))/length(ee$cp)
length(intersect(ee$cp, we$cp))/length(we$cp)
# 
# 
# 
# COUNT OF DIFFERENT TYPES OF VARIANTS
get_ro <- function(z, e, v, a) {
  s <- dh_sub[dh_sub$ZONE %in% z & dh_sub$ECOT==e & dh_sub$Variant_type %in% v & dh_sub$ANN %in% a, 'segsites']
  return(s)
}
levels(dh_sub$ANN)
vv <- 'INS'
aa <- c('nongenic', 'nonsyn', 'syn')
# Pie Chart with Percentages
# slices <- c(10, 12, 4, 16, 8)
slices <- get_ro(z = 'CZD', e = 'WAVE_LEFT', v = vv, a = aa)
# lbls <- c("US", "UK", "Australia", "Germany", "France")
# lbls <- c('SW','WS', 'WWSS')
lbls <- aa
pct <- round(slices/sum(slices)*100, digits = 1)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main=vv)

an <- c('nongenic', 'nonsyn', 'syn')
dh_sub <- dh_res[dh_res$ANN %in% an, ]

ggplot(data = dh_sub, aes(x = Variant_type, y = segsites, fill = ANN)) +
  facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_col()

ggplot(data = dh_sub, aes(x = Variant_type, y = segsites)) +
  # facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_col()

tt <- aggregate(x = dh_sub$segsites, by = list(dh_sub$Variant_type, dh_sub$ANN), sum)
slices <- tt$x
# lbls <- c("US", "UK", "Australia", "Germany", "France")
# lbls <- c('SW','WS', 'WWSS')
lbls <- paste(tt$Group.1, tt$Group.2)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)))

isl <- 'CZD'
eco <- 'CRAB'
dd <- dh_sub[dh_sub$ZONE==isl & dh_sub$ECOT==eco,]
dd <- dd[6:10,]
slices <- dd$segsites
# lbls <- c("US", "UK", "Australia", "Germany", "France")
# lbls <- c('SW','WS', 'WWSS')
lbls <- paste(dd$Variant_type, dd$ANN)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)), main = paste(isl, eco))

aggregate(x = dh_sub$segsites, by = list(dh_sub$ZONE, dh_sub$ECOT, dh_sub$Variant_type), sum)

isl <- 'CZD'
dd <- dh_sub[dh_sub$ZONE==isl,]
dd <- split(x = dd, f = dd$ECOT)
d1 <- merge(x = dd$WAVE_LEFT, y = dd$CRAB, by = c('Variant_type', 'ANN'))
d2 <- merge(x = dd$WAVE_LEFT, y = dd$WAVE_RIGHT, by = c('Variant_type', 'ANN'))
d3 <- merge(x = dd$CRAB, y = dd$WAVE_RIGHT, by = c('Variant_type', 'ANN'))
ggplot(data = d1, aes(x = segsites.x, y = segsites.y, col = ANN, shape = Variant_type)) +
  geom_abline(slope = 1) +
  geom_point(size = 3)
ggplot(data = d2, aes(x = segsites.x, y = segsites.y, col = ANN, shape = Variant_type)) +
  geom_abline(slope = 1) +
  geom_point(size = 3)
ggplot(data = d3, aes(x = segsites.x, y = segsites.y, col = ANN, shape = Variant_type)) +
  geom_abline(slope = 1) +
  geom_point(size = 3)
# 
# 
# 
# VCFTOOLS PI
del <- read.csv(file = 'results/Lsax_short_DEL_czs_daf_inv_findv.csv')
head(del)
oo <- data.frame(cp=unique(del$cp))
library(tidyr)
oo <- separate(data = oo, col = cp, into = c('CHROM', 'POS'), sep = '_')
head(oo)
write.table(x = oo, file = 'test/DEL_list.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
# 
# 
# 
# ANNOTATION
rm(list = ls())
ic <- unique(read.table(file = 'annotated/AN_CZD_INDEL.filt2.txt', header = TRUE))
ic$SIZE <- nchar(as.character(ic$REF)) - nchar(as.character(ic$ALT))
ic$SIZE <- abs(ic$SIZE)
head(ic)

# dt <- ic
dtp <- ic[ic$SIZE <= 50, ]
# dtp <- ic
snpeff <- read.csv(file = 'data/ANN_snpeff_classes.csv')
snpeff$impact <- as.character(snpeff$impact)
imp <- as.character(snpeff$impact)

dtp$IMP <- NA
dtp$FUN <- NA
dtp$CAT <- NA
for (i in 1:nrow(dtp)) {
  stsp <- strsplit(x = as.character(dtp$ANN[i]), split = '\\|')[[1]]
  reg <- stsp[2]
  fim <- stsp[3]
  uim <- paste(unique(stsp[stsp %in% imp]), collapse = ':')
  dtp$IMP[i] <- uim
  dtp$FUN[i] <- reg
  dtp$CAT[i] <- fim
}
head(dtp)
table(dtp$IMP)
table(dtp$CAT)

dtp$COD <- ifelse(test = dtp$CAT=='HIGH' | dtp$CAT=='MODERATE', yes = 1, no = 0)
table(dtp$COD)

write.table(x = dtp, file = 'annotated/noncod_cod/AN_CZD_INDEL.filt2.txt', quote = FALSE, sep = '\t', row.names = FALSE,
            col.names = TRUE)


