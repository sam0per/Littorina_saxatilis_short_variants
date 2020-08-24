rm(list=ls())

# Suppose we sample 5 individuals, then there are 10 chromosomes
# The SFS is usually expressed in terms of counts of the derived allele
# since variants must be polymorphic, this means 9 possible classes: 1 - 9 copies of the derived allele

# simple example of SNP data

snp <- c(1000,200,100,50,40,30,20,10,20) # counts of SNPs in each of the 9 classes

# and indels

indel <- c(180,30,15,10,8,8,8,9,12) # counts of indels in each of the 9 classes

# Nielsen et al Test 1 compares the likelihood of the indel data, given the frequencies expected from the SNP data
# to the likelihood of the indel data, based on the indel frequencies
# each likelihood is just a multinomial probability

snp_p <- snp/sum(snp)   # ML frequencies of the classes based on SNP data
indel_p <- indel/sum(indel)   # ML frequencies of the classes based on indel data


LL0 <- sum(indel*log(snp_p)) # log-likelihood given expectations based on SNPs

LL1 <- sum(indel*log(indel_p)) # log-likelihood given expectations based on indels

chi_app <- 2*(LL1-LL0) # chi-square approximation with n_classes-1 = 8 df

# can examine the contributions of different classes to the overall chisquare

class_contrib <- rep(0,9)
# i=1
for (i in 1:9){class_contrib[i] <- 2*((indel[i]*log(indel_p[i])+sum(indel[-i])*log(sum(indel_p[-i])))-((indel[i]*log(snp_p[i]))+sum(indel[-i])*log(sum(snp_p[-i]))))}

plot(1:9,class_contrib)


# alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)

joint_p <- (snp+indel)/sum(snp+indel)

LL0 <- sum(indel*log(joint_p)) + sum(snp*log(joint_p))

LL1 <- sum(indel*log(indel_p)) + sum(snp*log(snp_p))

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)

a <- rep(0,9)
for (i in 1:9){a[i] <- indel[i]*log(indel_p[i]) + snp[i]*log(snp_p[i]) + sum(indel[-i])*log(sum(indel_p[-i])) + sum(snp[-i])*log(sum(snp_p[-i]))}
b <- rep(0,9)
for (i in 1:9){b[i] <- (indel[i]*log(joint_p[i])) + (snp[i]*log(joint_p[i])) + sum(indel[-i])*log(sum(joint_p[-i])) + sum(snp[-i])*log(sum(joint_p[-i]))}

class_contrib <- rep(0,9)
for (i in 1:9){class_contrib[i] <- 2*(a[i]-b[i])}
plot(1:9,class_contrib)

# 
# 
# 

library(tools)
library(data.table)
library(ggplot2)

(count_fl <- list.files(path = 'summary/allele_count', pattern = '_count.txt', full.names = TRUE))
(miss_fl <- list.files(path = 'summary/allele_count', pattern = '_misscount.txt', full.names = TRUE))

make_dt <- function(vec_fl) {
  one_dt <- lapply(vec_fl, function(x) {
    read_dt <- read.table(x, sep = '\t')
    parts <- strsplit(file_path_sans_ext(basename(x)), split = "_")[[1]]
    
    if (ncol(read_dt) > 3) {
      colnames(read_dt) <- c('CHROM', 'POS', parts[length(parts)], 'TOT_INDV')
    } else {
      colnames(read_dt) <- c('CHROM', 'POS', parts[length(parts)])
    }
    read_dt$cp <- paste(read_dt$CHROM, read_dt$POS, sep = '_')
    read_dt$ZONE <- parts[2]
    if (parts[3] == 'WAVE') {
      read_dt$ECOT <- paste(parts[3], parts[4], sep = '_')
      read_dt$VTYPE <- parts[5]
    } else {
      read_dt$ECOT <- parts[3]
      read_dt$VTYPE <- parts[4]
    }
    
    return(unique(read_dt))
  })
  
  return(one_dt)
}
# in the column 'count' is recorded the number of alternative alleles.
count_dt <- make_dt(vec_fl = count_fl)
lapply(count_dt, head)
lapply(count_dt, nrow)

# in the column 'misscount' is recorded the number of individuals with missing genotypes.
miss_dt <- make_dt(vec_fl = miss_fl)
lapply(miss_dt, head)
lapply(miss_dt, nrow)

all_dt <- lapply(1:length(count_fl), function(x) {
  mg <- merge(x = count_dt[[x]], y = miss_dt[[x]])
  mg$EFF_INDV <- mg$TOT_INDV - mg$misscount
  return(mg)
})
lapply(all_dt, nrow)
lapply(all_dt, head)

# filter variants by number of individuals with called genotypes 'EFF_INDV'.
filt_indv <- function(lsdata, strcol, fltr) {
  invisible(lapply(X = 1:length(lsdata), FUN = function(x) {
    filt_one <- lsdata[[x]][lsdata[[x]][, strcol] > fltr, ]
    zn <- unique(lsdata[[x]]$ZONE)
    et <- unique(lsdata[[x]]$ECOT)
    vt <- unique(lsdata[[x]]$VTYPE)
    cat(zn, et, '\nBefore', fltr, 'filter:', nrow(unique(lsdata[[x]])), vt,
        '\nAfter', fltr, 'filter:', nrow(unique(filt_one)), vt, '\n')
    # hist(filt_one[, strcol], main = paste(zn, et, vt))
    return(filt_one)
  }))
}
filt_dt <- as.data.frame(rbindlist(filt_indv(lsdata = all_dt, strcol = 'EFF_INDV', fltr = 20)))
filt_dt$cp <- as.character(filt_dt$cp)
filt_dt$VTYPE <- as.character(filt_dt$VTYPE)
head(filt_dt)
filt_uni <- unique(filt_dt[, c('cp', 'VTYPE')])
table(filt_uni$VTYPE)
# length(unique(filt_uni[filt_uni$VTYPE=='INDEL', 'cp']))
# length(unique(filt_uni[filt_uni$VTYPE=='SNP', 'cp']))

df <- read.csv(file = 'results/Lsax_short_var_czs_daf_inv.csv')
df_uni <- unique(df[, c('cp', 'VTYPE')])
table(df_uni$VTYPE)
head(df)

dt <- merge(x = df, y = filt_dt)
head(dt)
write.csv(x = dt, file = 'results/Lsax_short_var_czs_daf_inv_findv.csv', quote = FALSE, row.names = FALSE)
dt_uni <- unique(dt[, c('cp', 'VTYPE')])
table(dt_uni$VTYPE)

nrow(dt[dt$DAF==0 | dt$DAF==1, ])
dtp <- dt[dt$DAF!=0 & dt$DAF!=1, ]
nrow(dtp)
dtp_uni <- unique(dtp[, c('cp', 'VTYPE')])
table(dtp_uni$VTYPE)
sum(dtp$DAF==0)
sum(dtp$DAF==1)
sum(dtp$DAF>1)
# sum(duplicated(dtp_uni$cp))

all_df <- as.data.frame(rbindlist(all_dt))
all_df$cp <- as.character(all_df$cp)
all_df$VTYPE <- as.character(all_df$VTYPE)
all_uni <- unique(all_df[, c('cp', 'VTYPE')])
table(all_uni$VTYPE)

# after the filter for a minimum number of individuals, the number of variants decreased slightly.
rm(list = ls())
