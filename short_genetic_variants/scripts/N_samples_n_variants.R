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

dt <- read.csv(file = 'results/Lsax_short_var_czs_daf_inv_findv.csv')
head(dt)
dt_uni <- unique(dt[, c('cp', 'VTYPE')])
table(dt_uni$VTYPE)
length(unique(dt_uni$cp))
dt$CAT <- paste(dt$ZONE, dt$ECOT, sep = ':')
ds <- split(dt, f = dt$CAT)
# 69 Littorina_saxatilis/short_genetic_variants/CRAB_CZA.txt
# 62 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZA.txt
# 23 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZA.txt

# 64 Littorina_saxatilis/short_genetic_variants/CRAB_CZB.txt
# 59 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZB.txt
# 45 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZB.txt

# 67 Littorina_saxatilis/short_genetic_variants/CRAB_CZD.txt
# 31 Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_CZD.txt
# 73 Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_CZD.txt
lapply(ds, function(x) table(x$TOT_INDV))
lapply(ds, function(x) table(x$EFF_INDV))

(ctg <- levels(as.factor(dt$CAT))[9])
head(ds[[ctg]])
# sum(ds[[ctg]]$DAF==0)
ds_uni <- unique(ds[[ctg]][, c('cp', 'VTYPE', 'EFF_INDV')])
tb <- table(ds_uni$EFF_INDV, ds_uni$VTYPE)
tb <- as.data.frame(cbind(tb, TOT=as.vector(table(ds_uni$EFF_INDV))))
c1 <- row.names(tb)[(nrow(tb)-5):nrow(tb)]
ot <- rbind(tb[nrow(tb), ] + tb[nrow(tb)-1, ] + tb[nrow(tb)-2, ] + tb[nrow(tb)-3, ] + tb[nrow(tb)-4, ] + tb[nrow(tb)-5, ],
            tb[nrow(tb), ] + tb[nrow(tb)-1, ] + tb[nrow(tb)-2, ] + tb[nrow(tb)-3, ] + tb[nrow(tb)-4, ],
            tb[nrow(tb), ] + tb[nrow(tb)-1, ] + tb[nrow(tb)-2, ] + tb[nrow(tb)-3, ],
            tb[nrow(tb), ] + tb[nrow(tb)-1, ] + tb[nrow(tb)-2, ],
            tb[nrow(tb), ] + tb[nrow(tb)-1, ],
            tb[nrow(tb), ])
# ot <- rbind(tb[nrow(tb), ] + tb[nrow(tb)-1, ] + tb[nrow(tb)-2, ],
#             tb[nrow(tb), ] + tb[nrow(tb)-1, ],
#             tb[nrow(tb), ])
ot <- data.frame(CAT = ctg, N = c1, ot)
ot
# of <- ot
(of <- rbind(of, ot))
# nrow(ds_uni[ds_uni$EFF_INDV==21,])
# sum(table(ds_uni$EFF_INDV, ds_uni$VTYPE)[,1])
# sum(table(ds_uni$EFF_INDV, ds_uni$VTYPE)[,2])
# length(ds_uni[ds_uni$VTYPE=='INDEL', 'cp'])
# length(ds_uni[ds_uni$VTYPE=='SNP', 'cp'])
# length(unique(ds_uni[ds_uni$VTYPE=='INDEL', 'cp']))
# length(unique(ds_uni[ds_uni$VTYPE=='SNP', 'cp']))
of[of$CAT!='CZA:WAVE_RIGHT', 'X'] <- rep(c('I', 'J', 'K', 'L', 'M', 'N'), 8)
of[of$CAT=='CZA:WAVE_RIGHT', 'X'] <- c('L', 'M', 'N')
str(of)
of$CAT <- as.character(of$CAT)
of$N <- as.character(of$N)
of <- rbind(of, c('ALL', 1000, as.integer(table(dt_uni$VTYPE)), length(unique(dt_uni$cp)), 'A'))
of$INDEL <- as.integer(of$INDEL)
of$SNP <- as.integer(of$SNP)
of$TOT <- as.integer(of$TOT)

write.csv(x = of, file = 'results/Lsax_short_var_czs_N_samples_n_variants.csv', quote = FALSE, row.names = FALSE)
# edit in excel
of <- read.csv(file = 'results/Lsax_short_var_czs_N_samples_n_variants.csv')

library(ggplot2)
# data.frame(INDEL="#1B9E77", SNP="#666666")
N_ind <- ggplot(data = of, aes(x = X, y = INDEL, alpha = CAT)) +
  geom_col(position = 'dodge', col = 'black', fill = '#1B9E77') +
  scale_x_discrete(labels = c('T', 'N-5', 'N-4', 'N-3','N-2', 'N-1', 'N')) +
  labs(x = '', alpha = '') +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
N_ind
ggsave(filename = 'figures/N_samples_n_INDEL.pdf', plot = N_ind)
N_snp <- ggplot(data = of, aes(x = X, y = SNP, alpha = CAT)) +
  geom_col(position = 'dodge', col = 'black', fill = '#666666') +
  scale_x_discrete(labels = c('T', 'N-5', 'N-4', 'N-3','N-2', 'N-1', 'N')) +
  labs(x = '', alpha = '') +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 16),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
N_snp
ggsave(filename = 'figures/N_samples_n_SNP.pdf', plot = N_snp)

# choose N-3
rm(list = ls())
dt <- read.csv(file = 'test/GM_CZA_CRAB_INDEL.filt2.66N.csv')
head(dt)
dt <- dt[dt$CHROM!='CHROM', ]
franc_dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv")
head(franc_dt)

dta <- merge(x = dt, y = franc_dt[, c('cp', 'ZONE', 'VTYPE', 'ECOT', )])

franc_dt$DAF <- NA
franc_dt$DAF <- ifelse(test = grepl(pattern = "alt", x = franc_dt$NE_W_Lcomp), yes = franc_dt$RFREQ, no = franc_dt$AFREQ)
sample_n(tbl = franc_dt, size = 30)