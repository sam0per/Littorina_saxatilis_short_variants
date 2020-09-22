dt <- unique(read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv"))
head(dt)
dti <- dt[dt$VTYPE=='INDEL', ]
table(dti$NE_W_Lcomp)
dti$SIZE <- nchar(as.character(dti$REF)) - nchar(as.character(dti$ALT))
head(dti)
hist(dti$SIZE, breaks = 50)

dti$CLASS <- NA
dti[dti$NE_W_Lcomp == 'ref_anc:ref_anc', 'CLASS'] <- ifelse(test = dti[dti$NE_W_Lcomp == 'ref_anc:ref_anc', 'SIZE'] < 0,
                                                            yes = 'INS', no = 'DEL')
head(dti[dti$CLASS == 'DEL', ])
head(dti[dti$NE_W_Lcomp == 'alt_anc:alt_anc', ])
dti[dti$NE_W_Lcomp == 'alt_anc:alt_anc', 'CLASS'] <- ifelse(test = dti[dti$NE_W_Lcomp == 'alt_anc:alt_anc', 'SIZE'] > 0,
                                                            yes = 'INS', no = 'DEL')

table(dti$CLASS)

library(ggplot2)
ggplot(data = dti, aes(x = DAF, fill = CLASS)) +
  geom_histogram(position = 'dodge', binwidth = 0.05)


write.table(x = dti, file = 'results/Lsax_short_ins_del_czs_daf_inv_findv.csv', quote = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)

dt <- unique(read.csv(file = 'results/Lsax_short_ins_del_czs_daf_inv_findv.csv'))
dt$VTYPE <- dt$CLASS
dt$CLASS <- NULL
head(dt)
write.table(x = dt, file = 'results/Lsax_short_ins_del_czs_daf_inv_findv.csv', quote = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)
# 
# 
# 
dt <- unique(read.csv(file = "results/Lsax_short_var_czs_daf_inv_findv.csv"))
dts <- dt[dt$VTYPE=='SNP', ]
write.table(x = dts, file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv', quote = FALSE, sep = ',',
            row.names = FALSE, col.names = TRUE)
