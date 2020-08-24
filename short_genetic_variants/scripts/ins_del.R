dt <- read.csv(file = "results/Lsax_short_var_czs_daf_inv.csv")
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
