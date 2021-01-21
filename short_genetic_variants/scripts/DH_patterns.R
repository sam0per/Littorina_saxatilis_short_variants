rm(list = ls())

dh_res <- read.csv('results/DH_estimates.csv')
head(dh_res)
str(dh_res)
dh_res$ECOT <- factor(dh_res$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))
levels(dh_res$Variant_type)
# dh_res$Variant_type <- factor(dh_res$Variant_type, levels = c("INDEL", "DELETION", "INSERTION", "SNP"))
levels(dh_res$ANN)
# dh_res$ANN <- factor(dh_res$ANN, levels = c("nongenic", "syn", "nonsyn", "INDEL", "SNP"))
# dh_res$ANN <- factor(dh_res$ANN, levels = c("all", "noncoding", "nongenic", "coding", "inframe",
#                                             "synonymous", "nonsynonymous", "frameshift",
#                                             "WW", "SW", "WS", "SS"))
# 
# 
# 
aggregate(x = dh_res[, c("D", "H")], by = list(dh_res$ZONE, dh_res$ECOT, dh_res$Variant_type, dh_res$ANN),
          FUN = mean)

vt <- c('DEL', 'INS')
vt <- c('INDEL', 'SNP')
vt <- c('SNP')
vt <- c('SW', 'WS', 'WWSS')

library(RColorBrewer)
# vpal <- c("#1B9E77", "#666666")
# display.brewer.pal(n = 8, name = "Dark2")
# vpal <- brewer.pal(n = 6, name = "Dark2")[4:5]
vpal <- brewer.pal(n = 3, name = "Dark2")
# display.brewer.pal(n = 12, name = "Paired")
# vpal <- brewer.pal(n = 12, name = "Paired")[5:8]
# an <- levels(dh_res$ANN)[1:6]
an <- c('noncoding', 'coding')
# an <- intersect(as.character(unique(dh_res[dh_res$Variant_type==vt[1], 'ANN'])),
#                 as.character(unique(dh_res[dh_res$Variant_type==vt[2], 'ANN'])))
# an <- levels(dh_res$ANN)[8:length(levels(dh_res$ANN))]
# an <- levels(dh_res$ANN)[9:10]
# an <- 'coding'

dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN %in% an, ]
table(dh_sub$Variant_type)
table(dh_sub$ANN)
# dh_sub[dh_sub$ANN=='nongenic',]

head(dh_sub)
aggregate(x = dh_sub$segsites, by = list(dh_sub$ZONE, dh_sub$ECOT, dh_sub$ANN), function(x) round(c(x[1], x[2], x[1]/x[2]), 2))
id <- aggregate(x = dh_sub$segsites, by = list(dh_sub$ZONE, dh_sub$ECOT), function(x) c(x[1], x[2], x[3], x[4]))
id <- data.frame(Population = paste(id$Group.1, id$Group.2), id$x)
colnames(id) <- c("Population", "Coding deletions", "Coding insertions", "Non-coding deletions", "Non-coding insertions")
colnames(id) <- c("Population", "Coding INDELs", "Coding SNPs", "Non-coding INDELs", "Non-coding SNPs")
apply(X = aggregate(x = dh_sub$segsites, by = list(dh_sub$ZONE, dh_sub$ECOT), function(x) cbind(x[1], x[2], x[3], x[4]))$x,
      MARGIN = 1, FUN = function(y) {
        ct <- data.frame(rbind(y[1], y[2]), rbind(y[3], y[4]))
        list(expected=round(chisq.test(ct)$expected),
             observed=round(chisq.test(ct)$observed))
      })
cqt <- apply(X = id[,-1], MARGIN = 1, FUN = function(y) {
  ct <- data.frame(rbind(y[1], y[2]), rbind(y[3], y[4]))
  cq <- chisq.test(ct)
  return(cq)
})
id
# i <- 1
for (i in seq_along(cqt)) {
  id[i, -1] <- paste0(id[i, -1], " (", round(as.vector(cqt[[i]]$expected)), ")")
  id[i, "X-squared"] <- round(cqt[[i]]$statistic, 3)
}
id

dh_sub[dh_sub$segsites==32,]
dh_sub[dh_sub$segsites==29,]
dh_sub[dh_sub$segsites==1332,]
dh_sub[dh_sub$segsites==1445,]

library(ggplot2)
# DHp <- ggplot(data = dh_sub, aes(x = H, y = D, col = ANN, shape = Variant_type)) +
#   facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = vpal) +
#   labs(col = '', shape = '') +
#   theme(axis.text = element_text(size = 11),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 10),
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = '#91bfdb', color = 'black'),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))
# DHp
# ggsave(filename = 'figures/DH_var_ann_czs.pdf', plot = DHp, width = 8, height = 6, dpi = "screen")
dad <- dh_sub[, 1:6]
colnames(dad)[length(colnames(dad))] <- "DH"
dah <- dh_sub[, c(1:5,8)]
colnames(dah)[length(colnames(dah))] <- "DH"
da <- rbind(cbind(dad, S='D'),
            cbind(dah, S='H'))
da$ANN <- factor(x = da$ANN, levels = c('noncoding', 'coding'))
# da$ZE <- paste(da$ZONE, da$ECOT, da$Variant_type, sep = ' ')
# table(da$ZE)
# head(da)

# DHp <- ggplot(data = da, aes(x = ANN, y = DH, group = ZE, col = Variant_type)) +
#   facet_wrap(~S) +
#   # geom_point(aes(shape = ZE), size = 3) +
#   geom_point() +
#   geom_line() +
#   # geom_line(aes(linetype = ZE)) +
#   # scale_color_manual(values = vpal) +
#   labs(x = '') +
#   theme(axis.text = element_text(size = 17),
#         axis.title = element_text(size = 17),
#         strip.text = element_text(size = 15),
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = 'white', color = 'black'),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))
# DHp
da$ZE <- paste(da$ZONE, da$ECOT, da$ANN, sep = ' ')
# da$Variant_type <- factor(x = da$Variant_type, levels = c("SW", "WWSS", "WS"))
library(Rmisc)
dm <- aggregate(x = da$DH, by = list(S = da$S, Variant_type = da$Variant_type, ANN = da$ANN), function(x) round(CI(x),2))
dm
dm <- as.data.frame(cbind(dm[,1:3], DH = dm$x[,2]))
dm$ZE <- paste(dm$Variant_type, dm$ANN, sep = ' ')
dm

DHp <- ggplot(data = da, aes(x = Variant_type, y = DH, group = ZE, col = ANN)) +
  facet_wrap(~S) +
  # geom_point(aes(shape = ZE), size = 3) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(size = 2, alpha = 0.7) +
  geom_point(data = dm, aes(x = Variant_type, y = DH), col = "black", size = 7) +
  geom_point(data = dm, aes(x = Variant_type, y = DH, col = ANN), size = 5) +
  # geom_line(aes(linetype = ZE)) +
  scale_color_manual(values = vpal) +
  labs(x = '', y = '', col = '') +
  theme(axis.text = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.position = 'top',
        # axis.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  guides(col = guide_legend(override.aes = list(size = 4)))
DHp
ggsave(filename = paste('figures/DH', paste(an, collapse = '_'), paste(vt, collapse = '_'),'czs.pdf', sep = '_'),
       plot = DHp, scale = 0.7, dpi = "screen")

library(Rmisc)
round(CI(da[da$Variant_type=='INDEL' & da$ANN=='coding' & da$S == 'D', 'DH']), 2)
round(CI(da[da$Variant_type=='INDEL' & da$ANN=='noncoding' & da$S == 'D', 'DH']), 2)
round(CI(da[da$Variant_type=='INDEL' & da$ANN=='coding' & da$S == 'H', 'DH']), 2)
round(CI(da[da$Variant_type=='INDEL' & da$ANN=='noncoding' & da$S == 'H', 'DH']), 2)

round(CI(da[da$Variant_type=='SNP' & da$ANN=='coding' & da$S == 'D', 'DH']), 2)
round(CI(da[da$Variant_type=='SNP' & da$ANN=='noncoding' & da$S == 'D', 'DH']), 2)
round(CI(da[da$Variant_type=='SNP' & da$ANN=='coding' & da$S == 'H', 'DH']), 2)
round(CI(da[da$Variant_type=='SNP' & da$ANN=='noncoding' & da$S == 'H', 'DH']), 2)

# ggsave(filename = 'figures/DH_gBGC_czs.pdf', plot = DHp, width = 8, height = 6)
#
# hist(dh_sub$D)
# summary(lm(formula = D~Variant_type, data = dh_sub, weights = sqrt(segsites)))
# summary(lm(formula = H~Variant_type, data = dh_sub, weights = sqrt(segsites)))
#
# ggplot(data = dh_sub, aes(x = Variant_type, y = D)) +
#   geom_point()
# ggplot(data = dh_sub, aes(x = Variant_type, y = H)) +
#   geom_point()

# BETWEEN
lapply(X = an, FUN = function(x) {

  dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN==x, ]

  Dlm <- summary(lm(formula = D~-1+Variant_type, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~-1+Variant_type, data = dh_sub, weights = sqrt(segsites)))

  return(list(x, Dlm, Hlm))

})

lapply(X = an, FUN = function(x) {

  dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN==x, ]

  Dlm <- summary(lm(formula = D~Variant_type*ECOT, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~Variant_type*ECOT, data = dh_sub, weights = sqrt(segsites)))

  return(list(x, Dlm, Hlm))

})

library(data.table)
dh_vt <- lapply(X = an, FUN = function(x) {

  dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN==x, ]

  Dlm <- summary(lm(formula = D~-1 + Variant_type, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~-1 + Variant_type, data = dh_sub, weights = sqrt(segsites)))

  # return(list(x, Dlm, Hlm))
  op <- list(as.data.frame(cbind(Par = row.names(coef(Dlm)), coef(Dlm), ANN = x, SS = 'D')),
             as.data.frame(cbind(Par = row.names(coef(Hlm)), coef(Hlm), ANN = x, SS = 'H')))

  return(op)

})
dh_vt <- rbind(data.frame(rbindlist(dh_vt[[1]])),
               data.frame(rbindlist(dh_vt[[2]])),
               data.frame(rbindlist(dh_vt[[3]])))

dh_vt <- rbind(data.frame(rbindlist(dh_vt[[1]])),
               data.frame(rbindlist(dh_vt[[2]])),
               data.frame(rbindlist(dh_vt[[3]])),
               data.frame(rbindlist(dh_vt[[4]])),
               data.frame(rbindlist(dh_vt[[5]])),
               data.frame(rbindlist(dh_vt[[6]])))
str(dh_vt)
# dh_vt$Pal <- ifelse(test = dh_vt$Par == levels(dh_vt$Par)[1], yes = 'INDELs', no = 'SNPs')
dh_vt$Pal <- substr(x = as.character(dh_vt$Par), start = 13, stop = nchar(as.character(dh_vt$Par)))
levels(dh_vt$ANN)
dh_vt$ANN <- factor(x = dh_vt$ANN, levels = c('nongenic', 'syn', 'nonsyn'))

dh_vt$Estimate <- as.numeric(as.character(dh_vt$Estimate))
dh_vt$Std..Error <- as.numeric(as.character(dh_vt$Std..Error))
btw <- ggplot(data = dh_vt, aes(x = ANN, y = Estimate, col = Pal)) +
  facet_grid(cols = vars(SS)) +
  geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error), width=.1, size = 1) +
  geom_point(size = 4) +
  # scale_color_manual(values = c("#1B9E77", "#666666")) +
  scale_color_manual(values = vpal) +
  labs(col = '', x = '', y = 'Fitted value') +
  theme(legend.text = element_text(size = 16),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 320, hjust = 0),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 14),
        panel.background = element_blank(),
        # strip.background = element_rect(fill = '#91bfdb', color = 'black'),
        strip.background = element_rect(fill = 'white', color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
btw
ggsave(filename = paste('figures/DH_noncoding', paste(vt, collapse = '_'),'fit.pdf', sep = '_'),
       plot = btw, width = 8, height = 6, dpi = "screen")
ggsave(filename = paste('figures/DH_annotated', paste(vt, collapse = '_'),'fit.pdf', sep = '_'),
       plot = btw, width = 8, height = 6, dpi = "screen")
ggsave(filename = 'figures/DH_between_variants.pdf', plot = btw, width = 8, height = 6, dpi = "screen")
#
#
#
# rm(list = setdiff(ls(), c('dh_res', 'an', 'vt')))
summary(lm(formula = D ~ Variant_type * ANN, data = dh_sub, weights = sqrt(segsites)))
summary(lm(formula = D ~ -1 + Variant_type + ANN, data = dh_sub, weights = sqrt(segsites)))
dh_one <- dh_sub[dh_sub$Variant_type=='WWSS', ]
dh_one <- dh_one[dh_one$ECOT=='CRAB', ]
dh_one <- split(x = dh_one, f = dh_one$ZONE)

ggplot(data = dh_one, aes(x = ANN, y = D, col = ZONE)) +
  geom_point(size = 3) +
  labs(col = '', x = '') +
  theme(axis.text = element_text(size = 18), legend.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))
ggplot(data = dh_one, aes(x = ANN, y = H, col = ZONE)) +
  geom_point(size = 3) +
  labs(col = '', x = '') +
  theme(axis.text = element_text(size = 18), legend.text = element_text(size = 16),
        axis.title.y = element_text(size = 18))

summary(lm(formula = D ~ -1 + ANN, data = dh_one$CZA, weights = sqrt(segsites)))
summary(lm(formula = D ~ -1 + ANN, data = dh_one$CZB, weights = sqrt(segsites)))
summary(lm(formula = D ~ -1 + ANN, data = dh_one$CZD, weights = sqrt(segsites)))
# WITHIN INDELS AND SNPS
lapply(X = vt, FUN = function(x) {

  dh_sub <- dh_res[dh_res$Variant_type==x & dh_res$ANN %in% an, ]

  Dlm <- summary(lm(formula = D~-1 + ANN, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~-1 + ANN, data = dh_sub, weights = sqrt(segsites)))

  return(list(x, Dlm, Hlm))
  # return(list(coef(Dlm), coef(Hlm)))

})
# -0.14284723+0.05209618
dh_vt <- lapply(X = vt, FUN = function(x) {

  dh_sub <- dh_res[dh_res$Variant_type==x & dh_res$ANN %in% an, ]

  Dlm <- summary(lm(formula = D~-1 + ANN, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~-1 + ANN, data = dh_sub, weights = sqrt(segsites)))

  op <- list(as.data.frame(cbind(Par = row.names(coef(Dlm)), coef(Dlm), VT = x, SS = 'D')),
             as.data.frame(cbind(Par = row.names(coef(Hlm)), coef(Hlm), VT = x, SS = 'H')))

  return(op)

})

dh_vt <- rbind(data.frame(rbindlist(dh_vt[[1]])),
               data.frame(rbindlist(dh_vt[[2]])))
str(dh_vt)
dh_vt$Pal <- substr(x = as.character(dh_vt$Par), start = 4, stop = nchar(as.character(dh_vt$Par)))
dh_vt$Estimate <- as.numeric(as.character(dh_vt$Estimate))
dh_vt$Std..Error <- as.numeric(as.character(dh_vt$Std..Error))
dh_vt$Pal <- factor(x = dh_vt$Pal, levels = c('nongenic', 'synonymous', 'nonsynonymous'))

btw <- ggplot(data = dh_vt, aes(x = Pal, y = Estimate, col = VT)) +
  facet_grid(cols = vars(SS)) +
  geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error), width=.1, size = 1) +
  geom_point(size = 3) +
  # scale_color_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  scale_color_manual(values = c("#1B9E77", "#666666")) +
  labs(col = '', x = '') +
  theme(legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        strip.background = element_rect(fill = '#91bfdb', color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
btw
ggsave(filename = 'figures/DH_within_variants.pdf', plot = btw, width = 8, height = 6, dpi = "screen")
