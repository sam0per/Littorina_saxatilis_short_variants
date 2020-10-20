rm(list = ls())

dh_res <- read.csv('results/DH_estimates.csv')
head(dh_res)
str(dh_res)
dh_res$ECOT <- factor(dh_res$ECOT, levels = c("WAVE_LEFT", "CRAB", "WAVE_RIGHT"))
dh_res$Variant_type <- factor(dh_res$Variant_type, levels = c("INDEL", "DELETION", "INSERTION", "SNP"))
dh_res$ANN <- factor(dh_res$ANN, levels = c("all", "noncoding", "coding", "inframe", "synonymous", "nonsynonymous", "frameshift",
                                            "WW", "SW", "WS", "SS"))

vt <- c('INDEL', 'SNP')
# vt <- c('SNP')
vpal <- c("#1B9E77", "#666666")
# display.brewer.pal(n = 12, name = "Paired")
# vpal <- brewer.pal(n = 12, name = "Paired")[5:8]
an <- intersect(as.character(unique(dh_res[dh_res$Variant_type==vt[1], 'ANN'])),
                as.character(unique(dh_res[dh_res$Variant_type==vt[2], 'ANN'])))
# an <- levels(dh_res$ANN)[8:length(levels(dh_res$ANN))]
# an <- levels(dh_res$ANN)[9:10]
# an <- 'coding'
lapply(X = an, FUN = function(x) {
  
  dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN==x, ]
  
  # library(ggplot2)
  # DHp <- ggplot(data = dh_sub, aes(x = H, y = D, col = Variant_type)) +
  #   facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  #   geom_point(size = 2) +
  #   scale_color_manual(values = vpal) +
  #   labs(col = '') +
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
  
  Dlm <- summary(lm(formula = D~Variant_type, data = dh_sub, weights = sqrt(segsites)))
  Hlm <- summary(lm(formula = H~Variant_type, data = dh_sub, weights = sqrt(segsites)))
  
  return(list(x, Dlm, Hlm))
  
  # ggplot(data = dh_sub, aes(x = Variant_type, y = D)) +
  #   geom_point() +
  #   labs(title = x)
  # ggplot(data = dh_sub, aes(x = Variant_type, y = H)) +
  #   geom_point() +
  #   labs(title = x)
  
})

dh_sub <- dh_res[dh_res$Variant_type %in% vt & dh_res$ANN %in% an, ]

library(ggplot2)
DHp <- ggplot(data = dh_sub, aes(x = H, y = D, col = Variant_type, alpha = ANN)) +
  facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_point(size = 3) +
  scale_color_manual(values = vpal) +
  labs(col = '', alpha = '') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        strip.background = element_rect(fill = '#91bfdb', color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
DHp
ggsave(filename = 'figures/DH_var_ann_czs.pdf', plot = DHp, width = 8, height = 6)

DHp <- ggplot(data = dh_sub, aes(x = H, y = D, col = ANN)) +
  facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
  geom_point(size = 3) +
  scale_color_manual(values = vpal) +
  labs(col = '') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 10),
        panel.background = element_blank(),
        strip.background = element_rect(fill = '#91bfdb', color = 'black'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
DHp
# ggsave(filename = 'figures/DH_SW_WS_czs.pdf', plot = DHp, width = 8, height = 6)
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
