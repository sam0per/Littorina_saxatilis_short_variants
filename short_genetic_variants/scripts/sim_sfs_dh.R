rm(list = ls())
library(RColorBrewer)
library(ggplot2)
library(patchwork)
# D < 0 - H > 0
pu_ac <- c(12,6,3,3,2,1,1,1,1)
sfs <- data.frame(N=1:9, C=pu_ac)
sum(pu_ac)
barplot(sfs$C)
sfs$T <- 30/sum(1/sfs$N)
sfs$E <- sfs$T/sfs$N
barplot(sfs$E)

# D > 0 - H < 0
po_ac <- c(9,5,3,2,2,2,2,2,3)
sum(po_ac)
sfs <- data.frame(N=1:9, C=po_ac)
barplot(sfs$C)
sfs$T <- 30/sum(1/sfs$N)
sfs$E <- sfs$T/sfs$N
barplot(sfs$E)

# D > 0 - H > 0
ba_ac <- c(9,5,3,4,4,2,1,1,1)
sum(ba_ac)
sfs <- data.frame(N=1:9, C=ba_ac)
barplot(sfs$C)
sfs$T <- 30/sum(1/sfs$N)
sfs$E <- sfs$T/sfs$N
barplot(sfs$E)

sfs <- data.frame(N=rep(1:9,3), C=c(pu_ac, po_ac, ba_ac), S=c(rep("Purifying",9), rep("Positive",9), rep("Balancing",9)))
sfs$N <- as.factor(sfs$N)
sfs$S <- factor(x = sfs$S, levels = c("Purifying", "Positive", "Balancing"))
sfs

sim_sfs <- ggplot(data = sfs, aes(x = N, y = C)) +
  facet_wrap(facets = ~S, nrow = 1) +
  geom_col(aes(fill = S), col = 'black') +
  scale_fill_manual(values = c("blue", "green", "red")) +
  labs(x = 'Derived allele class', y = "count", fill = "") +
  theme(legend.position = "top",
        legend.title = element_text(size = 12), legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
sim_sfs
ggsave(filename = 'figures/sim/DAC_sim_sfs_selection.pdf',
       plot = sim_sfs, height = 5, scale = 0.75, dpi = "screen")

# sim_sfs <- ggplot(data = sfs, aes(x = N, y = C)) +
#   facet_wrap(facets = ~S, nrow = 2) +
#   geom_col(aes(fill = S), col = 'black') +
#   scale_fill_manual(values = c("red", "green", "blue")) +
#   labs(x = 'Derived allele class', y = "count", fill = "") +
#   theme(legend.position = "top",
#         legend.title = element_text(size = 12), legend.text = element_text(size = 11),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 14),
#         strip.text = element_text(size = 12),
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = 'white', color = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.line = element_line(size = 0.2, linetype = "solid",
#                                  colour = "black"),
#         panel.grid = element_line(colour = "gray70", size = 0.2))
# sim_sfs

ddh <- data.frame(D=c(0.37, -0.06, -0.18), H=c(0.11, -0.83, 0.27), S=c("Balancing", " Positive", "Purifying"),
                  A=rep("Summary statistics", 3))
ddh <- ddh[order(ddh$S), ]
sim_ss <- ggplot(data = ddh, aes(x = D, y = H, col = S)) +
  facet_wrap(~A) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4) +
  scale_color_manual(values = c("green", "red", "blue")) +
  labs(x = "Tajima's D", y = "Fay and Wu's H", col = "") +
  theme(legend.position = "none",
        legend.title = element_text(size = 12), legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
sim_f <- sim_sfs / sim_ss
ggsave(filename = 'figures/sim/DAC_sim_sfs_DH_selection.pdf', plot = sim_f, width = 9, scale = 0.75, dpi = "screen")
