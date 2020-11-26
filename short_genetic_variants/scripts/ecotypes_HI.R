rm(list = ls())
library(data.table)
library(ggplot2)

HI <- read.csv(file = "data/HI_het_embryos_20190102.csv")
head(HI)
spa <- as.data.frame(rbindlist(lapply(list.files(path = "data", pattern = "spatial", full.names = TRUE), read.csv)))
head(spa)
# lapply(spa, head)
HI <- merge(x = HI, y = spa, by = "snail_ID")
HI$ZONE <- substr(x = HI$snail_ID, start = 1, stop = 3)

ggplot(data = HI, aes(x = ShoreH, y = HI_mod)) +
  facet_wrap(~ZONE) +
  geom_point()

ggplot(data = HI, aes(x = ShoreH, y = HI)) +
  facet_wrap(~ZONE) +
  geom_point()

ggplot(data = HI, aes(x = LCmeanDist, y = HI)) +
  facet_wrap(~ZONE) +
  geom_point()

HI$pal <- "black"
HI$pal <- ifelse(test = HI$HI < 0.25, yes = "red", no = HI$pal)
HI$pal <- ifelse(test = HI$HI > 0.5, yes = "blue", no = HI$pal)

hl <- data.frame(ZONE = c("CZA", "CZB", "CZD"), L = c(200,
                                                      mean(HI$LCmeanDist[HI$ZONE=="CZB" & HI$pal=="red"]),
                                                      mean(HI$LCmeanDist[HI$ZONE=="CZD" & HI$pal=="red"])))

hip <- ggplot(data = HI, aes(x = LCmeanDist, y = HI)) +
  facet_grid(rows = vars(ZONE), switch = "y") +
  geom_point(col = HI$pal, size = 2.5) +
  # geom_vline(data = hl, aes(xintercept = L), linetype = "dashed") +
  scale_y_continuous(position = "right") +
  ylab("Hybrid index\n") +
  xlab("Position along transect (m)") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        # axis.title.y = element_text(margin = margin(l = 20)),
        strip.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'white', color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
hip
ggsave(filename = '/Volumes/Seagate Remote Backup/3.indels/manuscript/figures/HI.pdf',
       plot = hip, height = 6, scale = 0.75, dpi = "screen")


CRAB <- HI[HI$HI<0.25, ]
CRAB$ECOT <- 'CRAB'
CRAB$SIDE <- 'C'
WAVE <- HI[HI$HI>0.5, ]
WAVE$ECOT <- 'WAVE'
head(WAVE)
zn <- unique(WAVE$ZONE)
lr <- c(200, 100, 100)
WAVE$SIDE <- unlist(lapply(X = seq_along(zn), FUN = function(x) {
  ifelse(test = WAVE[WAVE$ZONE==zn[x], "LCmeanDist"] < lr[x], yes = paste(WAVE$ECOT, "LEFT", sep = "_"),
         no = paste(WAVE$ECOT, "RIGHT", sep = "_"))
}))

CW <- rbind(CRAB,WAVE)
head(CW)


ggplot(data = CW, aes(x = LCmeanDist, y = HI, col = SIDE)) +
  facet_grid(rows = vars(ZONE)) +
  geom_point()

wdt <- split(CW, CW$ZONE)

lapply(unique(as.character(CW$ZONE)), FUN = function(x) {
  write.table(x = wdt[[x]][wdt[[x]]$ECOT=='CRAB', 1], file = paste0('Littorina_saxatilis/short_genetic_variants/CRAB_', x, '.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(x = wdt[[x]][wdt[[x]]$SIDE=='WAVE_LEFT', 1], file = paste0('Littorina_saxatilis/short_genetic_variants/WAVE_LEFT_', x, '.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(x = wdt[[x]][wdt[[x]]$SIDE=='WAVE_RIGHT', 1], file = paste0('Littorina_saxatilis/short_genetic_variants/WAVE_RIGHT_', x, '.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
})
