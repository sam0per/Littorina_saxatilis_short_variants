rm(list = ls())

HI <- read.csv(file = "data/HI_het_embryos_20190102.csv")
head(HI)
HI$ZONE <- substr(x = HI$snail_ID, start = 1, stop = 3)

ggplot(data = HI, aes(x = ShoreH, y = HI_mod)) +
  facet_wrap(~ZONE) +
  geom_point()

ggplot(data = HI, aes(x = ShoreH, y = HI)) +
  facet_wrap(~ZONE) +
  geom_point()

CRAB <- HI[HI$HI<0.25, ]
CRAB$ECOT <- 'CRAB'
WAVE <- HI[HI$HI>0.5, ]
WAVE$ECOT <- 'WAVE'
CW <- rbind(CRAB,WAVE)
head(CW)

ggplot(data = CW, aes(x = ShoreH, y = HI)) +
  facet_wrap(~ZONE) +
  geom_point()

wdt <- split(CW, CW$ZONE)
lapply(unique(as.character(CW$ZONE)), FUN = function(x) {
  write.table(x = wdt[[x]][,1], file = paste0('Littorina_saxatilis/short_genetic_variants/ecotypes_', x, '.txt'),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
})
