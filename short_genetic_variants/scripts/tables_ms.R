# TABLE Ya: proportions of unique variants in each population
# INDEL
rm(list = ls())
(ff <- list.files(path = 'results/marker_density', pattern = 'coding', full.names = TRUE))
(ff <- ff[grepl(pattern = "INDEL", x = ff)])
(ff <- ff[grepl(pattern = "SNP", x = ff)])
library(data.table)
ucp <- rbindlist(lapply(X = ff, FUN = read.table, header = TRUE))[, 'cp']
ucp <- as.character(ucp$cp)
ccp <- data.frame(table(ucp))
head(ccp)

table(ccp$Freq)

max(ccp$Freq)
sum(ccp$Freq)

ucp <- as.character(ccp[ccp$Freq==1, 'ucp'])
ccp <- as.character(ccp[ccp$Freq>1, 'ucp'])

tt <- vector(mode = "integer", length = length(ff))
tt
pq <- matrix(nrow = length(ff), ncol = 2)
pq
for (i in seq(from = 1, to = length(ff), by = 2)) {
  # i <- 15
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  oo <- unique(c(as.character(dd$cp), as.character(ee$cp)))
  tt[i] <- paste(round(length(intersect(oo, ucp))/(nrow(dd)+nrow(ee)) * 100), '%', sep = '')
  pq[i,] <- c(length(intersect(oo, ucp)),(nrow(dd)+nrow(ee)))
  # tt[i] <- paste(round(length(setdiff(oo, ccp))/(nrow(dd)+nrow(ee)) * 100), '%', sep = '')
  # pq[i,] <- c(length(setdiff(oo, ccp)),(nrow(dd)+nrow(ee)))
}
tt[tt!="0"]
pq <- pq[complete.cases(pq),]
pqi <- pq
pqi
rm(list = setdiff(ls(), "pqi"))

tt <- vector(mode = "integer", length = length(ff))
tt
pq <- matrix(nrow = length(ff), ncol = 2)
pq
for (i in seq(from = 1, to = length(ff), by = 2)) {
  # i <- 1
  dd <- read.table(file = ff[i], header = TRUE)
  ee <- read.table(file = ff[i+1], header = TRUE)
  oo <- unique(c(as.character(dd$cp), as.character(ee$cp)))
  tt[i] <- paste(round(length(intersect(oo, ccp))/(nrow(dd)+nrow(ee)) * 100), '%', sep = '')
  pq[i,] <- c(length(intersect(oo, ccp)),(nrow(dd)+nrow(ee)))
}
tt[tt!="0"]
pq <- pq[complete.cases(pq),]
pqi <- pq
pqi
rm(list = setdiff(ls(), "pqi"))

lapply(1:9, function(x) chisq.test(rbind(pqi[x,], pq[x,])))
# chisq.test(rbind(pqi[1,], pq[1,]))
# chisq.test(rbind(pqi[2,], pq[2,]))

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
colnames(ind) <- paste(colnames(ind), 'ind', sep = '_')
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
ind$nonc_snp <- nonc
# pop <- c('CZA_C', 'CZA_WL', 'CZA_WR', 'CZB_C', 'CZB_WL', 'CZB_WR', 'CZD_C', 'CZD_WL', 'CZD_WR')

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
ind$cod_snp <- cod
tb <- ind[,-3]

apply(X = tb, MARGIN = 1, FUN = function(x) {
  ct <- data.frame(INDEL = x[1:2], SNP = x[3:4])
  chisq.test(ct)
})

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