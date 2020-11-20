# INDEL SIZE
# NONCODING
rm(list = ls())
(ff <- list.files(path = 'results/marker_density', pattern = '_noncoding_DEL', full.names = TRUE))

library(data.table)
da <- data.frame(rbindlist(lapply(X = ff, FUN = read.table, header = TRUE)))
ue <- unique(as.character(da$cp))
colnames(da)
head(da)
da <- unique(da[, c('cp', 'SIZE', 'COD')])
median(da$SIZE)
mean(da$SIZE)
# das <- split(da, da$COD)
table(da$COD)
# load library
library(ggplot2)

da$category <- ifelse(test = da$SIZE>6, yes = "6+", no = as.character(da$SIZE))
# da$category <- ifelse(test = da$SIZE>9, yes = "9+", no = as.character(da$SIZE))
table(da$category)
data <- data.frame(table(da$category))
colnames(data) <- c('category', 'count')
data
# Compute percentages
data$fraction <- data$count / sum(data$count)
sum(data$fraction[2:9])
sum(data$count[2:9])
sum(data$count)
# test different proportions between insertions and deletions
tt <- data.frame(twonine=c(1251,1562), tot=c(2384,2580))
chisq.test(tt)
# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# data$labelPosition <- cumsum(data$fraction) - 0.5*data$fraction

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
# Make the plot
gg <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(col = 'black', alpha = 0.8) +
  # geom_text( x=3.5, aes(y=labelPosition, label=paste(round(fraction * 100), "%")), size=6) + # x here controls label position (inner / outer)
  geom_text(data = data[c(-8:-9),], x=3.5, aes(y=labelPosition, label=paste(round(fraction * 100), "%", sep = "")), size=7) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  # scale_color_manual(values = 'black') +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = c(.5, .5), legend.title = element_text(size = 17),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  labs(fill = "Non-coding\ninsertion size")
gg
ggsave(filename = '/Volumes/Seagate Remote Backup/3.indels/manuscript/figures/noncoding_ins_size.pdf',
       plot = gg, scale = 0.8, dpi = "screen")

noncins <- data
# INDEL SIZE
# CODING
rm(list = setdiff(ls(), 'noncins'))
(ff <- list.files(path = 'results/marker_density', pattern = '_coding_INS', full.names = TRUE))

library(data.table)
da <- data.frame(rbindlist(lapply(X = ff, FUN = read.table, header = TRUE)))
ue <- unique(as.character(da$cp))
colnames(da)
head(da)
da <- unique(da[, c('cp', 'SIZE', 'COD')])
median(da$SIZE)
mean(da$SIZE)
# das <- split(da, da$COD)
table(da$COD)
# load library
library(ggplot2)

da$category <- ifelse(test = da$SIZE>6, yes = "6+", no = as.character(da$SIZE))
table(da$category)
data <- data.frame(table(da$category))
colnames(data) <- c('category', 'count')
data
if (sum(data$category==4)==0) {
  # data$category <- as.character(data$category)
  data <- rbind(data, data.frame(category=factor("4", levels = "4"),count=0))
}
data
# str(noncins)
data$category <- factor(x = data$category, levels = c("1", "2", "3", "4", "5", "6", "6+"))
data <- data[order(data$category),]
# test proportions noncoding vs coding
tt <- t(data.frame(nonc=noncins$count, cod=data$count))
chisq.test(tt, simulate.p.value = TRUE)
cq <- chisq.test(tt)
cq
cq$expected
cq$observed
# 
rr <- data.frame(C=c(tt[1,], tt[2,]), A=c(rep('nonc', 7), rep('cod', 7)), S=1:7)
rr$AS <- paste0(rr$A, rr$S)
rr$S <- as.character(rr$S)
summary(glm(formula = C~-1+A, family = 'poisson', data = rr))
summary(glm(formula = C~-1+S, family = 'poisson', data = rr))
library(caret)
varImp(glm(formula = C~S, family = 'poisson', data = rr))
# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
# data$labelPosition <- cumsum(data$fraction) - 0.5*data$fraction

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
# Make the plot
gg <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(col = 'black', alpha = 0.8) +
  # geom_text( x=3.5, aes(y=labelPosition, label=paste(round(fraction * 100), "%")), size=6) + # x here controls label position (inner / outer)
  geom_text(x=3.5, aes(y=labelPosition, label=paste(round(fraction * 100), "%", sep = "")), size=7) +
  scale_fill_brewer(palette = "PiYG", direction = 1) +
  # scale_color_manual(values = 'black') +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = c(.5, .5), legend.title = element_text(size = 17),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  labs(fill = "Coding\ninsertion size")
gg
ggsave(filename = '/Volumes/Seagate Remote Backup/3.indels/manuscript/figures/coding_ins_size.pdf',
       plot = gg, scale = 0.8, dpi = "screen")
