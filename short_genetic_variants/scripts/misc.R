A_dir <- "/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115/"
txt_fl <- list.files(path = "CZCLI006_comp", pattern = "_CZ", full.names = TRUE)[c(1,7)]
tbl <- lapply(txt_fl, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
identical(tbl[[1]]$cp, tbl[[2]]$cp)

A_txt_fl <- list.files(path = paste0(A_dir, "CZCLI006_comp"), pattern = "_CZ", full.names = TRUE)[c(2,14)]
A_tbl <- lapply(A_txt_fl, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
identical(A_tbl[[1]]$cp, A_tbl[[2]]$cp)

head(tbl[[1]])
head(A_tbl[[1]])
cp_diff <- data.frame(cp = setdiff(A_tbl[[1]]$cp, tbl[[1]]$cp))
cp_diff <- data.frame(cp = setdiff(A_tbl[[1]][A_tbl[[1]]$sel==FALSE, "cp"], tbl[[1]][tbl[[1]]$sel==FALSE, "cp"]))
sum(A_tbl[[1]]$sel)
sum(tbl[[1]]$sel)
cp_diff <- data.frame(cp = setdiff(tbl[[1]]$cp, A_tbl[[1]]$cp))
head(cp_diff)
library(tidyr)
contig_diff <- separate(data = cp_diff, col = cp, into = c("contig", "pos"), sep = "_")
head(contig_diff)

contig_int <- read.table(file = "Littorina_saxatilis/short_genetic_variants/contig_intervals.list")
head(contig_int)
id <- round(runif(n = 1, min = 1, max = nrow(contig_diff)))
as.character(contig_int[, 1])[grepl(pattern = paste0(contig_diff$contig[id], ":"), x = as.character(contig_int[, 1]))]
contig_diff$pos[id]
contig_diff$contig[id]

# 
# 
# 

# MARKER DENSITY
tt <- mer[mer$ISL=='CZA' & mer$ECOT=='WAVE_LEFT', ]
head(tt)
tt <- tt[tt$x.x!=0 & tt$x.y!=0, ]

# perspective plot
library(MVN)
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "persp")
# contour plot
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "contour")

library(bestNormalize)
(BNobject.x <- bestNormalize(tt$prop.x))
(BNobject.y <- bestNormalize(tt$prop.y))

identical(BNobject.x$x, tt$prop.x)
identical(BNobject.x$x.t, BNobject.x$chosen_transform$x.t)

hist(tt$prop.x)
hist(tt$prop.y)
MASS::truehist(BNobject.x$chosen_transform$x.t)
MASS::truehist(BNobject.y$x.t)
result <- mvn(data.frame(BNobject.x$x.t, BNobject.y$x.t), mvnTest = "hz", multivariatePlot = "persp")
# contour plot
result <- mvn(tt[, c('prop.x', 'prop.y')], mvnTest = "hz", multivariatePlot = "contour")

sunflowerplot(tt$prop.x, tt$prop.y)

# install the package 
# install.packages("ggstatsplot")

# Load the package
library(ggstatsplot)

# Load the dataset 
data("warpbreaks")

# Create a boxplot of the dataset, outliers are shown as two distinct points
boxplot(tt$prop.x)$out

#Create a boxplot that labels the outliers  
head(tt)
# ggbetweenstats(tt, x.x, x.y, outlier.tagging = TRUE)
Q <- quantile(tt$prop.x, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(tt$prop.x)
up <- Q[2]+1.5*iqr # Upper Range  
low <- Q[1]-1.5*iqr # Lower Range
eliminated <- subset(tt, tt$prop.x > (Q[1] - 1.5*iqr) & tt$prop.x < (Q[2]+1.5*iqr))

outliers <- boxplot(tt$prop.x, plot=FALSE)$out
tt_notl <- tt[-which(tt$prop.x %in% outliers),]

# boxplot(eliminated$prop.x)$out
boxplot(tt_notl$prop.x)$out

outliers <- boxplot(tt_notl$prop.y, plot=FALSE)$out
tt_notl <- tt_notl[-which(tt_notl$prop.y %in% outliers),]

tt <- tt_notl
hist(tt$prop.x)

# 
# 
# 

library(tidyverse)
library(caret)
# Load the data
data("swiss")
# Inspect the data
sample_n(swiss, 3)
# Define training control
set.seed(123) 
train.control <- trainControl(method = "cv", number = 10)
# Train the model
model <- train(Fertility ~., data = swiss, method = "lm",
               trControl = train.control)
# Summarize the results
print(model)
train.control <- trainControl(method = "repeatedcv", 
                              number = 10, repeats = 3)
# Train the model
model <- train(Fertility ~., data = swiss, method = "lm",
               trControl = train.control)
# Summarize the results
print(model)

model$finalModel
summary(lm(formula = Fertility ~., data = swiss))

confint(model$finalModel)
confint(object = lm(formula = Fertility ~., data = swiss))
