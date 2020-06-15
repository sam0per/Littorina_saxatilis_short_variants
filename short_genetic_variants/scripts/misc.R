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
