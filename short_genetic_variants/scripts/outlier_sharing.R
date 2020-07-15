rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "reshape2", "tidyr", "tools", "data.table", "RColorBrewer", "dplyr", "textshape", "plotly",
              "devtools", "limma")
# source("https://bioconductor.org/biocLite.R")
# biocLite("snpStats")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

################################################################################################################
##### INPUT ####################################################################################################
################################################################################################################
# Get cline fits
(cl_fl <- list.files(path = "CZCLI006_comp", full.names = TRUE))
(cl_fl <- cl_fl[grep(pattern = "ANG|NoInv", x = cl_fl, invert = TRUE)])
cl_ls <- lapply(cl_fl, read.table, header = TRUE)
# lapply(cl_ls, head)

cl_dt <- as.data.frame(rbindlist(lapply(seq_along(cl_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][2]
  side <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][3]
  zone <- paste(island, side, sep = "_")
  vtype <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][4]
  # paste(zone, vtype)
  
  odt <- mutate(cl_ls[[x]], ZONE = zone, VTYPE = vtype)
  sorted <- odt[order(odt$Var.Ex, decreasing = TRUE), ]
  # head(sorted)
  n_outl <- round((5 * nrow(sorted)) / 100)
  # cat('5% of', nrow(sorted), 'is', n_outl, '\n')
  sorted$sel[1:n_outl] <- TRUE
  
  out_dt <- data.frame(table(sorted[sorted$sel==TRUE, 'invRui']),
                       ZONE = zone, VTYPE = vtype)
  
  return(out_dt)
})))
head(cl_dt)

table(cl_dt$ZONE, cl_dt$VTYPE)
table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)[1,,,]
table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)['CZA_left',,,]
table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)[1:2,,,]
table(cl_dt$ZONE, cl_dt$sel, cl_dt$VTYPE, cl_dt$invRui)

ggplot(data = cl_dt, aes(x = Var1, y = sqrt(Freq), fill = VTYPE)) +
  facet_wrap(~ZONE) +
  geom_col(position = position_dodge(preserve = "single")) +
  labs(x = '', y = 'square root count', fill = '')
