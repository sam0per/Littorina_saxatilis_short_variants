# exclude all individuals that were within 40m (Wave side) or within 15m (Crab side) of the transition

# setwd("/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115/CZCLI006_comp/")
source(file = "/Users/samuelperini/Documents/research/projects/3.indels/Littorina_saxatilis/short_genetic_variants/scripts/colour_text_hadley.R")
.packages = c("optparse", "dplyr", "data.table")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-i", "--inputdir"), type="character", default=NULL,
              help="Path to directory with input data.", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="Path to output directory.", metavar="character"),
  make_option(c("-t", "--transition"), type="character", default="Rock",
              help="Type of habitat transition [default: %default].", metavar="character"),
  make_option(c("-w", "--wave"), type="integer", default=40,
              help="Distance (m) from habitat transition into wave side [default: %default].", metavar="integer"),
  make_option(c("-c", "--crab"), type="integer", default=15,
              help="Distance (m) from habitat transition into crab side [default: %default].", metavar="integer"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$inputdir) | is.null(opt$outputdir)){
  print_help(opt_parser)
  stop("Input and output directories must be declared\n", call.=FALSE)
}

in_files = list.files(path = opt$inputdir, pattern = "spatial_LCP", full.names = TRUE)
in_dat = lapply(1:length(in_files), function(x) read.csv(file = in_files[x]))
# lapply(in_dat, colnames)

hab_tra = read.csv(file = "/Users/samuelperini/Documents/research/projects/3.indels/Anja/Anja_results/20200115/2_Habitat_transitions_new_centre.csv")
hab_tra = hab_tra[hab_tra$Island!="ANG" & hab_tra$Env==opt$transition, ]
hab_tra
# island = unique(hab_tra$Island)
# side = unique(hab_tra$Side)

cat("\n")
pure_dt = invisible(lapply(1:nrow(hab_tra), function(x) {
  PW = hab_tra$Step[x] + opt$wave
  PC = hab_tra$Step[x] - opt$crab
  pr_scr = paste(hab_tra$Island[x], hab_tra$Side[x], "exclude samples from", PC, "m to", PW, "m")
  cat(colourise(pr_scr, "blue"), "\n")
  out_idx = which(in_dat[[x]][, "DistAlongPath"] < PC | in_dat[[x]][, "DistAlongPath"] > PW)
  out_dat = in_dat[[x]][out_idx, ]
  return(out_dat)
}))
# pure_dt
cat(colourise(paste("Check out your output directory", opt$outputdir), "green"), "\n")
invisible(lapply(1:nrow(hab_tra), function(x) {
  write.table(x = pure_dt[[x]][, "snail_ID"], file = paste0(opt$outputdir, hab_tra$Island[x], "_", hab_tra$Side[x],
                                                            "_pure_ecotypes.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.csv(x = pure_dt[[x]], file = paste0(opt$outputdir, hab_tra$Island[x], "_", hab_tra$Side[x],
  #                                           "_pure_ecotypes.csv"), row.names = FALSE)
}))

