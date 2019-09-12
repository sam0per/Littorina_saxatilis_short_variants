rm(list = ls())

.packagesdev = "thomasp85/patchwork"
.packages = c("ggplot2", "dplyr", "reshape2", "parallel", "optparse", "tidyr", "splitstackshape", "data.table", "gdata",
              "adjclust", "matrixStats", "bbmle", "RCurl")
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

option_list = list(
  make_option(c("-V", "--vcf"), type="character", default=NULL,
              help="input VCF genotype table", metavar="character"),
  make_option(c("-D", "--directory"), type="character", default=NULL,
              help="DIR containing files for analysis", metavar="character"),
  make_option(c("-E", "--ecotype"), type="character", default=NULL,
              help="select either crab or wave", metavar="character"),
  make_option(c("-I", "--island"), type="character", default=NULL,
              help="select one or more islands: CZA, CZB and/or CZD", metavar="character"),
  make_option(c("-L", "--linkgrp"), type="integer", default=NULL,
              help="choose a linkage group for testing", metavar="integer"))

opt_parser = OptionParser(option_list=option_list,
                          description = "LD-based approach based on r2 to assign LG to variants in contigs without map position",
                          epilogue = "Example: Rscript scripts/LD_pure_inds.R -V data/genotype.table -D data -E wave -I CZA CZB")
opt = parse_args(opt_parser)

if (is.null(opt$vcf) | is.null(opt$directory) | is.null(opt$ecotype) | is.null(opt$island)){
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

######################
#######  TODO  #######
######################

# filter by minor allels frequency 0.1 after sampling pure Crab or Wave individuals
# plot final results for multiple islands and multiple variants of the same test contigs at different positions
# package snpStats: convert vcf to plink
# https://cran.r-project.org/web/packages/adjclust/vignettes/snpClust.html
# https://bioconductor.org/packages/release/bioc/html/snpStats.html

######################
tabreads = read.table(opt$vcf, header=TRUE, sep = "\t")
spat_dir = opt$directory
ecotype = opt$ecotype
island = strsplit(opt$island, split = " ")[[1]]
testLG = opt$linkgrp
# island = strsplit("CZA CZB", split = " ")[[1]]

#### read spatial data ####
CZ_LCP = lapply(island, function(i) {
  LCP_fl = list.files(path = spat_dir, pattern = paste0(i, "_spatial"), full.names = TRUE)
  LCP = read.csv(LCP_fl)
  LCP_shore = mutate(LCP, shore=i)
  return(LCP_shore)
})
names(CZ_LCP) = island
#### read linkage map ####
linkmap = read.table(list.files(path = spat_dir, pattern = "map", full.names = TRUE), header = TRUE)
linkmap = arrange(linkmap, LG, av)
#### read phenotype (length) data ####
CZ_len = lapply(island, function(i) {
  len_fl = list.files(path = spat_dir, pattern = paste0(i, "_length"), full.names = TRUE)
  len_df = read.csv(len_fl)
  len_shore =  mutate(len_df, shore=i)
  return(len_shore)
})
#### read dissection data ####
CZ_diss = lapply(island, function(i) {
  LCP_fl = list.files(path = spat_dir, pattern = paste0(i, "_dissections"), full.names = TRUE)
  LCP = read.csv(LCP_fl)
  return(LCP)
})
#### identify sex of each snail, using brood pouch and penis data ####
sex <- function(b, p) {
  if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
  if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
  if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N"))==F | b==p) y<-"NA"
  return(y)
}
CZ_sex = lapply(seq_along(island), function(i) {
  mutate(CZ_diss[[i]], sex=apply(CZ_diss[[i]][, c("brood", "penis")],
                                 MARGIN = 1, FUN = function(x) sex(x[1], x[2])),
         shore=island[i])[, c('snail_ID', 'sex', 'shore')]
})
#### merge data ####
CZ_data = lapply(seq_along(island), function(m) {
  joined_df = Reduce(function(x,y) merge(x = x, y = y, by = c("snail_ID", "shore")), list(CZ_LCP[[m]], CZ_len[[m]], CZ_sex[[m]]))
  fin_df = joined_df[joined_df$sex!="NA", ]
  return(fin_df)
})
#### phenotypic cline analysis ####
cline_2c4s <- function(phen,position,sex,cl,cr,lwl,lwr,crab,wave,zs_c,zs_w,sc,shl,sh,sw){
  wl = exp(lwl)
  wr = exp(lwr)
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(position-cl)/wl))  # decreasing
  z_xl <- crab+(wave-crab)*p_xl  # z_xl is expected phenotype for left cline
  z_xl[sex=="female"] <- z_xl[sex=="female"] + zs_c + (zs_w-zs_c)*p_xl[sex=="female"]
  s_xl <- sqrt(sc^2 + 4*p_xl*(1-p_xl)*shl^2 + (p_xl^2)*(sw^2-sc^2))

  # right cline
  p_x <- 1/(1+exp(0-4*(position-cr)/wr))  # increasing
  z_x <- crab+(wave-crab)*p_x  # z_x is expected phenotype for the right cline
  z_x[sex=="female"] <- z_x[sex=="female"] + zs_c + (zs_w-zs_c)*p_x[sex=="female"]
  s_x <- sqrt(sc^2 + 4*p_x*(1-p_x)*sh^2 + (p_x^2)*(sw^2-sc^2))

  # combined cline
  cond <- z_x < z_xl
  z_x[cond] <- z_xl[cond]
  s_x[cond] <- s_xl[cond]
  # z_x[z_x < z_xl] <- z_xl[z_x < z_xl]
  # s_x[z_x < z_xl] <- s_xl[z_x < z_xl]
  minusll <- -sum(dnorm(phen,z_x,s_x,log=T))
  if(crab > wave){minusll <- minusll+1000}
  if(cl > cr){minusll <- minusll+1000}
  # phen_cline = data.frame(phen_cline = z_x, sd_cline = s_x, sex = sex, position = position)
  # return(phen_cline)
  return(minusll)
}
theta.init = list(CZA=list(cl=130, cr=280, lwl=3, lwr=2.3, crab=-2.1, wave=-1.9, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2),
                  CZB=list(cl=70, cr=150, lwl=1.6, lwr=3.9, crab=-2.5, wave=-1.5, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2),
                  CZD=list(cl=80, cr=175, lwl=1.6, lwr=1.6, crab=-2.5, wave=-1.5, zs_c=-0.1, zs_w=-0.1,
                           sc=0.2, shl=0.2, sh=0.2, sw=0.2))
cline_pars = lapply(seq_along(island), function(c) {
  cat("Fitting cline for island", island[c], "...\n")
  mle.cline.2c4s = mle2(cline_2c4s, theta.init[[c]],
                        control=list(parscale=abs(unlist(theta.init[[c]]))),
                        data=list(phen=-log(CZ_data[[c]][CZ_data[[c]]$shore==island[c],]$length_mm),
                                  position=CZ_data[[c]][CZ_data[[c]]$shore==island[c],]$LCmeanDist,
                                  sex=CZ_data[[c]][CZ_data[[c]]$shore==island[c],]$sex))
  cline_est = round(coef(summary(mle.cline.2c4s)), 3)
  return(cline_est)
})
#### sample from crab or wave habitat (commented out) ####
#### dataset(s) is provided on GitHub ####
gitdir = getURL(paste0("https://raw.githubusercontent.com/The-Bioinformatics-Group/Littorina_saxatilis/master/LGid_r2_approach/data/CZA_",
                       ecotype, "50_LCP_ID.csv"))
df_eco = list(read.csv(text = gitdir))
# df_eco = lapply(seq_along(island), function(x) {
#   if (ecotype == "crab") {
#     cl = cline_pars[[x]]["cl", "Estimate"]
#     cr = cline_pars[[x]]["cr", "Estimate"]
#     wl = cline_pars[[x]]["lwl", "Estimate"]
#     wr = cline_pars[[x]]["lwr", "Estimate"]
#     crab = data.frame(LCmeanDist = seq(from = cl+wl, to = cr-wr), snail_ID=NA, stringsAsFactors = FALSE)
#     if (nrow(crab) > 50) {
#       crab = sample_n(crab, size = 50)
#     }
#     for (l in 1:nrow(crab)) {
#       crab[l, "snail_ID"] = as.character(CZ_LCP[[x]][, "snail_ID"])[which(as.integer(CZ_LCP[[x]][, "LCmeanDist"]) == as.integer(crab[l, "LCmeanDist"]))[1]]
#     }
#     crab = crab[complete.cases(crab), ]
#   } else {
#     czsp = arrange(CZ_LCP[[x]], desc(LCmeanDist))
#     wave = data.frame(snail_ID = czsp$snail_ID[c(1:25, (nrow(czsp)-25):nrow(czsp))],
#                       LCmeanDist = czsp$LCmeanDist[c(1:25, (nrow(czsp)-25):nrow(czsp))])
#   }
# })
###############################
# plot samples along transect #
###############################
# df_eco_spa = merge(CZ_data[[1]], df_eco[[1]], by = c("snail_ID"))
# pdf(paste0("../Littorina_saxatilis/LGid_r2_approach/figures/", island, "_", ecotype, nrow(df_eco_spa), "_spatial.pdf"))
# plot(x = CZ_data[[1]]$LCmeanDist, CZ_data[[1]]$length_mm, pch=19, cex=0.5,
#      xlab=paste0(island, " transect position"), ylab = "length (mm)")
# text(x = df_eco_spa$LCmeanDist.x, y = df_eco_spa$length_mm, labels = df_eco_spa$snail_ID, cex = 0.5, pos = 4)
# segments(x0 = df_eco_spa$LCmeanDist.x, y0 = df_eco_spa$length_mm, x1 = df_eco_spa$LCmeanDist.x+5, y1 = df_eco_spa$length_mm)
# dev.off()
###############################
#### find sample names in vcf ####
ecoID = lapply(seq_along(island), function(n) {
  intersect(colnames(tabreads), paste(df_eco[[n]]$snail_ID, "GT", sep = "."))
})
#### subset linkage map by desired linkage group ####
#### this is just for testing the accuracy of the approach ####
#### info about inversions required ####
dtmap = linkmap[linkmap$LG==testLG, ]
invRui = read.xls(list.files(path = spat_dir, pattern = "rearrangements", full.names = TRUE), sheet = 1, header = TRUE)
colnames(invRui) = c("LG", "cluster", "start", "end", "n_SNP", "n_contigs", "median_LD(r2)", "PC1_var")
conlg6 = unique(as.character(dtmap$contig))
conlg6invcf = sample(intersect(conlg6, as.character(tabreads$CHROM)), size = 1)
convcfNNlg6 = sample(intersect(unique(as.character(tabreads$CHROM)), unique(as.character(linkmap$contig))), size = 10)
#### checkpoint: which variants will be used for the test ####
cat("\nThe test variant in LG", testLG, "is:\n")
linkmap[linkmap$contig==conlg6invcf,]
cat("\nThe test variants NOT in LG", testLG, "are:\n")
linkmap[linkmap$contig==convcfNNlg6[1],]
linkmap[linkmap$contig==convcfNNlg6[2],]
linkmap[linkmap$contig==convcfNNlg6[3],]
#### find variants in vcf ####
contID = intersect(as.character(tabreads$CHROM), c(conlg6, convcfNNlg6))
nnGT_tab = rbindlist(lapply(unique(contID), function(c) {
  chr_id = which(tabreads$CHROM==c)
  tabreads[chr_id, !grepl("GT$", x = colnames(tabreads))]
}))
GT_tab = lapply(seq_along(island), function(i) {
  rbindlist(lapply(unique(contID), function(c) {
    chr_id = which(tabreads$CHROM==c)
    tabreads[chr_id, ecoID[[i]]]
  }))
})
#### substitute allele info (REF and ALT) into genotypes (0, 1, 2) ####
GT_gsub = lapply(seq_along(island), function(i) {
  data.frame(lapply(GT_tab[[i]], function(x) {
    gsub("/", "", x)
  }))
})
GT_gsub = lapply(seq_along(island), function(i) {
  data.frame(lapply(GT_gsub[[i]], function(x) {
    gsub("|", "", x, fixed = TRUE)
  }))
})
# uptab = cbind(nnGT_tab, GT_sub)
CZ_GT = lapply(seq_along(island), function(i) {
  cbind(nnGT_tab, GT_gsub[[i]])
})
GT_012 = function(REF, ALT, GT) {
  het1 = paste0(REF, ALT)
  hom0 = paste0(REF, REF)
  hom2 = paste0(ALT, ALT)
  if (GT == het1) {
    return(1)
  } else if (GT == hom0) {
    return(0)
  } else {
    return(2)
  }
}
dat_012_empty = lapply(seq_along(island), function(i) {
  matrix(nrow = nrow(CZ_GT[[i]]), ncol = ncol(GT_gsub[[i]]))
})
GT_col = lapply(seq_along(island), function(i) {
  seq_along(ecoID[[i]]) + 5
})
dat_012 = lapply(seq_along(island), function(isl) {
  cztmp = as.data.frame(CZ_GT[[isl]])
  coltmp = GT_col[[isl]]
  for (i in 1:nrow(cztmp)) {
    for (c in coltmp) {
      ref = cztmp$REF[i]
      alt = cztmp$ALT[i]
      gt = cztmp[i, c]
      dat_012_empty[[isl]][i, c-5] = GT_012(REF = ref, ALT = alt, GT = gt)
    }
  }
  return(dat_012_empty[[isl]])
})
dat_012 = lapply(seq_along(island), function(isl) {
  as.data.frame(cbind(nnGT_tab, dat_012[[isl]]))
})
#### remove invariant sites ####
dat_012_invind = lapply(seq_along(island), function(isl) {
  inv_ind_all = lapply(c(0,1,2), function(n) {
    TF_inv_ind = apply(dat_012[[isl]][, GT_col[[isl]]], MARGIN = 1, FUN = function(x) sum(x==n)!=length(ecoID[[isl]]))
    dat_012[[isl]][TF_inv_ind, ]
  })
  inv_ind_fin = Reduce(function(...) merge(..., by=colnames(dat_012[[isl]])), inv_ind_all)
  return(inv_ind_fin)
})
dat_012_var = lapply(seq_along(island), function(isl) {
  var_T = lapply(c(0,1,2), function(n) {
    var_TF = apply(dat_012_invind[[isl]][, GT_col[[isl]]], MARGIN = 2, FUN = function(x) sum(x==n)!=nrow(dat_012_invind[[isl]]))
    names(var_TF[var_TF==TRUE])
  })
  var_fin_nm = Reduce(intersect, var_T)
  var_fin = dat_012_invind[[isl]][, c(colnames(nnGT_tab), var_fin_nm)]
  return(var_fin)
})
dat_012_info = lapply(seq_along(island), function(isl) {
  info_tmp = dat_012_var[[isl]][, colnames(nnGT_tab)]
  unite(info_tmp, cp, CHROM, POS, sep = "_")
})
#### compute r2 ####
dat_012_cor = lapply(seq_along(island), function(isl) {
  corr_tmp = cor(t(dat_012_var[[isl]][ , -which(names(dat_012_var[[isl]]) %in% colnames(nnGT_tab))]))
  colnames(corr_tmp) = dat_012_info[[isl]]$cp
  rownames(corr_tmp) = dat_012_info[[isl]]$cp
  return(corr_tmp^2)
})
#### find mean and max r2 ####
contID_map = rbindlist(lapply(contID, function(c) {
  map_df = linkmap[which(linkmap$contig==c), ]
  av_rng = range(map_df[, 'av'])
  if (av_rng[2]-av_rng[1] < 3) {
    sample_n(map_df, size = 1)
  }
}))
r2_ukn_map = lapply(seq_along(island), function(isl) {
  rbindlist(lapply(seq_along(contID_map$contig), function(c) {
    unknown = dat_012_cor[[isl]][grepl(pattern = paste0(conlg6invcf, "_"), x = rownames(dat_012_cor[[isl]])), ]
    cor_ukn_kn = as.data.frame(unknown[, grepl(pattern = paste0(contID_map$contig[c], "_"), x = colnames(unknown))])
    r2_stats = data.frame(cont_ukn = rownames(cor_ukn_kn),
                          r2_mean = apply(X = cor_ukn_kn, MARGIN = 1, mean),
                          r2_max = apply(X = cor_ukn_kn, MARGIN = 1, max),
                          cont_map = contID_map$contig[c], av = contID_map$av[c], LG = contID_map$LG[c])
    return(r2_stats)
  }))
})
#### plot r2 against map position ####
#### NOTE: only for one island and one test variant in desired LG ####
dir.create(file.path(getwd(), "figures"))
repodir = getwd()
# repodir = "/Users/samuelperini/Documents/research/projects/Littorina_saxatilis/LGid_r2_approach/"
r2_1ukn_map = r2_ukn_map[[1]][r2_ukn_map[[1]]$cont_ukn==as.character(unique(r2_ukn_map[[1]]$cont_ukn)[1]), ]
r2_1ukn_map = r2_1ukn_map[r2_1ukn_map$cont_map!=conlg6invcf, ]
pdf(paste0(repodir, "figures/", island, "_", ecotype, "_", unique(r2_1ukn_map$cont_ukn), "_mean&max_r2_map.pdf"))
ggplot(data = r2_1ukn_map) +
  facet_wrap(~LG, ncol = length(unique(r2_1ukn_map$LG))) +
  geom_point(aes(x = av, y = r2_mean, col = 'mean r2')) +
  geom_point(aes(x = av, y = r2_max, col = 'max r2')) +
  scale_color_manual(values = c('red', 'black')) +
  labs(col='', x='map position', y=paste('r2 for', unique(r2_1ukn_map$cont_ukn), sep = ' ')) +
  theme(legend.position = 'top',
        strip.text = element_text(face="bold", size=11),
        strip.background = element_rect(fill="lightblue", colour="black",size=1))
dev.off()
cat("\n****Figure is saved as:\n", paste0(repodir, "figures/", island, "_", ecotype, "_",
                                          unique(r2_1ukn_map$cont_ukn), "_mean&max_r2_map.pdf\n"),
    "\nTHE JOB IS DONE!")
