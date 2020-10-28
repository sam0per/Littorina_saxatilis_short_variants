rm(list=ls())

pkgs <- c("tools", "ggplot2", "data.table", "optparse", "dplyr", "patchwork", "gridExtra", "RColorBrewer")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))
################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-i", "--vone"), type = "character", default = NULL,
              help = "csv file with allele count for a type of variant and population (e.g., CRAB INDELs or WAVE coding INDEL).",
              metavar = "character"),
  make_option(opt_str = "--by", type = "character", default = NULL,
              help = "variant name and annotation separated by colon (e.g., INDEL:coding).",
              metavar = "character"),
  make_option(c("-c", "--csv"), type = "character", default = NULL,
              help = "csv file with other info such as ancestral state.", metavar = "character"),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = paste("Rscript scripts/summary_stats.R -i summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv",
                                        "--by VTYPE:INDEL",
                                        "-c results/Lsax_short_var_czs_daf_inv_findv.csv"),
                          option_list=option_list,
                          description = "Generate txt file with haplotypes that can be used as input for the software DH.")
opt = parse_args(opt_parser)

if (is.null(opt$vone) | is.null(opt$by) | is.null(opt$csv)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_WAVE_LEFT_SNP_filt2_59N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_SNP_filt2_66N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_LEFT_SNP_filt2_42N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZD_WAVE_RIGHT_INDEL_filt2_70N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_WAVE_LEFT_INDEL_filt2_59N.csv'))
# head(ic)
# sum(duplicated(ic))
ic <- unique(read.csv(file = opt$vone))

# cnm <- 'INDEL:TEST'
# cnm <- 'INDEL:syn'
# cnm <- 'INDEL:nongenic'
# cnm <- 'SNP:nongenic'
# cnm <- 'ECOT'
# cnm <- 'ANN'
# cnm <- 'ANC'
cnm <- opt$by

# tv <- strsplit(c('A:C'), split = ":")[[1]]
# tv <- strsplit(c('inframe_DEL:inframe_INS'), split = ":")[[1]]
# tv <- strsplit(c('INDEL:SNP'), split = ":")[[1]]
tv <- strsplit(cnm, split = ":")[[1]]

# dt <- unique(read.csv(file = 'results/Lsax_short_var_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_ins_del_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_ins_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_WWSS_czs_daf_inv_findv.csv'))
# head(dt)
dt <- unique(read.csv(file = opt$csv))

if (tv[1] == 'INDEL') {
  
  dt$SIZE <- abs(dt$SIZE)
  dt <- dt[dt$SIZE <= 50, ]
  
}

if (tv[1]=='DEL' | tv[1]=='INS') {
  
  ann_fl <- list.files(path = 'annotated', pattern = paste(levels(ic$ZONE), 'INDEL',
                                                           sep = '_'), full.names = TRUE)
  ann <- read.table(file = ann_fl, header = TRUE)
  
  ann$VTYPE <- 'INDEL'
  ann$CLASS <- tv[1]
  
} else {
  
  ann_fl <- list.files(path = 'annotated', pattern = paste(levels(ic$ZONE), tv[1],
                                                           sep = '_'), full.names = TRUE)
  ann <- read.table(file = ann_fl, header = TRUE)
  
  ann$VTYPE <- tv[1]
  
}

ann$ZONE <- levels(ic$ZONE)
ann$cp <- paste(ann$CHROM, ann$POS, sep = '_')

dtn <- unique(merge(x = dt, y = ic))
dtn <- unique(merge(x = dtn, y = ann))

dtn$VAR <- dtn$VTYPE
dtn$VTYPE <- dtn$CLASS
  
# head(dtn)
# table(dtn$ZONE)
# table(dtn$CLASS)
# table(dtn$VTYPE)
# table(dtn$ECOT)
# table(dtn$CLASS)
# table(dtn$gBGC) 

if (opt$rminv) {
  cat('1.  Varinats inside inversions will be removed.\n')
  dtn <- dtn[dtn$invRui==FALSE, ]
}

# table(dtn$NE_W_Lcomp)
# nrow(unique(dtn))

dtn$DAC <- ifelse(test = dtn$NE_W_Lcomp=='ref_anc:ref_anc', yes = dtn$AC, no = (unique(dtn$N)*2)-dtn$AC)

dtp <- dtn[!is.na(dtn$av) & dtn$DAC!=0 & dtn$DAC!=(unique(dtn$N)*2), ]
# dtp[duplicated(dtp$cp), ]
dtp <- dtp[!duplicated(dtp$cp), ]

# 
# 
# 
# dtp$ANC <- ifelse(test = dtp$NE_W_Lcomp=='ref_anc:ref_anc', yes = as.character(dtp$REF), no = as.character(dtp$ALT))
# if (cnm == 'ANC') {
#   dtp <- dtp[dtp$ANC==tv[1] | dtp$ANC==tv[2], ]
#   # table(dtp$ECOT)
#   dtp$ECOT <- dtp$ANC
# }
# table(dtp$ANC)
# 
# 
# 

## READ TABLE WITH SNPEFF ANNOTATION CLASSES
snpeff <- read.csv(file = 'data/ANN_snpeff_classes.csv')
# str(snpeff)
# head(snpeff)

if (tv[1] != tv[2]) {
  
  if (sum(tv[2] %in% colnames(snpeff))==0) {
    
    effcat <- as.character(snpeff[grepl(pattern = tv[2], x = snpeff$eff), 1])
    
  } else {
    
    effcat <- as.character(snpeff[snpeff[, tv[2]]==TRUE, 1])
    
  }
  
  eff_tar <- as.data.frame(rbindlist(lapply(dtp$ANN, FUN = function(x) {
    stsp <- strsplit(x = as.character(x), split = '\\|')[[1]]
    
    if (sum(grepl(pattern = '&', x = stsp)) > 0) {
      
      tosplit <- stsp[grepl(pattern = '&', x = stsp)]
      stsp <- c(stsp, as.character(unlist(data.frame(strsplit(tosplit, "&")))))
      
    }
    
    if (tv[2] != 'nonsyn') {
      
      ancat <- stsp[stsp %in% as.character(snpeff[,1])]
      tan <- ancat %in% effcat
      # unique(c(effcat, as.character(snpeff[snpeff$nongenic==TRUE, 1])))
      
      if (tv[2] == 'nongenic') {
        
        an1 <- sum(tan) == length(tan)
        
      } else {
        
        if (sum(ancat %in% as.character(snpeff[snpeff$nonsyn==TRUE, 1])) > 0) {
          
          an1 <- FALSE
          
        } else if (sum(ancat %in% as.character(snpeff[snpeff$nongenic==TRUE, 1])) == length(ancat)) {
          
          an1 <- FALSE
          
        } else {
          
          an1 <- TRUE
          
        }
        
      }
      
    } else {
      
      an1 <- sum(stsp %in% effcat) > 0
      
    }
    
    # an1 <- sum(stsp %in% effcat) == 1
    an_dt <- data.frame(an1)
    colnames(an_dt) <- paste(tv[2], unique(dtn$VTYPE), sep = '_')
    return(an_dt)
  })))
  # head(eff_tar)
  # sum(eff_tar[,1])
  eff_tar$VTYPE <- dtp$VTYPE
  
  dtp <- dtp[eff_tar[, 1], ]
  
}

# eff_tar[eff_tar[,1]==FALSE,]

if (tv[2] %in% levels(dt$gBGC)) {
  
  if (nrow(dtp) != nrow(eff_tar)) {
    
    stop("There are annotation codes missing.\n", call.=FALSE)
    
  } else {
    
    dtp <- dtp[dtp$gBGC==tv[2], ]
    
  }
  
}

write.table(x = dtp, file = paste('results/marker_density/MD', unique(dtp$ZONE), levels(ic$ECOT),
                                  colnames(eff_tar)[1], 'count.txt', sep = '_'),
            quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

dac <- dtp[, 'DAC']
df <- data.frame(matrix(data = 0, nrow = unique(dtp$N)*2, ncol = length(dac),
                        dimnames = list(1:(unique(dtp$N)*2), as.character(dtp[, 'cp']))),
                 stringsAsFactors=F)

for (i in 1:ncol(df)) {
  df[, i] <- c(rep(1, dac[i]), rep(0, (unique(dtp$N)*2)-dac[i]))
}
# colSums(df)

fileConn <- file(paste('summary/haplotypes/HAP', levels(ic$ZONE), levels(ic$ECOT),
                       colnames(eff_tar)[1], 'DH.txt', sep = '_'))
writeLines(c('# ', colnames(eff_tar)[1], '\n//\nsegsites: ', length(dac),
             '\npositions: ', paste(round(1:length(dac)/length(dac), 4), collapse = ' '), '\n'),
           con = fileConn, sep = '')
close(fileConn)

write.table(x = df,
            file = paste('summary/haplotypes/HAP', levels(ic$ZONE), levels(ic$ECOT),
                         colnames(eff_tar)[1], 'DH.txt', sep = '_'),
            quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE, append = TRUE)
