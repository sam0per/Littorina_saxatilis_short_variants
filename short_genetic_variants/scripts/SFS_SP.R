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
  make_option(c("-j", "--vtwo"), type = "character", default = NULL,
              help = "csv file with allele count for another type of variant and population (e.g., WAVE INDELs, CRAB SNPs or WAVE non-coding INDELs).",
              metavar = "character"),
  make_option(opt_str = "--by", type = "character", default = NULL,
              help = "column name that defines what to compare (e.g., VTYPE, ECOT, ANN).",
              metavar = "character"),
  make_option(c("-t", "--types"), type = "character", default = NULL,
              help = "two types of comparison (e.g., INDEL:SNP, coding_INDEL:noncoding_INDEL or CRAB:WAVE).",
              metavar = "character"),
  make_option(c("-c", "--csv"), type = "character", default = NULL,
              help = "csv file with other info such as ancestral state.", metavar = "character"),
  make_option(opt_str = "--Test1", action = 'store_true', default = FALSE,
              help = "add this flag to perform Test 1 of Nielsen et al 2005."),
  make_option(opt_str = "--rminv", action = 'store_true', default = FALSE,
              help = "add this flag if variants inside inversions must be removed."))

opt_parser = OptionParser(usage = paste("Rscript scripts/SFS_SP.R -i summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv",
                                        "-j summary/allele_count/AC_CZA_CRAB_SNP_filt2_66N.csv",
                                        "--by ANN",
                                        "-t INDEL:SNP",
                                        "-c results/Lsax_short_var_czs_daf_inv_findv.csv"),
                          option_list=option_list,
                          description = "Generate SFS and test whether there is a significant difference.")
opt = parse_args(opt_parser)

if (is.null(opt$vone) | is.null(opt$vtwo) | is.null(opt$by) | is.null(opt$types) | is.null(opt$csv)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZD_WAVE_RIGHT_INDEL_filt2_70N.csv'))
# sc <- ic
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_LEFT_SNP_filt2_42N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZD_WAVE_RIGHT_INDEL_filt2_70N.csv'))
# ic <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_LEFT_INDEL_filt2_56N.csv'))
# head(ic)
# sum(duplicated(ic))
ic <- unique(read.csv(file = opt$vone))

# sc <- unique(read.table(file = 'annotated/AN_CZA_INDEL.filt2.txt', header = TRUE))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZA_CRAB_SNP_filt2_66N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZA_WAVE_LEFT_INDEL_filt2_59N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZB_WAVE_RIGHT_SNP_filt2_42N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZD_WAVE_RIGHT_SNP_filt2_70N.csv'))
# sc <- unique(read.csv(file = 'summary/allele_count/AC_CZB_CRAB_INDEL_filt2_56N.csv'))
# head(sc)
sc <- unique(read.csv(file = opt$vtwo))

if (unique(ic$N) != unique(sc$N)) {
  stop("The total number of individuals N must be the same.\n", call.=FALSE)
}

# cnm <- 'VTYPE'
# cnm <- 'ECOT'
# cnm <- 'ANN'
# cnm <- 'ANC'
cnm <- opt$by

# tv <- strsplit(c('A:C'), split = ":")[[1]]
# tv <- strsplit(c('frameshift_INS:inframe_INS'), split = ":")[[1]]
# tv <- strsplit(c('DEL:INS'), split = ":")[[1]]
tv <- strsplit(opt$types, split = ":")[[1]]
if (cnm == 'ANN') {
  tv2 <- unlist(strsplit(tv, split = "_"))
}

# dt <- unique(read.csv(file = 'results/Lsax_short_var_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_ins_del_czs_daf_inv_findv.csv'))
# dt <- unique(read.csv(file = 'results/Lsax_short_snp_czs_daf_inv_findv.csv'))
# head(dt)
dt <- unique(read.csv(file = opt$csv))

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZD_WAVE_RIGHT_INDEL_filt2_70N.csv')), split = "_")[[1]]
parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]


if (cnm == 'ANN') {
  
  ann_fl <- list.files(path = 'annotated', pattern = parts[2], full.names = TRUE)
  
  if (tv2[2]==tv2[4]) {
    
      ann <- read.table(ann_fl[grepl(pattern = substr(x = tv2[2], start = 1, stop = 2), x = ann_fl)], header = TRUE)
    
    if (tv2[2]=='DEL' | tv2[2]=='INS') {
      
      ann$VTYPE <- 'INDEL'
      ann$CLASS <- tv2[2]
      
    } else {
      
      ann$VTYPE <- tv2[2]
      
    }
    
    ann$ZONE <- parts[2]
    ann$cp <- paste(ann$CHROM, ann$POS, sep = '_')
    # head(ann)
    # dtn <- merge(x = dt, y = ic)
    
  } else {
    
    ann_ic <- read.table(ann_fl[grepl(pattern = tv2[2], x = ann_fl)], header = TRUE)
    ann_ic$VTYPE <- tv2[2]
    ann_sc <- read.table(ann_fl[grepl(pattern = tv2[4], x = ann_fl)], header = TRUE)
    ann_sc$VTYPE <- tv2[4]
    ann <- rbind(ann_ic, ann_sc)
    ann$ZONE <- parts[2]
    ann$cp <- paste(ann$CHROM, ann$POS, sep = '_')
    # dtn <- merge(x = dt, y = rbind(ic,sc))
    
  }
  dtn <- unique(merge(x = dt, y = rbind(ic,sc)))
  dtn <- unique(merge(x = dtn, y = ann))
  
} else {
  
  dtn <- unique(merge(x = dt, y = rbind(ic,sc)))
  
}

# head(dtn)
# table(dtn$ZONE)
# table(dtn$VTYPE)
# table(dtn$ECOT)
# table(dtn$CLASS)

# tv <- strsplit(c('CRAB:WAVE_LEFT'), split = ":")[[1]]
# tv <- strsplit(c('WAVE_LEFT:WAVE_RIGHT'), split = ":")[[1]]


if (opt$rminv) {
  cat('1.  Varinats inside inversions will be removed.\n')
  dtn <- dtn[dtn$invRui==FALSE, ]
}
# head(dtn)
# table(dtn$AC)
# table(dtn$NE_W_Lcomp)
# nrow(unique(dtn))

dtn$DAC <- ifelse(test = dtn$NE_W_Lcomp=='ref_anc:ref_anc', yes = dtn$AC, no = (unique(dtn$N)*2)-dtn$AC)
# head(dtn[dtn$NE_W_Lcomp=='alt_anc:alt_anc', ])
# head(dtn[dtn$DAC==0, ])
# head(dtn[dtn$DAC==(unique(dtn$N)*2), ])
dtp <- dtn[!is.na(dtn$av) & dtn$DAC!=0 & dtn$DAC!=(unique(dtn$N)*2), ]
# as.character(dtp$REF)
dtp$ANC <- ifelse(test = dtp$NE_W_Lcomp=='ref_anc:ref_anc', yes = as.character(dtp$REF), no = as.character(dtp$ALT))
if (cnm == 'ANC') {
  dtp <- dtp[dtp$ANC==tv[1] | dtp$ANC==tv[2], ]
  # table(dtp$ECOT)
  dtp$ECOT <- dtp$ANC
}
# table(dtp$ANC)
# str(dtp)
# dtp$ANN <- as.character(dtp$ANN)
# table(dtp$DAC)
# hist(dtp$DAC, breaks = 100)
# dtp[dtp$DAF==1,]

# ggplot(data = dtp, aes(x = DAC, fill = VTYPE)) +
#   facet_grid(rows = vars(ZONE), cols = vars(ECOT)) +
#   geom_histogram(col = 'black', position = 'dodge', binwidth = 2)
# 
# ggplot(data = dtp, aes(x = DAC, fill = VTYPE)) +
#   facet_wrap(~VTYPE, scales = 'free_y') +
#   geom_histogram(col = 'black', binwidth = 2)

## READ TABLE WITH SNPEFF ANNOTATION CLASSES
snpeff <- read.csv(file = 'data/ANN_snpeff_classes.csv')
# str(snpeff)
# head(snpeff)
if (exists('tv2')) {
  # colnames(snpeff) %in% tv2
  if (sum(tv2 %in% colnames(snpeff))==0) {
    
    snpeff1 <- as.character(snpeff[grepl(pattern = tv2[1], x = snpeff$eff), 1])
    snpeff2 <- as.character(snpeff[grepl(pattern = tv2[3], x = snpeff$eff), 1])
    
  } else {
    
    snpeff1 <- as.character(snpeff[snpeff[, tv2[1]]==TRUE, 1])
    snpeff2 <- as.character(snpeff[snpeff[, tv2[3]]==TRUE, 1])
    
  }
  
  # if (length(intersect(snpeff1, snpeff2))!=0) {
  #   stop("One class of annotation is shared between the two sets of comparison.\n", call.=FALSE)
  # }
  
  eff_tar <- as.data.frame(rbindlist(lapply(dtp$ANN, FUN = function(x) {
    stsp <- strsplit(x = as.character(x), split = '\\|')[[1]]
    an1 <- sum(stsp %in% snpeff1) > 0
    an2 <- sum(stsp %in% snpeff2) > 0
    an_dt <- data.frame(an1, an2)
    colnames(an_dt) <- tv2[c(1,3)]
    return(an_dt)
  })))
  # head(eff_tar)
  eff_tar$VTYPE <- dtp$VTYPE
  
  # rm_ann <- eff_tar[rowSums(x = eff_tar[, -3])>1, ]
  # write.table(x = dtp[as.integer(row.names(rm_ann)), ], file = paste('results/rm', parts[2],
  #                                                                    paste(unique(dtp$ECOT), collapse = '_'),
  #                                                                    paste(unique(dtp$VTYPE), collapse = '_'),
  #                                                                    'ANN.txt', sep = '_'),
  #             quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  # 
  # ## remove variants that appear in both classes of annotation
  # dtp <- dtp[-as.integer(row.names(rm_ann)), ]
  # eff_tar <- eff_tar[-as.integer(row.names(rm_ann)), ]
  
  if (tv2[1]==tv2[3]) {
    
    dtp$ANN <- ifelse(test = eff_tar[, 1]==TRUE, yes = paste(tv2[1], eff_tar$VTYPE, sep = '_'), no = NA)
    dtp <- dtp[!is.na(dtp$ANN), ]
    
  } else {
    dtp <- cbind(dtp, eff_tar)
    # sum(dtp$inframe)
    # sum(dtp$frameshift)
    dtp <- dtp[dtp[, tv2[1]] == TRUE | dtp[, tv2[3]] == TRUE, ]
    # dtp$ANN <- ifelse(test = eff_tar[, 1]==TRUE, yes = tv[1], no = tv[2])
    dtp$ANN <- ifelse(test = dtp[, tv2[1]]==TRUE, yes = tv[1], no = tv[2])
    
  }
  
}
# dtp$ANN <- gsub(pattern = ' ', replacement = '_', x = dtp$ANN)
# head(dtp)
# table(dtp$ANN)
# sum(eff_tar[eff_tar$VTYPE=='SNP',1])

# Suppose we sample 5 individuals, then there are 10 chromosomes
# The SFS is usually expressed in terms of counts of the derived allele
# since variants must be polymorphic, this means 9 possible classes: 1 - 9 copies of the derived allele

# simple example of SNP data

# snp <- c(1000,200,100,50,40,30,20,10,20) # counts of SNPs in each of the 9 classes

if (length(unique(dtp$ECOT)) > 1) {
  
  dtp <- split(x = dtp, f = dtp$ECOT)
  snp <- merge(x = data.frame(Var1=1:((unique(dtp[[2]]$N)*2)-1)),
               y = data.frame(table(dtp[[2]][dtp[[2]][, cnm]==tv[2], 'DAC'])),
               by = 'Var1', all.x = TRUE)
  indel <- merge(x = snp, y = data.frame(table(dtp[[1]][dtp[[1]][, cnm]==tv[1], 'DAC'])), by = 'Var1', all.x = TRUE)
  
} else {
  
  snp <- merge(x = data.frame(Var1=1:((unique(dtp$N)*2)-1)), y = data.frame(table(dtp[dtp[, cnm]==tv[2], 'DAC'])),
               by = 'Var1', all.x = TRUE)
  # ggplot(data = data.frame(N=1:length(snp), C=snp), aes(x = N, y = C)) +
  #   geom_col(col = 'black')
  
  # and indels
  # indel <- c(180,30,15,10,8,8,8,9,12) # counts of indels in each of the 9 classes
  indel <- merge(x = snp, y = data.frame(table(dtp[dtp[, cnm]==tv[1], 'DAC'])), by = 'Var1', all.x = TRUE)
  
}
indel <- apply(X = indel[,-1], MARGIN = 2, FUN = function(x) {
  ifelse(test = is.na(x), yes = 0, no = x)
})

# indel <- indel[order(indel$Var1), 3] + 1
snp <- indel[, 1] + 1
indel <- indel[, 2] + 1
# snp <- snp[order(snp$Var1), 2]


# snp <- ifelse(test = is.na(snp), yes = 0.00001, no = snp)
# ggplot(data = data.frame(N=1:length(indel), C=indel), aes(x = N, y = C)) +
#   geom_col(col = 'black')

## READ COLOUR PALETTES
ann_pal <- read.csv(file = 'data/ANN_colour_pal.csv')
# str(ann_pal)
# as.character(ann_pal$coding)
ann_pal <- mutate_all(.tbl = ann_pal, .funs = as.character)
row.names(ann_pal) <- ann_pal$VTYPE

if (exists('tv2')) {
  
  if (tv2[1]==tv2[3]) {
    
    fac_pal <- as.character(ann_pal[, tv2[1]])
    fac_pal <- fac_pal[fac_pal!='FALSE']
    
  } else {
    
    fac_pal <- as.character(ann_pal[tv2[2], tv2[c(1,3)]])
    
  }
  
} else {
  
  fac_pal <- data.frame(nb = c('A', 'C', 'G', 'T', 'DEL', 'INS'),
                        pal = c('green', 'blue', 'black', 'red', '#D95F02', '#7570B3'))
  fac_pal <- as.character(fac_pal[fac_pal$nb==tv, 'pal'])
  
}

# class(fac_pal)

ze_pal <- read.csv(file = 'data/ZONE_ECOT_pal.csv')
ze_pal <- mutate_all(.tbl = ze_pal, .funs = as.character)
row.names(ze_pal) <- ze_pal$ZONE
strip_pal <- as.character(ze_pal[parts[2], as.character(unique(dtn$ECOT))])

## ADD ZONE AND ECOTYPE TO FACET LABEL
ann.labs <- unique(c(paste(parts[2], as.character(unique(dtn$ECOT)), gsub(pattern = '_', replacement = ' ', x = tv[1])),
                     paste(parts[2], as.character(unique(dtn$ECOT)), gsub(pattern = '_', replacement = ' ', x = tv[2]))))
names(ann.labs) <- tv

# FIGURES
if (length(strip_pal) > 1) {
  
  EC <- as.character(unique(dtn$ECOT))
  dat <- data.frame(N=rep(1:length(indel), 2), C=c(indel,snp),
                    V=c(rep(tv[1], length(indel)), rep(tv[2], length(indel))),
                    E=c(rep(EC[1], length(indel)), rep(EC[2], length(indel))))
  dat$facet_fill_color <- strip_pal[dat$E]
  names(ann.labs) <- EC
  # head(dat)
  
  DAC_h <- ggplot(data = dat, aes(x = N, y = C-1)) +
    facet_wrap(facets = ~E, nrow = 2, scales = 'free', labeller = labeller(E=ann.labs)) +
    geom_col(aes(fill = V), col = 'black', position = 'dodge') +
    # scale_fill_manual(values = c("#1B9E77", "#666666")) +
    scale_fill_manual(values = fac_pal) +
    labs(x = '', y = 'Derived allele count') +
    theme(legend.position = 'none',
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 10),
          panel.background = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  dummy <- DAC_h
  dummy$layers <- NULL
  dummy <- dummy + geom_rect(data=dat, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf,
                             aes(fill = facet_fill_color)) +
    scale_fill_manual(values = unique(dat$facet_fill_color))
  
  library(gtable)
  
  g1 <- ggplotGrob(DAC_h)
  g2 <- ggplotGrob(dummy)
  
  gtable_select <- function (x, ...) 
  {
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
  }
  
  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip_t", g2$layout$name)
  stript <- grepl(pattern="strip-t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1
  
  new_strips <- gtable_select(g2, panels | strips | stript)
  grid.newpage()
  grid.draw(new_strips)
  
  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }
  ## ideally you'd remove the old strips, for now they're just covered
  new_plot <- gtable_stack(g1, new_strips)
  grid.newpage()
  grid.draw(new_plot)
  
} else {
  
  DAC_h <- ggplot(data = data.frame(N=rep(1:length(indel), 2), C=c(indel,snp), V=c(rep(tv[1], length(indel)),
                                                                                   rep(tv[2], length(indel)))),
                  aes(x = N, y = C-1)) +
    facet_wrap(facets = ~V, nrow = 2, scales = 'free', labeller = labeller(V=ann.labs)) +
    geom_col(aes(fill = V), col = 'black', position = 'dodge') +
    # scale_fill_manual(values = c("#1B9E77", "#666666")) +
    scale_fill_manual(values = fac_pal) +
    labs(x = '', y = 'Derived allele count') +
    theme(legend.position = 'none',
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 10),
          panel.background = element_blank(),
          strip.background = element_rect(fill = strip_pal, color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(size = 0.2, linetype = "solid",
                                   colour = "black"),
          panel.grid = element_line(colour = "gray70", size = 0.2))
  
}

# DAC_h
# ggsave(filename = 'test/figures/dac_hist.pdf', plot = DAC_h)


# Nielsen et al Test 1 compares the likelihood of the indel data, given the frequencies expected from the SNP data
# to the likelihood of the indel data, based on the indel frequencies
# each likelihood is just a multinomial probability

snp_p <- snp/sum(snp)   # ML frequencies of the classes based on SNP data
indel_p <- indel/sum(indel)   # ML frequencies of the classes based on indel data

# 'difference in proportion': (indels in category/all indels) - (SNPs in category/all SNPs).
dDAP_h <- ggplot(data = data.frame(N=1:length(indel), dP=indel_p-snp_p), aes(x = N, y = dP)) +
  geom_point() +
  # labs(x = '', y = paste(tv[1], 'prop. -', tv[2], 'prop.')) +
  labs(x = '', y = 'Difference in proportion') +
  theme(legend.position = 'none',
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2))
# dDAP_h

class_contrib <- rep(0,((unique(dtn$N)*2)-1))
if (opt$Test1) {
  
  LL0 <- sum(indel*log(snp_p)) # log-likelihood given expectations based on SNPs

  LL1 <- sum(indel*log(indel_p)) # log-likelihood given expectations based on indels

  # chi_app <- 2*(LL1-LL0) # chi-square approximation with n_classes-1 = 8 df
  # chi_app
  # can examine the contributions of different classes to the overall chisquare

  # class_contrib <- rep(0,9)
  # i=1
  for (i in 1:((unique(dtn$N)*2)-1)){class_contrib[i] <- 2*((indel[i]*log(indel_p[i])+sum(indel[-i])*log(sum(indel_p[-i])))-((indel[i]*log(snp_p[i]))+sum(indel[-i])*log(sum(snp_p[-i]))))}

  # plot(1:9,class_contrib)
  
} else {
  
  # alternative (better?) test version (testing whether indels and snps differ, rather than whether indels fit snp expectation)
  
  joint_p <- (snp+indel)/sum(snp+indel)
  
  LL0 <- sum(indel*log(joint_p)) + sum(snp*log(joint_p))
  
  LL1 <- sum(indel*log(indel_p)) + sum(snp*log(snp_p))
  
  a <- rep(0,(unique(dtn$N)*2)-1)
  for (i in 1:((unique(dtn$N)*2)-1)){a[i] <- indel[i]*log(indel_p[i]) + snp[i]*log(snp_p[i]) + sum(indel[-i])*log(sum(indel_p[-i])) + sum(snp[-i])*log(sum(snp_p[-i]))}
  b <- rep(0,(unique(dtn$N)*2)-1)
  for (i in 1:((unique(dtn$N)*2)-1)){b[i] <- (indel[i]*log(joint_p[i])) + (snp[i]*log(joint_p[i])) + sum(indel[-i])*log(sum(joint_p[-i])) + sum(snp[-i])*log(sum(joint_p[-i]))}
  
  # class_contrib <- rep(0,((unique(dtn$N)*2)-1))
  for (i in 1:((unique(dtn$N)*2)-1)){class_contrib[i] <- 2*(a[i]-b[i])}
  
}

chi_app <- 2*(LL1-LL0)  # now with 16 df (I think)
cat('Chi-square test statistic =', chi_app, '\n')
# qchisq(p = 0.05, df = 260)
ndf <- ((unique(dtn$N)*2)-1-1)*(length(tv)-1)
cat('Chi-square critical value at 0.01 with', ndf, 'df =', qchisq(p = 0.01, df = ndf), '\n')
cat('Test statistic - critical value =', chi_app-qchisq(p = 0.01, df = ndf), '\n')

chi_tb <- data.frame(chi2_stat=chi_app,
                     chi2_crit=qchisq(p = 0.01, df = ndf),
                     diff=chi_app-qchisq(p = 0.01, df = ndf))


# plot(1:((unique(dtn$N)*2)-1), class_contrib)

mytheme.b <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.7, col = 'black')),
  rowhead = list(fg_params=list(cex = 0.7)))
contr_p <- ggplot(data = data.frame(N=1:((unique(dtn$N)*2)-1), CC=class_contrib)) +
  labs(x = 'Derived allele class', y = 'Contribution') +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(size = 0.2, linetype = "solid",
                                 colour = "black"),
        panel.grid = element_line(colour = "gray70", size = 0.2)) +
  annotation_custom(tableGrob(round(chi_tb, 3), theme = mytheme.b, rows = NULL),
                    xmin=((unique(dtn$N)*2)-41), xmax=((unique(dtn$N)*2)-21),
                    ymin = max(class_contrib)-5, ymax=max(class_contrib)) +
  geom_point(aes(x = N, y = CC))

if (exists('new_plot')) {
  
  # grid.newpage()
  library(cowplot)
  da_cc <- plot_grid(new_plot, dDAP_h, contr_p, align = "v", nrow = 3, rel_heights = c(1/2, 1/4, 1/4))
  # da_cc <- grid.arrange(new_plot, dDAP_h, contr_p)
  # plot_layout(heights = c(2, 1, 1))
  
} else {
  
  da_cc <- DAC_h / dDAP_h / contr_p +
    plot_layout(heights = c(2, 1, 1))
  
}

# da_cc

# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_66N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename('summary/allele_count/AC_CZA_CRAB_INDEL_filt2_59N.csv')), split = "_")[[1]]
# parts <- strsplit(file_path_sans_ext(basename(opt$vone)), split = "_")[[1]]
if (cnm == 'VTYPE') {
  if (parts[3] == 'WAVE') {
    ECOT <- paste(parts[3], parts[4], sep = '_')
    # VTYPE <- parts[5]
  } else {
    ECOT <- parts[3]
    # VTYPE <- parts[4]
  }
  if (opt$rminv) {
    ggsave(filename = paste('figures/SFS', parts[2], ECOT, tv[1], tv[2] ,'chi_noinv.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  } else {
    ggsave(filename = paste('figures/SFS', parts[2], ECOT, tv[1], tv[2] ,'chi.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  }
} else if (cnm == 'ECOT') {
  if (parts[3] == 'WAVE') {
    # ECOT <- paste(parts[3], parts[4], sep = '_')
    VTYPE <- parts[5]
  } else {
    # ECOT <- parts[3]
    VTYPE <- parts[4]
  }
  if (opt$rminv) {
    ggsave(filename = paste('figures/SFS', parts[2], tv[1], tv[2], VTYPE,'chi_noinv.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  } else {
    ggsave(filename = paste('figures/SFS', parts[2], tv[1], tv[2], VTYPE,'chi.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  }
} else {
  if (opt$rminv) {
    ggsave(filename = paste('figures/SFS', parts[2], paste(unique(dtn$ECOT), collapse = '_'), tv[1], tv[2],
                            'chi_noinv.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  } else {
    ggsave(filename = paste('figures/SFS', parts[2], paste(unique(dtn$ECOT), collapse = '_'), tv[1], tv[2],
                            'chi.pdf', sep = "_"), plot = da_cc,
           width = 10, height = 7)
  }
}


# 
# 
# 