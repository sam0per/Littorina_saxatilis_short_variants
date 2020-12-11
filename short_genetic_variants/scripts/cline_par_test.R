rm(list = ls())

pkgs <- c("tools", "tidyr", "ggplot2", "data.table", "optparse", "dplyr", "patchwork", "gridExtra", "RColorBrewer")
# Install CRAN packages (if not already installed)
# .inst <- .packages %in% installed.packages()
# if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))

################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-i", "--island"), type = "character", default = NULL,
              help = "Name of the island (CZA, CZB, CZD).", metavar = "character"),
  make_option(c("-v", "--variant-type"), type = "character", default = NULL,
              help = "Names of the type of variants (INDEL:SNP, neu_INDEL:neu_SNP, noneu_INDEL:noneu_SNP).",
              metavar = "character"),
  make_option(opt_str = "--cdf", action = 'store_true', default = FALSE,
              help = "add this flag to plot empirical cumulative distribution."))

opt_parser = OptionParser(usage = paste("Rscript scripts/cline_par_test.R",
                                        "-i CZA",
                                        "-v neu_INDEL:neu_SNP"),
                          option_list=option_list,
                          description = "Test for difference in distribution of cline parameters between variant types.")
opt = parse_args(opt_parser)

if (is.null(opt$island) | is.null(opt$`variant-type`)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# Get cline fits
cl_fl <- list.files(path = "CZCLI006_comp", full.names = TRUE)
cl_fl <- cl_fl[grep(pattern = "NoInv", x = cl_fl)]

cl_ls <- lapply(cl_fl, read.table, header = TRUE)

cl_dt <- as.data.frame(rbindlist(lapply(seq_along(cl_fl), function(x) {
  island <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][2]
  side <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][3]
  zone <- paste(island, side, sep = "_")
  vtype <- strsplit(file_path_sans_ext(basename(cl_fl[[x]])), split = "_")[[1]][4]
  # paste(zone, vtype)
  thresh <- sort(cl_ls[[x]]$Var.Ex, decreasing=TRUE, na.last=TRUE)[round(length(cl_ls[[x]]$cp)*0.05)]
  cl_ls[[x]]$sel = FALSE
  cl_ls[[x]]$sel[cl_ls[[x]]$Var.Ex>=thresh & is.na(cl_ls[[x]]$Var.Ex)==FALSE] = TRUE
  odt <- mutate(cl_ls[[x]], ZONE = zone, VTYPE = vtype)
  return(odt)
})))

cl_dt$VT_sel <- ifelse(test = cl_dt$sel==TRUE, yes = paste0('noneu_', cl_dt$VTYPE), no = paste0('neu_', cl_dt$VTYPE))

vtype_pal <- data.frame(pal = c(brewer.pal(n = 8, name = 'Set2')[c(5, 3)], "#1B9E77", "#666666",
                                brewer.pal(n = 6, name = 'Paired')[c(6, 5)]),
                        vt = c(unique(cl_dt$VT_sel), unique(cl_dt$VTYPE)))

# head(cl_dt)
# data.frame(table(cl_dt$ZONE, cl_dt$VTYPE))
# round(1178*0.05)
# data.frame(table(cl_dt$ZONE, cl_dt$VTYPE, cl_dt$sel))
# data.frame(table(cl_dt$ZONE, cl_dt$VT_sel))

cl_dt <- separate(data = cl_dt, col = ZONE, into = c("ISL", "SIDE"), sep = "_", remove = FALSE)
# cl_dt$av <- round(cl_dt$av)
cl_dt$LGav <- paste(cl_dt$LG, cl_dt$av, sep = "_")
# data.frame(table(cl_dt$LGav))
ois <- cl_dt[cl_dt$sel, ]
# head(ois)
# table(ois$ZONE, ois$VT_sel)

mp_indel_snp <- vector(mode = "list", length = length(unique(ois$ISL)))
names(mp_indel_snp) <- unique(ois$ISL)
for (s in unique(ois$ISL)) {
  ti <- ois[ois$ISL==s, ]
  mp_indel <- unique(as.character(ti$LGav[ti$VT_sel=="noneu_INDEL"]))
  mp_snp <- unique(as.character(ti$LGav[ti$VT_sel=="noneu_SNP"]))
  # print(setdiff(mp_snp, mp_indel))
  # print(setdiff(mp_indel, mp_snp))
  # cat(s, "% unique INDELs", length(setdiff(mp_indel, mp_snp))/length(mp_indel), "\n")
  # cat(s, "% overlap", length(intersect(mp_indel, mp_snp))/length(mp_indel), "\n")
  cat(s, "% unique INDELs", length(setdiff(mp_indel, mp_snp))/length(unique(ti$LGav)), "\n")
  cat(s, "% unique SNPs", length(setdiff(mp_snp, mp_indel))/length(unique(ti$LGav)), "\n")
  cat(s, "% overlap", length(intersect(mp_snp, mp_indel))/length(unique(ti$LGav)), "\n")
  # mp_indel_snp[[s]] <- intersect(mp_snp, mp_indel)
  mp_indel_snp[[s]] <- setdiff(mp_indel, mp_snp)
}
data.frame(table(unlist(mp_indel_snp)))
ois[ois$LGav=="4_12.7", ]
mp_indel_snp <- unique(as.character(unlist(mp_indel_snp)))
ois[ois$LGav==mp_indel_snp[30], ]
mp_ovar <- unique(as.character(ois[ois$LGav %in% mp_indel_snp, "cp"]))
an_fl <- list.files(path = "/Volumes/Seagate Remote Backup/3.indels/annotated", full.names = TRUE)
an_dt <- lapply(an_fl, read.table, header = TRUE)
# lapply(an_dt, head)
mp_an <- lapply(an_dt, function(x) {
  ane <- mutate(x, cp=paste(CHROM, POS, sep = "_"))
  ane <- ane[, (ncol(ane)-3):ncol(ane)]
  ane_mp <- ane[ane$cp %in% mp_ovar, ]
  return(ane_mp)
})
# lapply(mp_an, head)
mp_cp <- as.data.frame(rbindlist(mp_an))
# mp_cp[mp_cp$cp=="Contig884_193337", ]
# ccp <- data.frame(table(mp_cp$cp))
data.frame(table(mp_cp$FUN))
himp <- mp_cp[mp_cp$CAT=="HIGH", ]
unique(himp)
# data.frame(table(ois[as.character(ois$cp) %in% himp$cp, "LGav"]))
ois[ois$cp=="Contig2938_3409", ]
ois[ois$cp=="Contig5500_56612", ]
ois[ois$cp=="Contig2865_93430", ]
ois[ois$cp=="Contig39720_47427", ]
ois[ois$cp=="Contig41173_9942", ]
ois[ois$LGav=="2_22.2", ]
mp_indel_snp == "5_23.9"
mp_indel_snp == "2_22.2"

# isl <- 'CZB'
isl <- opt$island
oisi <- ois[ois$ISL==isl, ]

# length(unique(oisi$LGav))
for (s in unique(oisi$SIDE)) {
  ti <- oisi[oisi$SIDE==s, ]
  mp_indel <- unique(as.character(ti$LGav[ti$VT_sel=="noneu_INDEL"]))
  mp_snp <- unique(as.character(ti$LGav[ti$VT_sel=="noneu_SNP"]))
  # print(setdiff(mp_snp, mp_indel))
  # print(setdiff(mp_indel, mp_snp))
  cat(isl, s, "% unique INDELs", length(setdiff(mp_indel, mp_snp))/length(mp_indel), "\n")
  # cat(isl, s, "% unique SNPs", length(setdiff(mp_snp, mp_indel))/length(unique(ti$LGav)), "\n")
  cat(isl, s, "% overlap", length(intersect(mp_indel, mp_snp))/length(mp_indel), "\n")
}
# oisi[oisi$LGav=="1_48.6", ]
# oisi[oisi$LGav=="5_30.5", ]

outl_c <- aggregate(x = cl_dt$cp, by = list(cl_dt$VT_sel, cl_dt$LGav), length)
outl_v <- split(x = outl_c, f = outl_c$Group.1)

mp_indel <- unique(outl_v$noneu_INDEL$Group.2)
mp_snp <- unique(outl_v$noneu_SNP$Group.2)
outl_indel <- cl_dt[cl_dt$LGav %in% setdiff(mp_indel, mp_snp), ]
length(setdiff(mp_indel, mp_snp))
# outl_indel <- cl_dt[cl_dt$LGav %in% setdiff(mp_snp, mp_indel), ]
# outl_indel <- cl_dt[cl_dt$LGav %in% mp_indel, ]
# outl_indel <- cl_dt[cl_dt$LGav %in% mp_snp, ]
# data.frame(table(outl_indel$LGav))
# table(outl_indel$VT_sel)
oi <- outl_indel[outl_indel$VT_sel=="noneu_INDEL", ]
# oi <- outl_indel[outl_indel$VT_sel=="noneu_SNP", ]
sh <- data.frame(table(as.character(oi$cp)))
sh <- data.frame(table(as.character(oi$LGav)))
# nrow(sh[sh$Freq==max(sh$Freq), ])/nrow(sh)
# nrow(sh[sh$Freq==(max(sh$Freq)-1), ])/nrow(sh)
# nrow(sh[sh$Freq==(max(sh$Freq)-2), ])/nrow(sh)
# nrow(sh[sh$Freq==(max(sh$Freq)-3), ])/nrow(sh)
# nrow(sh[sh$Freq==(max(sh$Freq)-4), ])/nrow(sh)
# nrow(sh[sh$Freq==(max(sh$Freq)-5), ])/nrow(sh)
# oi[oi$cp=="Contig1481_48552",]
# oi[oi$cp=="Contig21581_310",]
# oi[oi$cp=="Contig4084_125665",]

# mp_indel <- unique(outl_v$noneu_INDEL[outl_v$noneu_INDEL$Group.1==isl, 4])
# mp_snp <- unique(outl_v$noneu_SNP[outl_v$noneu_SNP$Group.1==isl, 4])
# outl_indel <- cl_dt[cl_dt$LGav %in% setdiff(mp_indel, mp_snp), ]
# table(outl_indel$LGav)
# table(outl_indel$VT_sel)
# outl_indel[outl_indel$VT_sel=="noneu_SNP", ]

isl_dt <- cl_dt[grepl(pattern = isl, x = cl_dt$ZONE), ]

zs <- split(x = isl_dt, f = isl_dt$ZONE)

# cpars <- c("Centre", "Width", "slope", "p_diff", "Var.Ex")
cpars <- c("Centre", "p_crab", "p_wave", "Var.Ex")

# vts <- "INDEL:SNP"
# vts <- "neu_INDEL:neu_SNP"
# vts <- strsplit(x = vts, split = ":")[[1]]
vts <- strsplit(x = opt$`variant-type`, split = ":")[[1]]
vtype_pal <- as.character(vtype_pal[vtype_pal$vt %in% vts, "pal"])

if (sum(unique(isl_dt$VTYPE)==vts)==2) {
  vcname <- "VTYPE"
} else {
  vcname <- "VT_sel"
}

lapply(X = seq_along(zs), FUN = function(x) {
  
  one_zs <- zs[[x]][zs[[x]][, vcname] %in% vts, ]
  
  d_comp <- vector(mode = "list", length = length(cpars))
  names(d_comp) <- cpars
  
  for (i in cpars) {
    # i <- cpars[3]
    one_p <- one_zs[, c("cp", i, "ZONE", vcname)]
    colnames(one_p) <- c("cp", "est", "ZONE", vcname)
    one_p$par <- i
    # head(one_p)
    d_comp[[i]] <- one_p
    nna_p <- one_p[!is.na(one_p$est), ]
    
    kt <- ks.test(x = unique(nna_p[nna_p[,vcname]==vts[1], "est"]),
                  y = unique(nna_p[nna_p[,vcname]==vts[2], "est"]))
    
    kt_tb <- data.frame(stats = kt$statistic, p.val = kt$p.value, alt = kt$alternative,
                        ZONE = unique(nna_p$ZONE), VTYPE = paste(vts, collapse = ":"), CPAR = i)
    print(kt_tb)
      
    write.table(x = kt_tb, file = paste("results/CPAR_comp/KS", unique(nna_p$ZONE), paste(vts, collapse = "_"),
                                        i, "stats.csv", sep = "_"), append = FALSE, quote = FALSE, sep = ",",
                row.names = FALSE, col.names = TRUE)
  }
  
  dt_comp <- as.data.frame(rbindlist(d_comp))
  # table(dt_comp$ZONE)
  # table(dt_comp$par)
  # table(dt_comp$VT_sel)
  # head(dt_comp)
  
  if (opt$cdf) {
    cpar_h <- ggplot(data = dt_comp, aes_string(x = "est", col = vcname)) +
      facet_wrap(facets = ~par, nrow = 2, scales = 'free') +
      stat_ecdf(geom = "step", pad = FALSE) +
      scale_color_manual(values = vtype_pal) +
      labs(x = '', y = '', col = unique(dt_comp$ZONE)) +
      theme(legend.position = "top",
            legend.title = element_text(size = 12), legend.text = element_text(size = 11),
            axis.text = element_text(size = 11),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 12),
            panel.background = element_blank(),
            strip.background = element_rect(fill = 'white', color = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.line = element_line(size = 0.2, linetype = "solid",
                                     colour = "black"),
            panel.grid = element_line(colour = "gray70", size = 0.2))
  } else {
    cpar_h <- ggplot(data = dt_comp, aes(x = est)) +
      facet_wrap(facets = ~par, nrow = 2, scales = 'free') +
      geom_histogram(aes_string(fill = vcname), bins = 25, col = 'black', position = 'dodge') +
      scale_fill_manual(values = vtype_pal) +
      labs(x = '', y = '', fill = unique(dt_comp$ZONE)) +
      theme(legend.position = c(0.85,0.25),
            legend.title = element_text(size = 12), legend.text = element_text(size = 11),
            axis.text = element_text(size = 11),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 12),
            panel.background = element_blank(),
            strip.background = element_rect(fill = 'white', color = "black"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.line = element_line(size = 0.2, linetype = "solid",
                                     colour = "black"),
            panel.grid = element_line(colour = "gray70", size = 0.2))
  }
  
  # cpar_h
  ggsave(filename = paste('figures/CPAR_comp/HIST', unique(dt_comp$ZONE), paste(vts, collapse = "_"), "cpars.pdf", sep = "_"),
         plot = cpar_h, scale = 0.75, dpi = "screen")
})


