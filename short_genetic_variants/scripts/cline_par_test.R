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
              metavar = "character"))

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
  
  odt <- mutate(cl_ls[[x]], ZONE = zone, VTYPE = vtype)
  return(odt)
})))

cl_dt$VT_sel <- ifelse(test = cl_dt$sel==TRUE, yes = paste0('noneu_', cl_dt$VTYPE), no = paste0('neu_', cl_dt$VTYPE))

vtype_pal <- data.frame(pal = c(brewer.pal(n = 8, name = 'Set2')[c(5, 3)], "#1B9E77", "#666666",
                                brewer.pal(n = 6, name = 'Paired')[c(6, 5)]),
                        vt = c(unique(cl_dt$VT_sel), unique(cl_dt$VTYPE)))

# isl <- 'CZA'
isl <- opt$island

isl_dt <- cl_dt[grepl(pattern = isl, x = cl_dt$ZONE), ]

zs <- split(x = isl_dt, f = isl_dt$ZONE)

# cpars <- c("Centre", "Width", "slope", "p_diff", "Var.Ex")
cpars <- c("Centre", "p_diff", "Var.Ex")

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
    # i <- cpars[1]
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
  # cpar_h
  ggsave(filename = paste('figures/CPAR_comp/HIST', unique(dt_comp$ZONE), paste(vts, collapse = "_"), "cpars.pdf", sep = "_"),
         plot = cpar_h, scale = 0.75, dpi = "screen")
})


