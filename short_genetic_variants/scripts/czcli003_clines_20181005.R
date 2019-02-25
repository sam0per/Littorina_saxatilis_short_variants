rm(list=ls())

# 23_2_18 CZ read cline fitting
#
# derived from simci02_czcli02_read_clines_20180125_Roger_no_profiling_no_jittering.R
#
# key modifications:
#
# this version splits transect and fits clines separately for the two contacts
#
# adds spline fit to generate another comparison for var.ex
#
# 25 Jan 2018
#
# NB Jittering and CI disabled while testing other things
#
# Trying - symmetrical cline fit to compare with asymmetrical fit
#        - 'residual variance' option
#        - switch to log(width)
#        - measure 'distance from centroid'
#
# Anja 13 Feb 2018:
# 1. correction in linear_ll
# 2. change of outfile names throughout
# profiling and jittering temporarily switched off
# currently only works for CZA - other zones not included yet
#
# Anja 8 March 2018:
# 1. boundaries for centre estimation for the first cline fit slightly changed, 
# in order to avoid estimates at the boundary that make second fit fail
# 2. 3 jittered replicates
# 3. if centre estimate is at the boundary, Type is now "Stuck" rather than "Cline"
#
# Roger 9 March 2019:
# left and right fits now fully separated

# List of packages for session
.packages = c("bbmle", "rms", "stringr", "tidyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

args = commandArgs(trailingOnly=TRUE)
# install.packages("bbmle")
# library(bbmle)
# install.packages("rms")
# library(rms)


##########################################################################
#
# input shore-specific starting values and boundaries
#
##########################################################################
zone <- args[2]
# zone = "CZA"

if (zone == "ANG"){
  left_transition <- 72 #190 # Bezier path distance of left habitat boundary
  right_transition <- 200 # 280 # Bezier path distance of right habitat boundary
  snail <- read.csv("data/ANG_spatial_LCP_20190121.csv", sep = ",", row.names = TRUE)
  snail <- merge(snail,read.csv("data/ANG_dissections_20150723.csv"),by="snail_ID", all=T)
}

if (zone == "CZA"){
  left_transition <- 200 #190 # Bezier path distance of left habitat boundary
  right_transition <- 277 # 280 # Bezier path distance of right habitat boundary
  snail <- read.csv("data/CZA_spatial_LCP_20190121.csv")
  colnames(snail) = c("idx", "snail_ID", "X", "Y", "Z", "Env", "DistAlongPath")
  snail <- merge(snail, read.csv("data/CZA_dissections_20150506.csv"), by="snail_ID", all=T)
}

if (zone == "CZB"){
  left_transition <- 76 #190 # Bezier path distance of left habitat boundary
  right_transition <- 113 # 280 # Bezier path distance of right habitat boundary
  snail <- read.csv("data/CZB_spatial_LCP_20190121.csv")
  colnames(snail) = c("idx", "snail_ID", "X", "Y", "Z", "Env", "DistAlongPath")
  snail <- merge(snail,read.csv("data/CZB_dissections_20150601.csv"),by="snail_ID", all=T)
}

if (zone == "CZD"){
  left_transition <- 72 #190 # Bezier path distance of left habitat boundary
  right_transition <- 150 # 280 # Bezier path distance of right habitat boundary
  snail <- read.csv("data/CZD_spatial_LCP_20190110.csv")
  colnames(snail) = c("idx", "snail_ID", "X", "Y", "Z", "Env", "DistAlongPath")
  snail <- merge(snail,read.csv("data/CZD_dissections_20150601.csv"),by="snail_ID", all=T)
}



###########################################################################
# read input with line distance and dissection info
###########################################################################

tabreads <- read.table(args[1], header=T, sep = "\t")
# tabreads <- read.table("2016_CZA_capture/tables/all.vcf.depth.table", header=TRUE, sep = "\t")
cz_cols = grepl(pattern = zone, x = colnames(tabreads))
czreads = cbind(tabreads[, 1:5], tabreads[, cz_cols])
# head(czreads[, 1:10])

# reads[is.na(reads)] <- "0,0"

sep <- function(...) {
  dots <- list(...)
  n <- stringr::str_count(dots[[1]][[dots[[2]]]], "\\d+")
  separate_(..., into = sprintf("%s.%d", dots[[2]], 1:n))
}

ad_cols = grepl(pattern = "AD$", x = colnames(czreads))
allreads = czreads %>% Reduce(f = sep, x = colnames(czreads[, ad_cols]))
# head(allreads[,1:11])
# head(tabreads[,1:8])
#reads <- read.table("SIM_CZB_1003.txt", header=T)
#reads <- read.table("2016_CZA_capture/test/CZA_indels_test.allele", sep = "",header = TRUE)
#reads <- read.table("2016_CZA_capture/test/clines/CZA_indels_all_101", sep = "",header = TRUE)

vartype <- args[3]
# vartype = "INDEL"
reads = allreads[allreads$TYPE == vartype, ]
# print(taskid)

# set up distance and start values
path_offset <- left_transition + (right_transition - left_transition)/2 # distance from transect start to new zero point in middle of crab habitat 
snail$position <- snail$DistAlongPath - path_offset 
lt <- left_transition - path_offset
rt <- right_transition - path_offset

# identify sex of each snail, using brood pouch and penis data
sex <- function(b, p) { 
  if(b=="Y" & p=="N" & (is.na(b)==F)) y <- "female"
  if(b=="N" & p=="Y" & (is.na(b)==F)) y <- "male"
  if((b %in% c("Y", "N"))==F | (p %in% c("Y", "N"))==F | b==p) y<-"NA"
  return(y)
}
snail$sex = apply(snail[, c("brood", "penis")], MARGIN = 1, FUN = function(x) sex(x[1], x[2]))



# in what follows, the 'reference' allele is the allele that is more common in the wave
# ecotype. That way, clines are always descending for the left contact, ascending for the right contact
# to make this easier, Bezier path distance is centred in the middle of the crab zone
# therefore, negative values to left, positive to right


# now two separate functions for log-likelihood of cline parameters
# cl=left centre, cr=right centre, wl= left width, wr= right width, [parameters now log(width)]
# crab freq = lpc, left wave freq = lpwl, right wave freq = lpwr
# error (for the locus in question, le, log scale)
# NB frequency parameters on logit scale!!
# given positions (x), read counts (n_w, n_c) 
# this incorporates the probability of reads

cline_left <- function(x,n_w,n_c,cl,lwl,lpcl,lpwl,le){
  pcl <- exp(lpcl)/(1+exp(lpcl)) # get ends back on frequency scale
  pwl <- exp(lpwl)/(1+exp(lpwl))
  
  wl <- exp(lwl) # and width back on natural scale
  
  
  # left cline
  p_xl <- 1-1/(1+exp(0-4*(x-cl)/wl))  # NB '1-' so that it declines
  f_xl <- pcl+(pwl-pcl)*p_xl
  
  
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-f_xl)^2*exp(llAA)+2*f_xl*(1-f_xl)*exp(llAR)+f_xl^2*exp(llRR)  
  Prp[f_xl>1 | f_xl<0] <- 0.0001 # this line should now be superfluous but it does no harm
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  
  return(minusll)
}

cline_right <- function(x,n_w,n_c,cr,lwr,lpcr,lpwr,le){
  pcr <- exp(lpcr)/(1+exp(lpcr)) # get ends back on frequency scale
  pwr <- exp(lpwr)/(1+exp(lpwr))
  
  wr <- exp(lwr) # and width back on natural scale
  
  
  # right cline
  p_xr <- 1/(1+exp(0-4*(x-cr)/wr))  # NB ascending
  f_xr <- pcr+(pwr-pcr)*p_xr  
  
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-f_xr)^2*exp(llAA)+2*f_xr*(1-f_xr)*exp(llAR)+f_xr^2*exp(llRR)  
  Prp[f_xr>1 | f_xr<0] <- 0.0001 # this line should now be superfluous but it does no harm
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  
  return(minusll)
}

# linear effect of distance on frequency, so that we can estimate end frequencies before cline fitting
# here a 'split regression' - one value for pc at position = 0, plus pwl and pwr estimates
# logit input parameters as a way to avoid p outside 0,1
linear_ll <- function(x,n_w,n_c,lpc,lpw,le){
  pc <- exp(lpc)/(1+exp(lpc)) # get ends back on frequency scale
  pw <- exp(lpw)/(1+exp(lpw))
  
  fx <- rep(0.5,length(x))
  fx <- pc+(pw-pc)*abs(x)/max(abs(x)) # straight lines between ends and centre, centre=0, end=1
  
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR)  
  Prp[fx>1 | fx<0] <- 0.0001
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  return(minusll)
}

# and the simplest model - just fixed fx
NC_ll <- function(x,n_w,n_c,lp_all,le){
  p_all <- exp(lp_all)/(1+exp(lp_all)) # get back on frequency scale
  
  fx <- p_all 
  
  # get read probability given frequency
  lf <- lfactorial(n_w+n_c)-(lfactorial(n_w)+lfactorial(n_c))
  llAA <- lf+log(1-exp(le))*n_c+le*n_w
  llAR <- lf+log(0.5)*(n_w+n_c)
  llRR <- lf+le*n_c+log(1-exp(le))*n_w
  Prp <- (1-fx)^2*exp(llAA)+2*fx*(1-fx)*exp(llAR)+fx^2*exp(llRR) #probability of genotype given local frequency, fx
  Prp[fx>1 | fx<0] <- 0.0001
  # sum log(P) over individuals
  minusll <- -sum(log(Prp[Prp>0]))
  return(minusll)
}


# take line i (contig,position,base1,base2,ANGXXX_1_count...,ANGXXX_2_count) and make two columns 
# for read counts and one with line distance, saving contig and position info

#first get snail IDs
snail_ID <- gsub(".AD.1", "", colnames(reads)[seq(6, length(reads[1,]), 2)], fixed=T)
reads = reads[, -c(3:5)]
#start a data frame for output (Type = SL[sex linked], Dup[duplicated], NC[ep_diff<0.1], Linear, Cline[deltaAIC>4 vs linear]), NV[maf<0.05]
rm(out,outci) # clear old output if necessary
clines_out <- args[4]
# clines_out = paste0("2016_CZA_capture/clines/", vartype, "_", zone, ".txt")
out <- data.frame(Index=numeric(),Contig=character(),Position=numeric(),Side=character(),Type=character(),Wave=character(),
                  Centre_left=numeric(),SEcentre_left=numeric(),
                  Centre_right=numeric(),SEcentre_right=numeric(),
                  Width_left=numeric(),SEwidth_left=numeric(),
                  Width_right=numeric(),SEwidth_right=numeric(),
                  p_crab=numeric(),SEp_crab=numeric(),
                  p_wave_left=numeric(),SEp_wave_left=numeric(),
                  p_wave_right=numeric(),SEp_wave_right=numeric(),
                  error_rate=numeric(),Cross_freq=numeric(),Var.Ex=numeric(),Var.Ex.Ex=numeric(),
                  Dev.Ex=numeric(),Dev.Ex.Fit=numeric(),Spline.Ex=numeric(),
                  AIC_NC=numeric(),AIC_linear=numeric(),AIC_S_Cline=numeric(),AIC_cline=numeric())
#write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=T, append=F)
write.table(out, clines_out, quote=F, row.names=F, col.names=T, append=F)
# zone = "CZA" 
# outci <- data.frame(Index=numeric(),Contig=character(),Position=numeric(),Type=character(),Wave=character(),
#                     cl_lo=numeric(),cl_hi=numeric(),
#                     cr_lo=numeric(),cr_hi=numeric(),
#                     wl_lo=numeric(),wl_hi=numeric(),
#                     wr_lo=numeric(),wr_hi=numeric(),
#                     lpc_lo=numeric(),lpc_hi=numeric(),
#                     lpwl_lo=numeric(),lpwl_hi=numeric(),
#                     lpwr_lo=numeric(),lpwr_hi=numeric())
# write.table(outci,paste("CZCLI003_cline_snps_CI_", zone, "_", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=T, append=F)


# start SNP loop, or focus on one or a few SNPs 
for (i in 1:length(reads[,1])) {
  #for (i in 1:500) {  # for first 500
  #for (i in 3473){  # to focus on a specific SNP
  
  
  type <- "NC" # default SNP type is no change
  wave <- "a1" # start by assuming allele 1 is more common in wave than crab
  
  contig <- as.character(reads[i,1])
  SNPpos <- as.numeric(reads[i,2])
  
  a1_count <- as.numeric(reads[i,seq(3,length(reads[1,]),2)]) #number of reads for allele 1
  a2_count <- as.numeric(reads[i,seq(4,length(reads[1,]),2)]) #number of reads for allele 2
  
  wf <- merge(data.frame(snail_ID,a1_count,a2_count),snail,by="snail_ID",all.x=T,all.y=F) #data frame to work with
  wf <- wf[is.na(wf$DistAlongPath)==F,] #remove snails with no position info
  
  
  #start by identifying 'problem' SNPs
  wf$ratio <- wf$a1_count/(wf$a1_count+wf$a2_count)
  wf$reads <- wf$a1_count+wf$a2_count
  
  
  ##### exclude snails with read count <3 for this SNP #####
  wf <- wf[wf$reads>2, ]
  
  #get crude counts of hets and homs, and so rough freq
  h1 <- length(wf$a1_count[wf$a2_count<2])
  h2 <- length(wf$a1_count[wf$a1_count<2])
  het <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1])
  ep <- (2*h1+het)/(2*length(wf$snail_ID))
  f_hets <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1 & wf$sex=='female'])
  m_hets <- length(wf$a1_count[wf$a2_count>1 & wf$a1_count>1 & wf$sex=='male'])
  f_tot <- length(wf$a1_count[wf$sex=='female'])
  m_tot <- length(wf$a1_count[wf$sex=='male'])
  
  #are there more hets than expected under HWE
  exp_tot <- 2*ep*(1-ep)*length(wf$ratio)
  chi_sq_tot <- (het-exp_tot)^2/exp_tot
  chi_sq_tot[het<exp_tot] <- 1 # set to 1 if obs hets < 2pq as expected
  if((f_hets==0) & (m_hets==0)){chi_sq_tot<-0}
  if (pchisq(chi_sq_tot,1,lower.tail = F)<0.01){type <- "Dup>HWE"}
  
  #do hets differ between males and females
  exp_f <- (f_hets+m_hets)*(f_tot/(f_tot+m_tot))
  exp_m <- (f_hets+m_hets)*(m_tot/(f_tot+m_tot))
  chi_sq_m_f <- (f_hets-exp_f)^2/exp_f + (m_hets-exp_m)^2/exp_m
  if((f_hets==0) & (m_hets==0)){chi_sq_m_f<-0}
  if (pchisq(chi_sq_m_f,1,lower.tail = F)<0.01){type <- "SL"}
  
  # is overall minor allele freq too low
  if (ep<0.05 | ep>0.95){type <- "NV"}
  
  # at this point 'good' loci still have type=NC
    # output info for all other types
  if (type!="NC"){
    
    
    out <- data.frame(Index=i,Contig=contig,Position=SNPpos,Side=NA,Type=type,Wave=wave,
                      Centre_left=NA,SEcentre_left=NA,
                      Centre_right=NA,SEcentre_right=NA,
                      Width_left=NA,SEwidth_left=NA,
                      Width_right=NA,SEwidth_right=NA,
                      p_crab=NA,SEp_crab=NA,
                      p_wave_left=NA,SEp_wave_left=NA,
                      p_wave_right=NA,SEp_wave_right=NA,
                      error_rate=NA,Cross_freq=NA,Var.Ex=NA,Var.Ex.Ex=NA,
                      Dev.Ex=NA,Dev.Ex.Fit=NA,Spline.Ex=NA,
                      AIC_NC=NA,AIC_linear=NA,AIC_S_Cline=NA,AIC_cline=NA)
    write.table(out, clines_out, quote=F, row.names=F,
                col.names=F, append=T)
    #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
  }
  
  
  ### loop to do cline fitting etc. 1x with real data, 2x with positions jittered, 'good' loci only
  ############################################################################################################
  
  # for good loci, start with null fit and linear regression, now for left contact and then right contact
  if(type=="NC"){
    
    ## get distribution information needed by rms
    ddist <- datadist(wf[,c("a1_count","a2_count","position")])
    options("datadist" = "ddist")
    
    ## first separate contacts
    wfl <- wf[wf$position<0,] # left side
    wfr <- wf[wf$position>0,] # right side
    
    ## then do left side
    
    ####### left contact #########
    side <- "left"
    
    ### keep original data frame, because wf will be changed
    wf_original <- wfl
    
    for (number1 in seq(1,3)){  ### jitter loop (now 3 runs)
      
      type <- "NC" # reset
      
      lep <- log(ep/(1-ep)) # starting frequency estimate on logit scale
      
      #null
      theta.init <- list(lp_all=lep,le=-5)
      
      mle.NC <- mle2(NC_ll, theta.init, method="L-BFGS-B",
                     upper=list(lp_all=10,le=-1),
                     lower=list(lp_all=-10,le=-10),
                     data=list(x=wfl$position,n_w=wfl$a1_count,n_c=wfl$a2_count))
      AIC.nc <- AIC(mle.NC)
      
      
      #linear
      theta.init <- list(lpc=lep,lpw=lep,le=-5)
      
      mle.linear <- mle2(linear_ll, theta.init, method="L-BFGS-B",
                         upper=list(lpc=10,lpw=10,lle=-1),
                         lower=list(lpc=-10,lpw=-10,le=-10),
                         data=list(x=wfl$position,n_w=wfl$a1_count,n_c=wfl$a2_count))
      AIC.linear<-AIC(mle.linear)
      
      # get fitted end frequencies, logit scale
      pars.lin <- coef(mle.linear)
      se.lin <- summary(mle.linear)@coef[1:3,2]
      lepc <- as.numeric(pars.lin[1])
      lepw <- as.numeric(pars.lin[2])
      
      # and so change in freq (natural scale)
      ep_diff <- exp(lepw)/(1+exp(lepw)) - exp(lepc)/(1+exp(lepc))   
      
      # if necessary, swap alleles to ensure that clines are ascending towards wave, 
      # i.e a1_count is for the allele that is more common in wave
      if (ep_diff<0) {
        temp <- wfl$a1_count
        wfl$a1_count <- wfl$a2_count
        wfl$a2_count <- temp
        lepw <- 0-lepw
        lepc <- 0-lepc
        wave <- "a2"
      }
      
      
      # if difference between ends is large enough, and consistent direction, fit cline models
      if (abs(ep_diff) > 0.1){
        
        # type stays NC unless cline fit is an improvement
        
        # for the cline fit, starting values matter. 
        # here centres are started at the habitat transition
        # widths are started at half distance from transition to mid-point of crab habitat
        # need to output SEs as well because some fits are clearly better than others
        
        
        # using cline.lle
         
        # use estimates to initialize
        # just 1 fit, then repeat with fitted values
        
        theta.init <- list(cl=lt,lwl=log((rt-lt)/2),lpcl=lepc,lpwl=lepw,le=-5)
        mle.cline <- mle2(cline_left, theta.init, method="L-BFGS-B",
                          upper=list(cl=-0.001, # ANJA CHANGED CL HERE FROM 0, 3.3.2018 # clines are between their end of the transect and the centre of the crab habitat
                                     lwl=log(1.5*min(abs(min(wfl$position)-lt),(rt-lt)-1)),
                                     lpcl=10,lpwl=10,le=-1),
                          lower=list(cl=min(wfl$position)+0.001, # ANJA CHANGED CL HERE FROM 0
                                     lwl=-5,
                                     lpcl=-10,lpwl=-10,le=-10),      
                          control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                          data=list(x=wfl$position,n_w=wfl$a1_count,n_c=wfl$a2_count))
        pars <- coef(mle.cline)
        AIC.cline <- AIC(mle.cline)
        
        
        #  re-run using output as starting values (tends to improve fit)
        theta.init <- list(cl=pars[1],lwl=pars[2],lpcl=pars[3],lpwl=pars[4],le=pars[5]) 
        
        mle.cline <- mle2(cline_left, theta.init, method="L-BFGS-B",
                          upper=list(cl=0, # clines are between their end of the transect and the centre of the crab habitat
                                     lwl=log(1.5*min(abs(min(wfl$position)-pars[1]),(rt-lt)-1)),
                                     lpcl=10,lpwl=10,le=-1),
                          lower=list(cl=min(wfl$position),
                                     lwl=-5,
                                     lpcl=-10,lpwl=-10,le=-10),
                          control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                          data=list(x=wfl$position,n_w=wfl$a1_count,n_c=wfl$a2_count))
        
        AIC.cline <- AIC(mle.cline)
        
        

        if(AIC.nc-AIC.cline>4){
          
          type <- "Cline"   #  cline is better than NC
          
          # calculate % var for simple cline
        
          pars <- coef(mle.cline) # replace pars with the estimates from the best fit
          se <- summary(mle.cline)@coef[1:5,2]
          
          if (pars[1] < min(wfl$position)+1) {type="Stuck"}
          
          #to get at se, extract the coef matrix using summary(mle.cline)@coef - estimates are the first column, ses are the second column
        
          # get fitted cline
          fpc <- exp(pars[3])/(1+exp(pars[3])) # get fitted ends back on frequency scale
          fpw <- exp(pars[4])/(1+exp(pars[4]))
          
          fwl <- exp(pars[2]) # and widths back to natural scale
         
        
          # left cline
          fp_xl <- 1-1/(1+exp(0-4*(wfl$position-pars[1])/fwl))  # NB '1-' so that it declines
          ff_xl <- fpc+(fpw-fpc)*fp_xl
          
          # estimate goodness of cline fit (rough % var explained)
          var.res <- var((wfl$a1_count/(wfl$a1_count+wfl$a2_count))-ff_xl)
          var.tot <- var(wfl$a1_count/(wfl$a1_count+wfl$a2_count))
          R2 <- (var.tot-var.res)*100/var.tot
          fitglm <- glm(cbind(wfl$a1_count,wfl$a2_count)~ff_xl,family=binomial(link="logit"))
          res_dev <- fitglm$deviance
          dev_ex <- 100*(fitglm$null.deviance-res_dev)/fitglm$null.deviance
          
          # how much variance would we expect this cline to explain?
          # need a few iterations to get a mean
          R2_fit <- 0
          res_dev_fit <- 0
          dev_fit <- 0
          for (cc in 1:100){
            # first sample genotypes
            fit_gen <- rbinom(length(ff_xl),2,ff_xl)/2
            # allow for error
            fit_gen[fit_gen==0] <- fit_gen[fit_gen==0] + exp(pars[5])
            fit_gen[fit_gen==1] <- fit_gen[fit_gen==1] - exp(pars[5])
            # then read counts
            fit_a1 <- rbinom(length(ff_xl),(wfl$a1_count+wfl$a2_count),fit_gen)
            # then get var.ex
            var.res <- var((fit_a1/(wfl$a1_count+wfl$a2_count))-ff_xl)
            var.tot <- var(fit_a1/(wfl$a1_count+wfl$a2_count))
            R2_fit <- R2_fit+(var.tot-var.res)*100/var.tot
            res_dev_fit <- res_dev_fit + glm(cbind(fit_a1,(wfl$a1_count+wfl$a2_count-fit_a1))~ff_xl,family=binomial(link="logit"))$deviance
            dev_fit <- dev_fit + glm(cbind(fit_a1,(wfl$a1_count+wfl$a2_count-fit_a1))~ff_xl,family=binomial(link="logit"))$null.deviance
          } #end of var.ex.ex loop
          
          R2_fit <- R2_fit/cc
          res_dev_fit <- res_dev_fit/cc
          dev_fit <- dev_fit/cc
          dev_ex_fit <- 100*(dev_fit-res_dev_fit)/dev_fit
          
          # then do spline fit - 6 knots to give 5 df like the cline fit
          lspline <- glm(formula=cbind(a1_count,a2_count)~rcs(position,6), family=binomial(link="logit"), data=wfl)
          spline.varex <- 100*(lspline$null.deviance-lspline$deviance)/lspline$null.deviance
        
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type=type,Wave=wave,
                            Centre_left=as.numeric(pars[1]),SEcentre_left=as.numeric(se[1]),
                            Centre_right=NA,SEcentre_right=NA,
                            Width_left=as.numeric(pars[2]),SEwidth_left=as.numeric(se[2]),
                            Width_right=NA,SEwidth_right=NA,
                            p_crab=as.numeric(pars[3]),SEp_crab=as.numeric(se[3]),
                            p_wave_left=as.numeric(pars[4]),SEp_wave_left=as.numeric(se[4]),
                            p_wave_right=NA,SEp_wave_right=NA,
                            error_rate=as.numeric(pars[5]),Cross_freq=NA,Var.Ex=R2,Var.Ex.Ex=R2_fit,
                            Dev.Ex=dev_ex,Dev.Ex.Fit=dev_ex_fit,Spline.Ex=spline.varex,
                            AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=AIC.cline)
          #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
          write.table(out, clines_out, quote=F, row.names=F,
                      col.names=F, append=T)
          
        
          # if (number1 == 1){ # changed from 1 to 2 to disable during trials
          #   ci <- NA
          #   try(ci <- confint(mle.cline,c(1:7)))
          #   
          #   
          #   if (typeof(ci)=="double"){
          #     outci <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Type=type,Wave=wave,
          #                       cl_lo=ci[1,1],cl_hi=ci[1,2],
          #                       cr_lo=ci[2,1],cr_hi=ci[2,2],
          #                       wl_lo=ci[3,1],wl_hi=ci[3,2],
          #                       wr_lo=ci[4,1],wr_hi=ci[4,2],
          #                       lpc_lo=ci[5,1],lpc_hi=ci[5,2],
          #                       lpwl_lo=ci[6,1],lpwl_hi=ci[6,2],
          #                       lpwr_lo=ci[7,1],lpwr_hi=ci[7,2])
          #     write.table(outci,paste("CZCLI003_cline_snps_CI_", zone, "_", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append=T)
          #   }
          # }
        } # end of cline output loop
        
        if(type=="NC"){
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type="NC>0.1",Wave=wave,
                            Centre_left=NA,SEcentre_left=NA,
                            Centre_right=NA,SEcentre_right=NA,
                            Width_left=NA,SEwidth_left=NA,
                            Width_right=NA,SEwidth_right=NA,
                            p_crab=as.numeric(pars.lin[1]),SEp_crab=as.numeric(se.lin[1]),
                            p_wave_left=as.numeric(pars.lin[2]),SEp_wave_left=as.numeric(se.lin[2]),
                            p_wave_right=as.numeric(pars.lin[3]),SEp_wave_right=as.numeric(se.lin[3]),
                            error_rate=as.numeric(pars.lin[4]),Cross_freq=NA,Var.Ex=NA,Var.Ex.Ex=NA,
                            Dev.Ex=NA,Dev.Ex.Fit=NA,Spline.Ex=NA,
                            AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=AIC.cline)
          write.table(out, clines_out, quote=F, row.names=F,
                      col.names=F, append=T)
          #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
        } # end of NC>0.1 loop
        
      } #end of "if diff>0.1"
      
      else {
        if(type=="NC"){
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type="NC",Wave=wave,
                            Centre_left=NA,SEcentre_left=NA,
                            Centre_right=NA,SEcentre_right=NA,
                            Width_left=NA,SEwidth_left=NA,
                            Width_right=NA,SEwidth_right=NA,
                            p_crab=as.numeric(pars.lin[1]),SEp_crab=as.numeric(se.lin[1]),
                            p_wave_left=as.numeric(pars.lin[2]),SEp_wave_left=as.numeric(se.lin[2]),
                            p_wave_right=as.numeric(pars.lin[3]),SEp_wave_right=as.numeric(se.lin[3]),
                            error_rate=as.numeric(pars.lin[4]),Cross_freq=NA,Var.Ex=NA,Var.Ex.Ex=NA,
                            Dev.Ex=NA,Dev.Ex.Fit=NA,Spline.Ex=NA,
                            AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=NA)
          write.table(out, clines_out, quote=F, row.names=F,
                      col.names=F, append=T)
          #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
        }
      } # end of NC output
      
      
      
      ### jittering of positions
      # first set to orginal data frame again
      wfl = wf_original
      # generate jittered distances
      wfl$position = rnorm(length(wf_original$position), mean = wf_original$position, sd = 0.5)
      
    } # end of "repeat 3x with jittered positions" loop
    
        # end of left side analysis
        
        ####### right contact #########
        
    side <- "right"
    
    ### keep original data frame, because wf will be changed
    wf_original <- wfr
    
    for (number1 in seq(1,3)){  ### jitter loop (now 3 runs)
    
    type <- "NC" # reset
        
        
        
        #null
        theta.init <- list(lp_all=lep,le=-5)
        
        mle.NC <- mle2(NC_ll, theta.init, method="L-BFGS-B",
                       upper=list(lp_all=10,le=-1),
                       lower=list(lp_all=-10,le=-10),
                       data=list(x=wfr$position,n_w=wfr$a1_count,n_c=wfr$a2_count))
        AIC.nc <- AIC(mle.NC)
        
        
        #linear
        theta.init <- list(lpc=lep,lpw=lep,le=-5)
        
        mle.linear <- mle2(linear_ll, theta.init, method="L-BFGS-B",
                           upper=list(lpc=10,lpw=10,lle=-1),
                           lower=list(lpc=-10,lpw=-10,le=-10),
                           data=list(x=wfr$position,n_w=wfr$a1_count,n_c=wfr$a2_count))
        AIC.linear<-AIC(mle.linear)
        
        # get fitted end frequencies, logit scale
        pars.lin <- coef(mle.linear)
        se.lin <- summary(mle.linear)@coef[1:3,2]
        lepc <- as.numeric(pars.lin[1])
        lepw <- as.numeric(pars.lin[2])
        
        # and so change in freq (natural scale)
        ep_diff <- exp(lepw)/(1+exp(lepw)) - exp(lepc)/(1+exp(lepc))   
        
        # if change is negative, swap alleles to ensure that clines are ascending towards wave, 
        # i.e a1_count is for the allele that is more common in wave
        if (ep_diff<0) {
          temp <- wfr$a1_count
          wfr$a1_count <- wfr$a2_count
          wfr$a2_count <- temp
          lepw <- 0-lepw
          lepc <- 0-lepc
          wave <- "a2"
        }
        
        
        # if difference between ends is large enough, and consistent direction, fit cline models
        if (abs(ep_diff) > 0.1){
          
          # type stays NC unless cline fit is an improvement
          
          # for the cline fit, starting values matter. 
          # here centres are started at the habitat transition
          # widths are started at half distance from transition to mid-point of crab habitat
          # need to output SEs as well because some fits are clearly better than others
          
          
          # using cline.lle
          
          # use estimates to initialize
          # just 1 fit, then repeat with fitted values
          
          theta.init <- list(cr=lt,lwr=log((rt-lt)/2),lpcr=lepc,lpwr=lepw,le=-5)
          mle.cline <- mle2(cline_right, theta.init, method="L-BFGS-B",
                            upper=list(cr=max(wfr$position)-0.001, # ANJA CHANGED CL HERE FROM 0 # clines are between their end of the transect and the centre of the crab habitat
                                       lwr=log(1.5*min(abs(max(wfr$position)-rt),(rt-lt)-1)),
                                       lpcr=10,lpwr=10,le=-1),
                            lower=list(cr=0.001, # ANJA CHANGED CL HERE FROM 0
                                       lwr=-5,
                                       lpcr=-10,lpwr=-10,le=-10),      
                            control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                            data=list(x=wfr$position,n_w=wfr$a1_count,n_c=wfr$a2_count))
          pars <- coef(mle.cline)
          AIC.cline <- AIC(mle.cline)
          
          
          #  re-run using output as starting values (tends to improve fit)
          theta.init <- list(cr=pars[1],lwr=pars[2],lpcr=pars[3],lpwr=pars[4],le=pars[5]) 
          
          mle.cline <- mle2(cline_right, theta.init, method="L-BFGS-B",
                            upper=list(cr=max(wfr$position), # clines are between their end of the transect and the centre of the crab habitat
                                       lwr=log(1.5*min(abs(max(wfr$position)-pars[1]),(rt-lt)-1)),
                                       lpcr=10,lpwr=10,le=-1),
                            lower=list(cr=0,
                                       lwr=-5,
                                       lpcr=-10,lpwr=-10,le=-10),
                            control=list(parscale=abs(unlist(theta.init))),            # parscale means that all step sizes are relative to starting values
                            data=list(x=wfr$position,n_w=wfr$a1_count,n_c=wfr$a2_count))
          
          AIC.cline <- AIC(mle.cline)
          
          
          
          if(AIC.nc-AIC.cline>4){
            
            type <- "Cline"   #  cline is better than NC
            
            
            # calculate % var for simple cline
            
            pars <- coef(mle.cline) # replace pars with the estimates from the best fit
            se <- summary(mle.cline)@coef[1:5,2]
            
            if (pars[1] > max(wfr$position)-1) {type="Stuck"}
            
            #to get at se, extract the coef matrix using summary(mle.cline)@coef - estimates are the first column, ses are the second column
            
            # get fitted cline
            fpc <- exp(pars[3])/(1+exp(pars[3])) # get fitted ends back on frequency scale
            fpw <- exp(pars[4])/(1+exp(pars[4]))
            
            fwr <- exp(pars[2]) # and widths back to natural scale
            
            
            # left cline
            #fp_xl <- 1-1/(1+exp(0-4*(wfr$position-pars[1])/fwl))  # NB '1-' so that it declines
            #ff_xl <- fpc+(fpw-fpc)*fp_xl
            
            # right cline
            fp_xr <- 1/(1+exp(0-4*(wfr$position-pars[1])/fwr))  # NB ascending
            ff_xr <- fpc+(fpw-fpc)*fp_xr  
            
            # estimate goodness of cline fit (rough % var explained)
            var.res <- var((wfr$a1_count/(wfr$a1_count+wfr$a2_count))-ff_xr)
            var.tot <- var(wfr$a1_count/(wfr$a1_count+wfr$a2_count))
            R2 <- (var.tot-var.res)*100/var.tot
            fitglm <- glm(cbind(wfr$a1_count,wfr$a2_count)~ff_xr,family=binomial(link="logit"))
            res_dev <- fitglm$deviance
            dev_ex <- 100*(fitglm$null.deviance-res_dev)/fitglm$null.deviance
            
            # how much variance would we expect this cline to explain?
            # need a few iterations to get a mean
            R2_fit <- 0
            res_dev_fit <- 0
            dev_fit <- 0
            for (cc in 1:100){
              # first sample genotypes
              fit_gen <- rbinom(length(ff_xr),2,ff_xr)/2
              # allow for error
              fit_gen[fit_gen==0] <- fit_gen[fit_gen==0] + exp(pars[5])
              fit_gen[fit_gen==1] <- fit_gen[fit_gen==1] - exp(pars[5])
              # then read counts
              fit_a1 <- rbinom(length(ff_xr),(wfr$a1_count+wfr$a2_count),fit_gen)
              # then get var.ex
              var.res <- var((fit_a1/(wfr$a1_count+wfr$a2_count))-ff_xr)
              var.tot <- var(fit_a1/(wfr$a1_count+wfr$a2_count))
              R2_fit <- R2_fit+(var.tot-var.res)*100/var.tot
              res_dev_fit <- res_dev_fit + glm(cbind(fit_a1,(wfr$a1_count+wfr$a2_count-fit_a1))~ff_xr,family=binomial(link="logit"))$deviance
              dev_fit <- dev_fit + glm(cbind(fit_a1,(wfr$a1_count+wfr$a2_count-fit_a1))~ff_xr,family=binomial(link="logit"))$null.deviance
            } # end of var.ex.ex loop
            
            R2_fit <- R2_fit/cc
            res_dev_fit <- res_dev_fit/cc
            dev_fit <- dev_fit/cc
            dev_ex_fit <- 100*(dev_fit-res_dev_fit)/dev_fit
            
            # then do spline fit - 6 knots to give 5 df like the cline fit
            rspline <- glm(formula=cbind(a1_count,a2_count)~rcs(position,6), family=binomial(link="logit"), data=wfr)
            spline.varex <- 100*(rspline$null.deviance-rspline$deviance)/rspline$null.deviance
            
            out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type=type,Wave=wave,
                              Centre_left=NA,SEcentre_left=NA,
                              Centre_right=as.numeric(pars[1]),SEcentre_right=as.numeric(se[1]),
                              Width_left=NA,SEwidth_left=NA,
                              Width_right=as.numeric(pars[2]),SEwidth_right=as.numeric(se[2]),
                              p_crab=as.numeric(pars[3]),SEp_crab=as.numeric(se[3]),
                              p_wave_left=NA,SEp_wave_left=NA,
                              p_wave_right=as.numeric(pars[4]),SEp_wave_right=as.numeric(se[4]),
                              error_rate=as.numeric(pars[5]),Cross_freq=NA,Var.Ex=R2,Var.Ex.Ex=R2_fit,
                              Dev.Ex=dev_ex,Dev.Ex.Fit=dev_ex_fit,Spline.Ex=spline.varex,
                              AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=AIC.cline)
            write.table(out, clines_out, quote=F, row.names=F,
                        col.names=F, append=T)
            #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
            
            
            # if (number1 == 1){ # changed from 1 to 2 to disable during trials
            #   ci <- NA
            #   try(ci <- confint(mle.cline,c(1:7)))
            #   
            #   
            #   if (typeof(ci)=="double"){
            #     outci <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Type=type,Wave=wave,
            #                       cl_lo=ci[1,1],cl_hi=ci[1,2],
            #                       cr_lo=ci[2,1],cr_hi=ci[2,2],
            #                       wl_lo=ci[3,1],wl_hi=ci[3,2],
            #                       wr_lo=ci[4,1],wr_hi=ci[4,2],
            #                       lpc_lo=ci[5,1],lpc_hi=ci[5,2],
            #                       lpwl_lo=ci[6,1],lpwl_hi=ci[6,2],
            #                       lpwr_lo=ci[7,1],lpwr_hi=ci[7,2])
            #     write.table(outci,paste("CZCLI003_cline_snps_CI_", zone, "_", taskid, ".txt", sep=""), quote=F, row.names=F, col.names=F, append=T)
            #   }
            # }
          
        } #end of cline output
        
        if(type=="NC"){
        out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type="NC>0.1",Wave=wave,
                          Centre_left=NA,SEcentre_left=NA,
                          Centre_right=NA,SEcentre_right=NA,
                          Width_left=NA,SEwidth_left=NA,
                          Width_right=NA,SEwidth_right=NA,
                          p_crab=as.numeric(pars.lin[1]),SEp_crab=as.numeric(se.lin[1]),
                          p_wave_left=as.numeric(pars.lin[2]),SEp_wave_left=as.numeric(se.lin[2]),
                          p_wave_right=as.numeric(pars.lin[3]),SEp_wave_right=as.numeric(se.lin[3]),
                          error_rate=as.numeric(pars.lin[4]),Cross_freq=NA,Var.Ex=NA,Var.Ex.Ex=NA,
                          Dev.Ex=NA,Dev.Ex.Fit=NA,Spline.Ex=NA,
                          AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=AIC.cline)
        write.table(out, clines_out, quote=F, row.names=F,
                    col.names=F, append=T)
        #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
        } # end of NC>0.1 output
        
        
      } #end of "if diff>0.1"
      
      else {
        if(type=="NC"){
          out <- data.frame(Index=paste(i,number1,sep="_"),Contig=contig,Position=SNPpos,Side=side,Type="NC",Wave=wave,
                            Centre_left=NA,SEcentre_left=NA,
                            Centre_right=NA,SEcentre_right=NA,
                            Width_left=NA,SEwidth_left=NA,
                            Width_right=NA,SEwidth_right=NA,
                            p_crab=as.numeric(pars.lin[1]),SEp_crab=as.numeric(se.lin[1]),
                            p_wave_left=as.numeric(pars.lin[2]),SEp_wave_left=as.numeric(se.lin[2]),
                            p_wave_right=as.numeric(pars.lin[3]),SEp_wave_right=as.numeric(se.lin[3]),
                            error_rate=as.numeric(pars.lin[4]),Cross_freq=NA,Var.Ex=NA,Var.Ex.Ex=NA,
                            Dev.Ex=NA,Dev.Ex.Fit=NA,Spline.Ex=NA,
                            AIC_NC=AIC.nc,AIC_linear=AIC.linear,AIC_S_Cline=NA,AIC_cline=NA)
          write.table(out, clines_out, quote=F, row.names=F,
                      col.names=F, append=T)
          #write.table(out,"2016_CZA_capture/test/clines/CZA_cline_indels_test.txt", quote=F, row.names=F, col.names=F, append=T)
        }
      } # end of NC output
      
      
      
      
      ### jittering of positions
      # first set to orginal data frame again
      wfr = wf_original
      # generate jittered distances
      wfr$position = rnorm(length(wf_original$position), mean = wf_original$position, sd = 0.5)
      
    } # end of "repeat 10x with jittered positions" loop
  
  # end of right side analysis
    
  } # end of if type=NC
  
  print(c(as.character(i),ep,type)) # print just to monitor progress
  
} # end of snp loop

