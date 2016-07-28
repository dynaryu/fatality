# Modified from Stan-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
#graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#require(rstan)
#source("DBDA2E-utilities.R")
#library(MASS)
#require(ggplot2)
#library(scales)

# copied from DBDA2E-utilities.R
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

# quantile function for truncated distribution
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

# random generaten for truncated distribution
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}

# read USGS exposure catalogue
read_expo_cat <- function(fileName = '/Users/hyeuk/Downloads/EXPO_CAT_2007_12.csv') {

  # estimate using the PAGER database
  expo_cat <- read.csv(fileName)

  # names(expo_cat)
  # [11] "ISO_code"
  # [17] "PAGER_prefShakingDeaths"

  # filter only ID and remove NaN value
  sel_dat <- expo_cat[expo_cat$ISO_code == 'ID' & complete.cases(expo_cat[,17]),]
  #nrow(sel_dat) # 71

  # remove two events to avoid any duplication
  sel_dat <- sel_dat[ (sel_dat$eqID != 200503281609) & (sel_dat$eqID != 200605262253), ]

  # observed fatality 
  observed <- sel_dat$PAGER_prefShakingDeaths
  eqID <- sel_dat$eqID

  mag_observed <- floor(log10(observed)) + 2
  mag_observed[is.infinite(mag_observed)] = 1 # 0-1 range

  # exposed population including rural and urban
  # U090: number of people exposed to MMI 9.0 +/- 0.25, or MMI 8.75 to MMI 9.25
  mmi_list <- seq(4.0, 10.0, 0.5)
  pop <- matrix(NA,ncol=length(mmi_list),nrow=length(observed))

  for (val in mmi_list) {
    str_R <- paste('R',sprintf("%03d", val*10),sep="")
    str_U <- paste('U',sprintf("%03d", val*10),sep="")
    icol <- which(mmi_list==val)
    pop[,icol] <- sel_dat[,str_R] + sel_dat[,str_U]
  }

  expo_cat <- data.frame(eqID=eqID, obs=observed, mag_obs=mag_observed)

  return(list(pop, expo_cat, mmi_list))
} 


compute_HDI_prob <- function(fat_by_event){

  # HDI of estimated fatality
  nevents <- nrow(fat_by_event)
  nfat <- ncol(fat_by_event)

  max_by_event <- apply(fat_by_event, 1, max)

  order_mag <- ceiling(log10(max_by_event))

  fatHDI <- matrix( NA, nrow = nevents, ncol = 2)
  prob_mag <- matrix(NA, nrow = nevents, ncol = max(order_mag) + 1)

  for (i in 1:nevents) {

    temp = fat_by_event[i, ]
    fatHDI[i,] = HDIofMCMC(temp)

    cprob <- numeric(length = order_mag[i] + 1)
    for (j in 1:order_mag[i]) {
        cprob[j] <- sum(temp < 10^(j-1)) / nfat
    } 
    cprob[order_mag[i]+1] <- 1.0

    prob_mag[i,1:(order_mag[i]+1)] <- c(cprob[1], diff(cprob))

  }

  return(list(fatHDI, prob_mag))
}

compute_fat_corr <- function( mcmcMat , pop, mmi_list, nsamples, saveName ) {

  # posterior prediction
  #mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = nrow( mcmcMat )
  chain_a = mcmcMat[,"a"]
  chain_b = mcmcMat[,"b"]
  chain_s = mcmcMat[,"s"]

  # Posterior prediction:
  nevents = nrow(pop)
  nmmi = ncol(pop)
  nfat = nsamples * chainLength

  fat_rate <- matrix( NA, nrow = nmmi, ncol = nfat)

  for ( i in 1:chainLength) {
      mu0 <- exp(chain_b[i]*mmi_list + chain_a[i])
      rate_ = chain_s[i]/mu0
      shape_ = chain_s[i] + 0*mu0

      for (j in 1:nsamples) {
        fat_rate[, nsamples*(i-1)+j]  = rtrunc(1,"gamma",b=1,rate=rate_, shape=shape_)   
      } 
  }

  fat_by_event = pop %*% fat_rate # (nevents, nmmi) * (nmmi, nfat) = (nevents, nfat)

  # HDI of estimated fatality
  dummy = compute_HDI_prob(fat_by_event)

  if ( !is.null(saveName) ) {
    saveRDS(fat_by_event, paste('fat_by_event_corr_', saveName,'.RDS',sep=""))
    saveRDS(dummy[[1]], paste('fatHDI_corr_', saveName,'.RDS',sep=""))
    saveRDS(dummy[[2]], paste('prob_mag_corr_', saveName,'.RDS',sep=""))
  }
}


compute_fat_no_corr <- function( mcmcMat , pop, mmi_list, nsamples, saveName ) {

  # posterior prediction
  #mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = nrow( mcmcMat )
  chain_a = mcmcMat[,"a"]
  chain_b = mcmcMat[,"b"]
  chain_s = mcmcMat[,"s"]

  # Posterior prediction:
  nevents = nrow(pop)
  nmmi = ncol(pop)
  nfat = nsamples * chainLength

  fat_rate <- matrix(NA, nrow = nmmi, ncol = nfat) 

  for ( i in 1:nmmi ) {
    mu0 <- exp(chain_b*mmi_list[i] + chain_a)
    rate_ = chain_s/mu0
    shape_ = chain_s + 0*mu0

    for (j in 1:nsamples) {
      fat_rate[i, seq(1+chainLength*(j-1),chainLength*j)]  = rtrunc(1,"gamma",b=1,rate=rate_, shape=shape_) 
    }    
  }
  fat_by_event <- pop %*% fat_rate # nevents x (nsamples*chainLength) 

  # HDI of estimated fatality
  dummy = compute_HDI_prob(fat_by_event)

  if ( !is.null(saveName) ) {
    saveRDS(fat_by_event, paste('fat_by_event_no_corr_', saveName,'.RDS',sep=""))
    saveRDS(dummy[[1]], paste('fatHDI_no_corr_', saveName,'.RDS',sep=""))
    saveRDS(dummy[[2]], paste('prob_mag_no_corr_', saveName,'.RDS',sep=""))
  }

}

#-------------------------------------------------------------------------------
# Main procedure
#------------------------------------------------------------------------------- 

#fileName = 'regression_wald_n40000mcmcMat.RDS'
fileName = 'regression_wald_19FebmcmcMat.RDS'
mcmcMat <- readRDS(fileName)

# read expo cat for reference
dummy <- read_expo_cat()
pop <- dummy[[1]]
#expo_cat <- dummy[[2]]
mmi_list <- dummy[[3]]

# fatality prediction for each of event
# 1. assuming perfect correlation
#compute_fat_corr( mcmcMat , pop, mmi_list, nsamples = 1, saveName = 'n1')

#compute_fat_corr( mcmcMat , pop, mmi_list, nsamples = 100, saveName = 'n100')

#compute_fat_no_corr( mcmcMat , pop, mmi_list, nsamples = 1, saveName = 'n1')

compute_fat_no_corr( mcmcMat , pop, mmi_list, nsamples = 100, saveName = 'n100')



