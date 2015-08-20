# collection of useful functions
# mainly copied from DBDA2E-utilities.R

library(qmap)
library(rstan)
library(coda)

# (0,1) -> (-inf, inf)
logit <- function(p) {
  log(p/(1-p))
}

# [-inf, inf] -> (0, 1)
inv_logit <- function(x) {
  1/(1+exp(-x))
}

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

summarizePost = function( paramSampleVec , 
                          compVal=NULL , ROPE=NULL , credMass=0.95 ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round( effectiveSize( paramSampleVec ) , 1 )
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                    / length( paramSampleVec ) )
  } else {
    compVal=NA
    pcgtCompVal=NA
  }
  if ( !is.null(ROPE) ) {
    pcltRope = ( 100 * sum( paramSampleVec < ROPE[1] ) 
                 / length( paramSampleVec ) )
    pcgtRope = ( 100 * sum( paramSampleVec > ROPE[2] ) 
                 / length( paramSampleVec ) )
    pcinRope = 100-(pcltRope+pcgtRope)
  } else { 
    ROPE = c(NA,NA)
    pcltRope=NA 
    pcgtRope=NA 
    pcinRope=NA 
  }  
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam , 
             ESS=mcmcEffSz ,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] , 
             CompVal=compVal , PcntGtCompVal=pcgtCompVal , 
             ROPElow=ROPE[1] , ROPEhigh=ROPE[2] ,
             PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) )
}

read_expo_cat <- function(fileName = '/Users/hyeuk/Project/fatality/data/EXPO_CAT_2007_12.csv') {

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

  # 0-10(1), 10-100(2), 100-1000(3), 1000-10000(4)
  mag_observed <- floor(log10(observed)) + 1
  mag_observed[is.infinite(mag_observed)] = 1 # 0 range

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
  # 0-10(1), 10-100(2), 100-1000(3), 1000-10000(4)

  prob_mag <- matrix(NA, nrow = nevents, ncol = max(order_mag))

  for (i in 1:nevents) {

    temp = fat_by_event[i, ]
    fatHDI[i,] = HDIofMCMC(temp)

    cprob <- numeric(length = order_mag[i])
    for (j in 1:order_mag[i]-1) {
        #cprob[j] <- sum(temp < 10^(j-1)) / nfat
        cprob[j] <- sum(temp < 10^j) / nfat
    } 
    cprob[order_mag[i]] <- 1.0

    prob_mag[i,1:(order_mag[i])] <- c(cprob[1], diff(cprob))

  }

  return(list(fatHDI, prob_mag))
}

estimate_fat_rate_HDI_bernlnorm <- function(nsamples) {

  nSlice = length(mmi_list)
  nfat = nsamples * chainLength  
  fat_rate <- matrix(NA, nrow = nSlice, ncol = nsamples * chainLength  )
  fat_rate_HDI <- matrix(NA, nrow=nSlice, ncol=2)

  for ( j in 1:nSlice ) {

    mmi <- mmi_list[j]

    for (k in 1:chainLength) {
      prob <- inv_logit(chain$c[k] + chain$d[k]*mmi)

      meanlog <- chain$a[k] + chain$b[k]*mmi
      sdlog <- chain$s[k]
      y <- rbernlnorm(nsamples, 1-prob, meanlog, sdlog)   
      fat_rate[j, seq(1+nsamples*(k-1), nsamples*k)] = y

      #y <- rtrunc(nsamples, "bernlnorm", b=1, meanlog=meanlog, sdlog=sdlog, prob=1-prob)
      #mu <- inv_logit(chain_a[k] + chain_b[k]*mmi_list[j])
      #mu <- exp(chain_a[k] + chain_b[k]*mmi_list[j])
      #shape <- chain_s[k]
      #scale <- mu*shape
      #y <- rtrunc(nsamples, "berngamma", b=1, scale=scale, shape=shape, prob=1-prob)
    }
    fat_rate_HDI[j,] <- HDIofMCMC(fat_rate[j,])  
  }

  # comparison with data
  fat_rate_HDI <- as.data.frame(fat_rate_HDI)
  names(fat_rate_HDI) <- c("L", "U")
  fat_rate_HDI$mmi_bin <- mmi_list

  return(list(fat_rate, fat_rate_HDI))
}

estimate_fat_rate_HDI_berngamma_log <- function(nsamples) {

  nSlice = length(mmi_list)
  nfat = nsamples * chainLength  
  fat_rate <- matrix(NA, nrow = nSlice, ncol = nsamples * chainLength  )
  fat_rate_HDI <- matrix(NA, nrow=nSlice, ncol=2)

  for ( j in 1:nSlice ) {

    mmi <- mmi_list[j]

    for (k in 1:chainLength) {
      prob <- inv_logit(chain$c[k] + chain$d[k]*mmi)

      #mu <- inv_logit(chain_a[k] + chain_b[k]*mmi_list[j])
      mu <- exp(chain$a[k] + chain$b[k]*mmi)
      shape <- chain$s[k]
      scale <- mu/shape
      y <- rberngamma(nsamples, scale=scale, shape=shape, prob=1-prob)
      fat_rate[j, seq(1+nsamples*(k-1), nsamples*k)] = y
    }
    fat_rate_HDI[j,] <- HDIofMCMC(fat_rate[j,])  
  }

  # comparison with data
  fat_rate_HDI <- as.data.frame(fat_rate_HDI)
  names(fat_rate_HDI) <- c("L", "U")
  fat_rate_HDI$mmi_bin <- mmi_list

  return(list(fat_rate, fat_rate_HDI))
}

estimate_fat_rate_HDI_berngamma_logit <- function(nsamples) {

  nSlice = length(mmi_list)
  nfat = nsamples * chainLength  
  fat_rate <- matrix(NA, nrow = nSlice, ncol = nsamples * chainLength  )
  fat_rate_HDI <- matrix(NA, nrow=nSlice, ncol=2)

  for ( j in 1:nSlice ) {

    mmi <- mmi_list[j]

    for (k in 1:chainLength) {
      prob <- inv_logit(chain$c[k] + chain$d[k]*mmi)
      mu <- inv_logit(chain$a[k] + chain$b[k]*mmi)
      shape <- chain$s[k]
      scale <- mu/shape
      y <- rberngamma(nsamples, scale=scale, shape=shape, prob=1-prob)
      fat_rate[j, seq(1+nsamples*(k-1), nsamples*k)] = y
    }
    fat_rate_HDI[j,] <- HDIofMCMC(fat_rate[j,])  
  }

  # comparison with data
  fat_rate_HDI <- as.data.frame(fat_rate_HDI)
  names(fat_rate_HDI) <- c("L", "U")
  fat_rate_HDI$mmi_bin <- mmi_list

  return(list(fat_rate, fat_rate_HDI))
}

count_order_match <- function(prob_mag, expo_cat) {

  idx_max <- apply(prob_mag, 1, function(x) max(which(x == max(x, na.rm = TRUE))))

  #k <- apply(prob_mag, 1, function(x) max(which(x == min(x, na.rm = TRUE))))

  tf = expo_cat$mag_obs == idx_max

  nmatch <- numeric(length=max(expo_cat$mag_obs))
  for (i in 1:length(nmatch)) {
    nmatch[i] <- sum( tf & (expo_cat$mag_obs == i) )
  } 
  return(nmatch)
}
