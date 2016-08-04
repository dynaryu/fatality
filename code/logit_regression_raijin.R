# Modified from Stan-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
source("DBDA2E-utilities.R")
require('qmap')
library(MASS)
#require(ggplot2)
#library(scales)
#library(reshape2)

# set.seed(999)
# N <- 100
# x <- 1:N
# z <- 0.01 * x - 5
# y <- sapply(1 / (1 + exp(-z)), function(p) {rbinom(1, 1, p)})

# dataList = list(
#   x = x, 
#   y = y,
#   Ndata = N
# )

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


model_berngamma = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // shape parameter
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    c ~ normal(0, 1.0e+2);
    d ~ normal(0, 1.0e+2);
    s ~ gamma(0.001, 0.001);

    for (n in 1:Ndata) {
    (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
    if (y[n] > 0)
       //y[n] ~ gamma(s, s / exp(a + b*x[n]));
       y[n] ~ gamma(s, s / inv_logit(a + b*x[n]));
    }    
}
"# close quote for modelString

saveName = 'wald_berngamma_n15k_logit'; saveType='eps'
#saveName = 'worden_berngamma_n15k'; saveType='eps'
#marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE
marPlot = FALSE; pairsPlot = FALSE; comPlot= FALSE
#dat <- read.csv('/Users/hyeuk/Project/fatality/DATA_WALD_COR_ROUND_12_Feb_2013.csv',header=0)
dat <- read.csv('~/scratch/fatality/DATA_WALD_COR_ROUND_12_Feb_2013.csv',header=0)
#dat <- read.csv('~/scratch/fatality/DATA_WORDEN_COR_ROUND_12_Feb_2013.csv',header=0)
names(dat) <- c("pop","fat","mmi","mmi_bin","id")
dat$rat <- dat$fat/dat$pop

dataList = list(
  x = dat$mmi, 
  y = dat$rat,
  Ndata = length(dat$mmi)
)


dummy <- read_expo_cat(fileName='~/scratch/fatality/EXPO_CAT_2007_12.csv')
pop <- dummy[[1]]
expo_cat <- dummy[[2]]
mmi_list <- dummy[[3]]

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_berngamma)   

nChains = 4
iterSteps = 15000
burnInSteps =5000

  # Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = 1 )

#print(stanFit, digits_summary=3, pars=c('alpha','beta'), probs=c(0.025, 0.5, 0.975))

#traceplot(stanFit, pars=c('alpha','beta'))

codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                 function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

#parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in c('a','b','s','c','d') ) {
diagMCMC( codaObject=codaSamples , parName=parName , 
          saveName= saveName , saveType= 'png')
}

mcmcMat = as.matrix(codaSamples,chains=TRUE)

if ( !is.null(saveName) ) {
  saveRDS(mcmcMat, file=paste(saveName,'mcmcMat.RDS',sep=''))
}

chainLength = nrow( mcmcMat )

summaryInfo = NULL

chain_a = mcmcMat[,"a"]
chain_b = mcmcMat[,"b"]
chain_c = mcmcMat[,"c"]
chain_d = mcmcMat[,"d"]
chain_s = mcmcMat[,"s"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
summaryInfo = rbind( summaryInfo , 
                   "s" = summarizePost( chain_s ) )
summaryInfo = rbind( summaryInfo , 
                   "c" = summarizePost( chain_c ) )
summaryInfo = rbind( summaryInfo , 
                   "d" = summarizePost( chain_d ) )

#show(summaryInfo)

if ( !is.null(saveName) ) {
  write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
}

if (marPlot) { 
  # marginal distribution
  openGraph(width=15,height=3)
  layout( matrix( 1:5 , nrow=1) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(alpha)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_c , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(c)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_d , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(d)) #, main=paste("a=",parameters$a) )

  if ( !is.null(saveName) ) {
  saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }
}

#-----------------------------------------------------------------------------
if ( pairsPlot ) {

  # Plot the parameters pairwise, to see correlations:
  nPtToPlot = 1000
  plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
  openGraph()
  panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
  }
  pairs( cbind( chain_a , chain_b , chain_s, chain_c, chain_d )[plotIdx,] ,
         labels=c( expression(a) , expression(b) , 
                   expression(s) , expression(c) , expression(d) ) , 
         lower.panel=panel.cor , col="skyblue" )

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
  }
}

#nSlice = 13
#min_x = 4; max_x = 10
#mmi_list = seq(min_x,max_x,length=nSlice)

nSlice = length(mmi_list)
nsamples = 100

# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
#for ( chainIdx in 1:chainLength ) {

nfat = nsamples * chainLength  
fat_rate <- matrix(NA, nrow = nSlice, ncol = nsamples * chainLength  )
fat_rate_HDI <- matrix(NA, nrow=nSlice, ncol=2)

for ( j in 1:nSlice ) {
  for (k in 1:chainLength) {
    
    lin1 <- chain_c[k] + chain_d[k]*mmi_list[j]
    p1 <- 1/(1+exp(-lin1))
    
    #y <- rbinom(nsamples, size = 1, prob = 1-p1)
        
    lin0 <- chain_a[k] + chain_b[k]*mmi_list[j]
    mu0 <- 1/(1+exp(-lin0))
    #mu0 <- exp(chain_a[k] + chain_b[k]*mmi_list[j])
    shape_ = chain_s[k]
    scale_ = mu0/chain_s[k]
    y <- rberngamma(nsamples, 1-p1, scale=scale_, shape=shape_)
    
    fat_rate[j, seq(1+nsamples*(k-1),nsamples*k)] = y
  }
  fat_rate_HDI[j,] <- HDIofMCMC(fat_rate[j,])  
}

fat_by_event <- pop %*% fat_rate # nevents x (nsamples*chainLength) 

fat_by_event0 <- fat_by_event[42,] # 200006041628, Bengkulu, Indonesia) INDEX: 42

# HDI of estimated fatality
dummy = compute_HDI_prob(fat_by_event)
fatHDI = dummy[[1]]
prob_mag = dummy[[2]]

count_order <- count_order_match(prob_mag, expo_cat)
show(count_order)

if ( !is.null(saveName) ) {
  saveRDS(fat_rate, file=paste(saveName,'fat_rate.RDS',sep=''))
  saveRDS(fat_rate_HDI, file=paste(saveName,'fat_rate_HDI.RDS',sep=''))

  #saveRDS(fat_by_event, paste(saveName, 'fat_by_event.RDS',sep=""))
  saveRDS(fat_by_event0, paste(saveName, 'fat_by_event_42.RDS',sep=""))
  saveRDS(fatHDI, paste(saveName, 'fatHDI.RDS',sep=""))
  saveRDS(prob_mag, paste(saveName, 'prob_mag.RDS',sep=""))
}

# comparison with data
if (comPlot) {

  y = dat$rat
  x = dat$mmi

  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  openGraph()
  par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0))
  # Plot data values:
  postPredHDImass = 0.95
  xRang = max(x)-min(x)
  yRang = max(y)-min(y)
  xLimMult = 0.25
  yLimMult = 0.45
  xLim= c( 4.0 , 10.5 )
  yLim= c( 0 , 0.15 )
  plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
        xlab='MMI' , ylab='Fatality rate' , cex.lab=1.5 ,
        main=NULL) #paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
        #cex.main=1.33  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=50
  xComb = seq(xLim[1],xLim[2],length=501)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    lines( xComb , exp(chain_a[i] + chain_b[i]*xComb) , col="skyblue" )
  }

  for (i in 1:nSlice) {
    yComb = seq(yHDI_mat[i,1],yHDI_mat[i,2],length=75)
    lines( mmi_list[i] - 0*yComb , yComb , col="red" , lwd=2 ) 
  }

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"comPlot",sep=""), type=saveType)
  }
}


