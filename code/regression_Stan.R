# Modified from Stan-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
source("DBDA2E-utilities.R")
library(MASS)
#require(ggplot2)
#library(scales)

# compute rate and shape parameters for gamma distribution
#calculate_gamma_parameters <- function(mu0, sd0){
#    tau = 1/sd0^2
#    r <- mu0*tau # rate
#    s <- mu0^2*tau # shape
#return(list("r"=r,"s"=s))
#}

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

genMCMC <- function(adaptSteps = 5000, burnInSteps = 10000, nChains = 4, thinSteps = 1, 
                    numSavedSteps = 20000, modelObject = stanDso, dataList = dataList, saveName = saveName, saveType = 'png') {

  #writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let Stan do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  #parameters = c( "a" ,  "b" ,  "k" )
  #adaptSteps = 5000  # Number of steps to "tune" the samplers
  #burnInSteps = 10000
  #nChains = 4 
  #thinSteps = 1
  #numSavedSteps = 50000

  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                  + adaptSteps + burnInSteps ) , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = thinSteps )

  codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                     function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

  #parameterNames = varnames(codaSamples) # get all parameter names
  for ( parName in c('a','b','s') ) {
    diagMCMC( codaObject=codaSamples , parName=parName , 
              saveName= saveName , saveType= saveType )
  }
  return( codaSamples )
}

# summary MCMC
smryMCMC = function(  codaSamples , 
                      compValA=NULL , ropeA=NULL , 
                      compValB=NULL , ropeB=NULL , 
                      compValS=NULL , ropeS=NULL , 
                      saveName=NULL ) {

  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = nrow( mcmcMat )
  chain_a = mcmcMat[,"a"]
  chain_b = mcmcMat[,"b"]
  chain_s = mcmcMat[,"s"]

  # summaryInfo
  summaryInfo = NULL
  summaryInfo = rbind( summaryInfo , 
                     "a" = summarizePost( chain_a , 
                                              compVal=compValA , 
                                              ROPE=ropeA ) )
  summaryInfo = rbind( summaryInfo , 
                     "b" = summarizePost( chain_b , 
                                              compVal=compValB , 
                                              ROPE=ropeB ) )
  summaryInfo = rbind( summaryInfo , 
                     "s" = summarizePost( chain_s , 
                                              compVal=compValS , 
                                              ROPE=ropeS ) )
  #show(summaryInfo)

  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

# Plot MCMC
plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     compValA=NULL , ropeA=NULL , 
                     compValB=NULL , ropeB=NULL , 
                     compValS=NULL , ropeS=NULL , 
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="png" ) {

  y = data[,yName]
  x = data[,xName]
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = nrow( mcmcMat )
  chain_a = mcmcMat[,"a"]
  chain_b = mcmcMat[,"b"]
  chain_s = mcmcMat[,"s"]

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
  pairs( cbind( chain_a , chain_b , chain_s )[plotIdx,] ,
         labels=c( expression(a) , expression(b) , 
                   expression(s)  ) , 
         lower.panel=panel.cor , col="skyblue" )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
  }
  }

  #-----------------------------------------------------------------------------
  # Marginal histograms:
  # Set up window and layout:
  #nPtToPlot = 1000
  #plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
  openGraph(width=9,height=3)
  layout( matrix( 1:3 , nrow=1) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE ,
                     compVal=NULL , ROPE=NULL ,
                     xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE ,
                     compVal=NULL , ROPE=NULL ,
                     xlab=bquote(b)) #, main=paste("b=",parameters$b) )
  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE ,
                     compVal=NULL , ROPE=NULL ,
                     xlab=bquote(alpha)) #, main=paste("k=",parameters$s) )
  if ( !is.null(saveName) ) {
  saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }

  # compare HDIofCDF vs HDIofMCMC

    
  #HDI_CDF <- matrix( NA, nrow = chainLength, ncol = 2)
  #HDI_MCMC <- matrix( NA, nrow = chainLength, ncol = 2)
  
  #for ( i in 1:chainLength ) {

  #    mu0 <- exp(chain_b[i]*x_val+chain_a[i])

  #    yHDI = HDIofICDF( qgamma , credMass=0.95, rate=chain_s[i]/mu0, shape=chain_s[i])
  #    HDI_CDF[i,] = yHDI

      #tmp = rtrunc(5000,"gamma",b=1,rate=chain_s[i]/mu0,shape=chain_s[i])
      #HDI_MCMC[i, ] = HDIofMCMC(tmp)
  #} 

  #mu0 <- exp(chain_b*x_val+chain_a)

  #yHDI = HDIofICDF( qgamma , credMass=0.95, rate=chain_s/mu0, shape=chain_s)



      #yComb = seq(yHDI[1],yHDI[2],length=75)



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

  # dotted line of Sengara a= 0.622 and b=8.033
  points(xComb, 10^(0.622*xComb-8.033), lty=2, lwd = 2)

  #nSlice = 13
  #min_x = 4; max_x = 10
  #mmi_list = seq(min_x,max_x,length=nSlice)
  #nsamples = 100
  #nPredCurves=30

  #PostPred <- matrix( NA, nrow = nSlice, ncol = 2)
  # # Generate posterior predicted y values.
  # # This gets only one y value, at each x, for each step in the chain.
  # #for ( chainIdx in 1:chainLength ) {

  # nfat = nsamples * chainLength  

  # for ( j in 1:nSlice ) {

  #   fat_rate <- numeric(length = nfat)

  #   #for ( i in seq(from=1,to=chainLength) ) {

  # #  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {

  #   mu0 <- exp(chain_b*mmi_list[j] + chain_a)
  #   rate_ = chain_s/mu0
  #   shape_ = chain_s + 0*mu0

  #   for (k in 1:nsamples) {
  #     fat_rate[seq(1+chainLength*(k-1),chainLength*k)]  = rtrunc(1,"gamma",b=1,rate=rate_, shape=shape_) 
  #   }

  #   yHDI = HDIofMCMC(fat_rate)

  #   yComb = seq(yHDI[1],yHDI[2],length=75)

  #   lines( mmi_list[j] - 0*yComb , yComb , col="red" , lwd=2 )

  # }
  
  if ( !is.null(saveName) ) {
     saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  } 

  # # # Superimpose some vertical distributions to indicate spread:
  # # nSlice = 13
  # # min_x = 4; max_x = 10
  # # curveXpos = seq(min_x,max_x,length=nSlice)
  # # curveWidth = (max_x-min_x)/(nSlice+2)
  # # for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
  # #   for ( j in 1:length(curveXpos) ) {

  # #     mu0 <- exp(chain_b[i]*curveXpos[j]+chain_a[i])
  # #     #sig0 <- chain_k[i]*mu0

  # #     #r0 <- mu0/(sig0^2)               
  # #     #s0 <- (mu0^2)/(sig0^2)

  # #     yHDI = HDIofICDF( qgamma , credMass=postPredHDImass, rate=chain_s[i]/mu0, shape=chain_s[i])

  # #     yComb = seq(yHDI[1],yHDI[2],length=75)
  # #     #xVals = dgamma( yComb , rate=chain_s[i]/mu0, shape=chain_s[i])
  # #     #xVals = curveWidth * xVals / max(xVals)
  # #     #lines( curveXpos[j] - xVals , yComb , col="skyblue" )
  # #     #lines( curveXpos[j] - 0*xVals , yComb , col="skyblue" , lwd=2 )
  # #     lines( curveXpos[j] - 0*yComb , yComb , col="skyblue" , lwd=2 )
  # #   }
  # # }
  # # replot the data, in case they are obscured by lines:
  # points( x , y , cex=1.5 )
}

# generate synthetic data
gen_synthetic_data <- function(a, b, s, npts){

  # mu ~ exp(a+bx)
  # sig ~ k*mu (shape is 1/k2)
  # y ~ gamma(1/k2, exp(a+b*x))

  # MMI
  xmin = 4
  xmax = 10  
  
  set.seed(999)
  
  x <- runif(npts, min = xmin, max = xmax)
  y <- numeric(length=npts)

  mu <- exp( b*x + a)
  shape <- s #1/k^2
  #sig <- k * mu
  #tmp <- calculate_gamma_parameters(mu,sig)
  y <- rtrunc(npts,"gamma", b=1, shape= shape, rate= shape / mu) 

  dataList = list(
    x = x, 
    y = y,
    Ndata = npts
  )
  return(dataList)
}

exploratory_analysis <- function(dat){

  # MMI vs fatality rate
  # box plot bin 
  #dat$mmi_bin <- cut(dat$mmi, seq(3.75, 10.5, 0.5), labels=seq(4.0, 10.0, 0.5))
  #
  # exploratory analysis
  #
  # given.n <- function(x) {
  #   return(c(y=mean(x), label=length(x)))
  # }

  # X11()
  # ggplot(dat, aes(factor(mmi_bin), rat)) + 
  # geom_boxplot() +
  # stat_summary(fun.data = given.n, geom = "text") +
  # xlab("MMI") + 
  # ylab("Fatality rate")

  # saveGraph(file='boxplot_dat',type='png')

  x11()
  #boxplot(dat$rat ~ dat$mmi_bin)
  #means <- tapply(dat$rat, dat$mmi_bin, mean)

  #points(means, col='red', pch=18)

  # dat8.5 <- dat$rat[dat$mmi_bin == 8.5]

  sd_ <- tapply(dat$rat, dat$mmi_bin, sd)

  mean_ <- tapply(dat$rat, dat$mmi_bin, mean)

  X11()
  qplot(names(sd_), log(mean_)) + xlab("MMI") + ylab("log mean of fatality rate")
  #dev.copy2eps(file='log_mean_fatality.eps')

  X11()
  qplot(names(sd_), sd_/mean_) + xlab("MMI") + ylab("COV of fatality rate")
  # dev.copy2eps(file='cov_fatality.eps')

}

develop_model_using_data <- function(csvfileName, fileNameRoot, graphFileType) {

  #------------------------------------------------------------------------------
  #dat <- read.csv('/Users/hyeuk/Project/fatality/DATA_WALD_COR_ROUND_12_Feb_2013.csv',header=0)
  dat <- read.csv(csvfileName, header=0)
  names(dat) <- c("pop","fat","mmi","mmi_bin","id")
  dat$fat[dat$fat==0] = 1.0e-3 # very small number instead of 0 fatality
  dat$rat <- dat$fat/dat$pop

  dataList = list(
      x = dat$mmi, 
      y = dat$rat,
      Ndata = length(dat$mmi)
    )

  dat.df <- data.frame(dat$mmi, dat$rat)
  names(dat.df) <- c("MMI","rate")

  # exploratory analysis
  #exploratory_analysis(dat)

  # GLM/ gamma
  gm <- glm(rate ~ MMI, Gamma(link=log), dat.df, start=c(0,0))
  x11(); par(mfrow=c(2,2)); plot(gm)

  myshape <- gamma.shape(gm) # more accurate estimation of shape parameter 
  # see http://civil.colorado.edu/~balajir/CVEN6833/lectures/GammaGLM-01.pdf
  gampred <- predict(gm, type = "response", se = T, dispersion = 1/myshape$alpha)
  summary(gm, dispersion=1/myshape$alpha)

  #fileNameRoot = 'regression_wald'
  #graphFileType = 'eps'

  codaSamples <- genMCMC(modelObject = stanDso, dataList = dataList, 
                        saveName = fileNameRoot, saveType = 'png') 
  #saveRDS(codaSamples, file=paste(fileNameRoot,'codaSamples.RDS',sep=''))
  
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  
  saveRDS(mcmcMat, file=paste(fileNameRoot,'mcmcMat.RDS',sep=''))

  #------------------------------------------------------------------------------- 
  # Get summary statistics of chain:
  summaryInfo = smryMCMC( codaSamples, 
                          saveName=fileNameRoot )
  show(summaryInfo)

  # Display posterior information:
  plotMCMC( codaSamples , data=dat.df , xName="MMI" , yName="rate" , 
            pairsPlot=TRUE , showCurve=FALSE ,
            saveName=fileNameRoot , saveType=graphFileType )

  #openGraph()
  #traceplot(stanFit,pars=c("a"))
  #openGraph()
  #plot(stanFit,pars=c("a","b","k"))

  #------------------------------------------------------------------------------- 
  # estimate fatality
  #estimate_fatality_wald = estimate_fatality_HDI( codaSamples_wald , pop , mmi_list )
  #plotCI(estimate_fatality_wald, observed, eqID, saveName = fileNameRoot, saveType=graphFileType)

  # comparison with Kishor's model
  #params_PAGER = c(14.05, 0.17, 1.74) # USGS OpenFile report 2009-1136
  #params_PAGER = c(14.05, 0.17, 2.15) # EQ spectra (Jaiswal and Wald, 2010)
  #params_PAGER = c(13.249,  0.151, 1.641) # http://pubs.usgs.gov/of/2009/1136/pdf/PAGER%20Implementation%20of%20Empirical%20model.xls(Jaiswal, 2009)
  #estimate_PAGER <- estimate_fatality_PAGER(params_PAGER, pop, mmi_list)
  #plotCI(estimate_PAGER, observed)

  #med_kishor = 13.249; beta_kishor = 0.151;  # v2.0
  #x_trev, y_trev = 1.03, 10.75 # Trevo Allen

  #14.05 0.17 2.15

  #------------------------------------------------------------------------------
}

#-------------------------------------------------------------------------------
# Main procedure
#------------------------------------------------------------------------------- 

#-----------------------------------------------------------------------------
# THE MODEL.
fatality_model = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a ; // parameter of linear predictor
  real b ; // (a+bx)
  real<lower=0> s ; // shape parameter 
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    s ~ gamma(0.001, 0.001);

    /*
    for (i in 1:Ndata) {
        y[i] ~ gamma(s[i], r[i]) T[0,1];
    }
    */

    y ~ gamma(s, s ./ exp(a+ b*x)); // T[0,1];
}
"# close quote for modelString

fatality_model2 = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a ; // parameter of linear predictor
  real b ; // (a+bx)
  /*
  real<lower=0> s ; // shape parameter
  */ 
  vector<lower=0>[Ndata] s ; // fatality rate
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    s ~ gamma(0.001, 0.001);

    for (i in 1:Ndata) {
        y[i] ~ gamma( s[i], s[i]/(a+b*x[i])) T[0,1];
    }
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = fatality_model) 

flag_synthetic = FALSE
flag_data = TRUE

if (flag_synthetic) {

  # synthetic data for model validation
  param <- data.frame(-15.0, 1.0, 1.0)
  names(param) <- c("a","b","s")

  dataList <- gen_synthetic_data(a = param$a, b = param$b, s = param$s, npts = 50)
  data_synth <- data.frame(dataList$x, dataList$y)
  names(data_synth) <- c("MMI","rate")

  # GLM/ gamma
  #gm_synth_gamma <- glm(rate ~ MMI, Gamma(link=log), data_synth, start=c(0,0))
  #x11(); par(mfrow=c(2,2)); plot(gm_synth_gamma)

  #myshape <- gamma.shape(gm_synth_gamma) # more accurate estimation of shape parameter 
  # see http://civil.colorado.edu/~balajir/CVEN6833/lectures/GammaGLM-01.pdf
  #gampred <- predict(gm_synth_gamma, type = "response", se = T, dispersion = 1/myshape$alpha)
  #summary(gm_synth_gamma, dispersion=1/myshape$alpha)
  # shape parameter is 1/K^2
  # dispersion parameter is K 

  #gm_synth <- glm(rate ~ MMI, quasi(link=log, variance="mu^2"), data_synth, start = c(0,0))
  #x11(); par(mfrow=c(2,2)); plot(gm_synth)

  # GLM/ normal
  #gm_synth_normal <- glm(rate ~ MMI, gaussian(link=log), data_synth)
  #gm_synth_norm <- glm(rate ~ MMI, gaussian(link=log), data_synth, start=c(0,0))
  #x11(); par(mfrow=c(2,2)); plot(gm_synth_normal)

  codaSamples <- genMCMC(modelObject = stanDso, dataList = dataList, 
                        saveName = 'synthetic', saveType = 'png') 

  saveRDS(codaSamples, file='codaSamples.RDS')

  #------------------------------------------------------------------------------- 
  # Get summary statistics of chain:
  fileNameRoot = 'regression_Stan'
  graphFileType = 'png'
  #summaryInfo = smryMCMC( codaSamples , 
  #                        saveName=fileNameRoot )
  #show(summaryInfo)

  # Display posterior information:
  #plotMCMC( codaSamples , data=data_synth , xName="MMI" , yName="rate" , 
  #          pairsPlot=TRUE , showCurve=FALSE ,
  #          saveName=fileNameRoot , saveType=graphFileType )

#estimate_fatality = estimate_fatality_HDI( codaSamples , pop , mmi_list )
#plotCI(estimate_fatality, observed, eqID, saveName = fileNameRoot, saveType='eps')


}

# #------------------------------------------------------------------------------
# # Another data set:
# dataList_pt100 <- gen_synthetic_data(param$a, param$b, param$k, 100)

if (flag_data) {

csvfileName_wald = '/Users/hyeuk/Project/fatality/DATA_WALD_COR_ROUND_12_Feb_2013.csv'
develop_model_using_data(csvfileName_wald,'regression_wald_19Feb','eps')

csvfileName_worden = '/Users/hyeuk/Project/fatality/DATA_WORDEN_COR_ROUND_12_Feb_2013.csv'
develop_model_using_data(csvfileName_worden, 'regression_worden_19Feb', 'eps')

}
