graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
source("DBDA2E-utilities.R")

library(qmap)

genMCMC <- function(burnInSteps = 10000, nChains = 4, thinSteps = 1, 
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
                                  + burnInSteps ) , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = thinSteps )

  print(stanFit, digits_summary=3, pars=c('a','b','s'), 
  	    probs = c(0.025, 0.5, 0.975))


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
  #chain_s = mcmcMat[,"s"]

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

  #summaryInfo = rbind( summaryInfo , 
  #                   "s" = summarizePost( chain_s , 
  #                                            compVal=compValS , 
  #                                            ROPE=ropeS ) )
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

  #y = data[,yName]
  #x = data[,xName]
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = nrow( mcmcMat )
  chain_a = mcmcMat[,"a"]
  chain_b = mcmcMat[,"b"]
  #chain_s = mcmcMat[,"s"]

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



  # #-----------------------------------------------------------------------------
  # # Data with superimposed regression lines and noise distributions:
  # openGraph()
  # par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0))
  # # Plot data values:
  # postPredHDImass = 0.95
  # xRang = max(x)-min(x)
  # yRang = max(y)-min(y)
  # xLimMult = 0.25
  # yLimMult = 0.45
  # xLim= c( 4.0 , 10.5 )
  # yLim= c( 0 , 0.15 )
  # plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
  #       xlab='MMI' , ylab='Fatality rate' , cex.lab=1.5 ,
  #       main=NULL) #paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
  #       #cex.main=1.33  )
  # # Superimpose a smattering of believable regression lines:
  # nPredCurves=50
  # xComb = seq(xLim[1],xLim[2],length=501)
  # for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
  #   lines( xComb , exp(chain_a[i] + chain_b[i]*xComb) , col="skyblue" )
  # }

  # # dotted line of Sengara a= 0.622 and b=8.033
  # points(xComb, 10^(0.622*xComb-8.033), lty=2, lwd = 2)

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
  
  #if ( !is.null(saveName) ) {
  #   saveGraph( file=paste(saveName,"PostPred",sep=""), type=saveType)
  #} 

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

# synthetic data
# copied from Gamma Hurdle Models by Sean Anderson

set.seed(1)
N = 500
p=0.7
shape0 = 0.256
rate0 = 2.0

x <- 1:N
y <- rbinom(N, size = 1, prob = p)
y <- y * rgamma(N, shape = shape0, rate=rate0)
non_zero <- ifelse(y > 0, 1, 0)

d <- data.frame(x, y, non_zero)

# glm 
glm.fit1 <- glm(y~1, data=subset(d, non_zero==1), family=Gamma(log))
#glm.fit2 <- glm(y~1, data=subset(d, non_zero==1), family=Gamma())

shape_ <- gamma.shape(glm.fit1)
#gamma.shape(glm.fit2)

rate_ <- 1/mean(d$y[d$non_zero==1])*gamma.shape(glm.fit2)$alpha

#x <- 1:30
#y = rgamma(length(x), shape = 0.256, rate=3704.0)

# glm 
#d <- data.frame(x, y)
#m1 <- glm( y ~ 1, data = d, family = Gamma(link = log))

# MASS fitdistr
#fitdistr(y, 'gamma')

# glm
#fit <- glm(y ~ 1, family=Gamma())
#gamma.shape(fit)

#1/mean(y)*gamma.shape(fit)$alpha

#y <- rbinom(length(x), size = 1, prob = 0.7)
#y <- y * rgamma(length(x), shape = 0.256, rate=3704.0)
#non_zero <- ifelse(y > 0, 1, 0)
#d <- data.frame(x, y, non_zero)


#library(ggplot2)
#p <- ggplot(d, aes(x, y, colour = as.factor(non_zero))) + geom_point()
#print(p)


# N=30
# p = 0.7 # prob. of non-zero events 
# MMI = 6.0
# a = -21.1; b= 1.92 #alpha = 0.256 (shape)
# shape = 0.256
# rate = shape/exp(a+b*MMI) 
# scale = 1/rate
# set.seed(999)

# 

# shape = 0.4; rate = 1.0; p= 0.7
# rv <- rberngamma(n=50, prob=p, scale=1/rate, shape=shape)
# rv_nonzero = rv[rv!=0.0]
# fitdistr(rv_nonzero, "gamma")
# d <- data.frame(rv_nonzero)
# glm.fit <- glm( rv_nonzero ~ 1, data = d, family = Gamma())
# shape_ <- gamma.shape(glm.fit)
# rate_ <- gamma.shape(glm.fit)$alpha/mean(rv_nonzero)

dataList = list(
   y = d$y,
   Ndata = length(d$y)
  )

#saveRDS(dataList, 'dataList.RDS')

modelString = "
data {
  int<lower=1> Ndata ;
  vector<lower=0, upper=1>[Ndata] y ; // rate
}
parameters {
  real<lower=0> a ; // shape
  real<lower=0> b ; // rate
  real<lower=0, upper=1> s ; // non-zero prob.
}
model {
    for (n in 1:Ndata) {
    (y[n] == 0) ~ bernoulli(1-s);
    if (y[n] > 0)
      y[n] ~ gamma(a, b);
    }
}
"

stanDso <- stan_model( model_code = modelString) 

codaSamples <- genMCMC(modelObject = stanDso, dataList = dataList, burnInSteps = 20000, numSavedSteps = 5000, 
                        saveName = NULL, saveType = 'png') 
#saveRDS(codaSamples, file=paste(fileNameRoot,'codaSamples.RDS',sep=''))

nChains = 4
numSavedSteps = 000
thinSteps = 1
burnInSteps = 500

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                  + burnInSteps ) , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = thinSteps )

  
mcmcMat = as.matrix(codaSamples,chains=TRUE)
  
#saveRDS(mcmcMat, file=paste(fileNameRoot,'mcmcMat.RDS',sep=''))

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( codaSamples, 
                          saveName=NULL )
show(summaryInfo)

# Display posterior information:
plotMCMC( codaSamples , data=dat.df , xName="MMI" , yName="rate" , 
            pairsPlot=TRUE , showCurve=FALSE ,
            saveName=NULL , saveType=NULL )

