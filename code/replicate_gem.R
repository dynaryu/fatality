# Modified from Stan-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
source("DBDA2E-utilities.R")
require('qmap')
library(MASS)
require(ggplot2)
#library(scales)
#library(reshape2)


# copied from GEM_Exercise7.1.r
Irpinia.data <- read.csv("/Users/hyeuk/Project/GEM Physical Vulnerability Development Course documents/Session 7 - Uncertainty in deriving fragility functions/Worksheets & handouts/Data_Irpinia1980-Masonry.csv")
head(Irpinia.data)          # Take a quick look to see what variables are present
# remove zero total entry
row_to_remove <- Irpinia.data$Total == 0
Irpinia.data <- Irpinia.data[!row_to_remove,]

attach(Irpinia.data)        # "attach" allows us to access the columns by name ...
#> names(Irpinia.data)
#[1] "DS0"   "DS1"   "DS2"   "DS3"   "DS4"   "DS5"   "Total" "IM"    "class"
Irpinia.data$Mod.Or.Severe <- (DS3 + DS4 + DS5)
Irpinia.data$Mod.Or.Severe.Rate <- (DS3 + DS4 + DS5)/Total
#Irpinia.data$IMgal <- IM*100*9.81
Irpinia.data$logIMgal <- log(IM*100*9.81) # g -> gal
Irpinia.data$bin <- cut(Irpinia.data$logIMgal, breaks=seq(4.25,7.25,by=0.5), labels=seq(4.5,7.0,by=0.5))

detach(Irpinia.data)        # ... but should be used sparingly to avoid confusion!

# exploratory analysis
plot(Irpinia.data$logIMgal, Irpinia.data$Mod.Or.Severe.Rate)
plot(Irpinia.data$bin, Irpinia.data$Mod.Or.Severe.Rate) #boxplot
#table(Irpinia.data$bin) # showing frequency
#
#4.5   5 5.5   6 6.5   7 
# 20  36  52  36   8  12 

# the probit regression model
Irpinia.Model1A <- glm(Mod.Or.Severe.Rate ~ logIMgal,
                      family=binomial(link="probit"),weights=Total,
                      data=Irpinia.data)

logIM.limits <- extendrange(Irpinia.data$logIM)
logIM.grid <- seq(logIM.limits[1],logIM.limits[2],length.out=500)

#
#   Extract the predicted curve and intervals from the fitted model, 
#   on the probit scale.
#   
Model1A.predictions <- 
  predict(Irpinia.Model1A,
          newdata=data.frame(logIMgal=logIM.grid),
          se.fit=TRUE)

Model1A.predictions$pred <- pnorm(Model1A.predictions$fit)
df.Model1A <- data.frame(logIM.grid, Model1A.predictions$pred)
names(df.Model1A) <- c("logIM.grid","pred")

#
#   And plot (convert back to original PGA scale for compatibility with
#   plots in Session 7.2)
#
x11()
ggp0 <- ggplot(Irpinia.data)
ggp0 <- ggp0 + geom_point(aes(x=logIMgal, y=Mod.Or.Severe.Rate))
ggp0 <- ggp0 + geom_point(data = df.Model1A, aes(x=logIM.grid, y=pred)) 

# Bayesian approach
dataList = list(
  x = Irpinia.data$logIMgal, 
  y = Irpinia.data$Mod.Or.Severe,
  N = Irpinia.data$Total,
  Ndata = length(Irpinia.data$Mod.Or.Severe)
)

model_string = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  int<lower=0> y[Ndata] ; // no. of head 
  int<lower=1> N[Ndata] ; // no. of flip
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    for (n in 1:Ndata) {
       y[n] ~ binomial(N[n],Phi(a + b*x[n]));
    }    
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string)   

saveName = 'replicate_gem'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE


nChains = 4
iterSteps = 10000
burnInSteps =500

  # Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = 1 )

#print(stanFit, digits_summary=3, pars=c('a','b'), probs=c(0.025, 0.5, 0.975))
#traceplot(stanFit, pars=c('a','b'))

codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                 function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

#parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in c('a','b') ) {
diagMCMC( codaObject=codaSamples , parName=parName , 
          saveName= saveName , saveType= 'png')
}

mcmcMat = as.matrix(codaSamples,chains=TRUE)

#if ( !is.null(saveName) ) {
#  saveRDS(mcmcMat, file=paste(saveName,'mcmcMat.RDS',sep=''))
#}

chainLength = nrow( mcmcMat )
summaryInfo = NULL

chain_a = mcmcMat[,"a"]
chain_b = mcmcMat[,"b"]
#chain_c = mcmcMat[,"c"]
#chain_d = mcmcMat[,"d"]
#chain_s = mcmcMat[,"s"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
# summaryInfo = rbind( summaryInfo , 
#                    "s" = summarizePost( chain_s ) )
# summaryInfo = rbind( summaryInfo , 
#                    "c" = summarizePost( chain_c ) )
# summaryInfo = rbind( summaryInfo , 
#                    "d" = summarizePost( chain_d ) )

#show(summaryInfo)

#if ( !is.null(saveName) ) {
#  write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
#}

if (marPlot) { 
  # marginal distribution
  openGraph(width=6,height=3)
  layout( matrix( 1:2 , nrow=1) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
#  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(alpha)) #, main=paste("a=",parameters$a) )
#  histInfo = plotPost( chain_c , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(c)) #, main=paste("a=",parameters$a) )
#  histInfo = plotPost( chain_d , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(d)) #, main=paste("a=",parameters$a) )

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
  pairs( cbind( chain_a , chain_b )[plotIdx,] ,
         labels=c( expression(a) , expression(b) ) , 
         lower.panel=panel.cor , col="skyblue" )

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
  }
}

#nSlice = 13
#min_x = 4; max_x = 10
#mmi_list = seq(min_x,max_x,length=nSlice)

# comparison with data
if (comPlot) {

  y = dataList$y/dataList$N
  x = dataList$x

  #-----------------------------------------------------------------------------
  # Data with superimposed regression lines and noise distributions:
  openGraph()
  par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0))
  # Plot data values:
  postPredHDImass = 0.95
  xRang = max(x)-min(x)
  yRang = max(y)-min(y)
  xLimMult = 0.25
  yLimMult = 0.25
  xLim= c( 4.0 , 7.5 )
  yLim= c( 0 , 1.0 )
  plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
        xlab='IM' , ylab='Probability' , cex.lab=1.5 ,
        main=NULL) #paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
        #cex.main=1.33  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=50
  xComb = seq(xLim[1],xLim[2],length=501)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    lines( xComb , pnorm(chain_a[i] + chain_b[i]*xComb) , col="skyblue" )
  }

  # for (i in 1:nSlice) {
  #   yComb = seq(yHDI_mat[i,1],yHDI_mat[i,2],length=75)
  #   lines( mmi_list[i] - 0*yComb , yComb , col="red" , lwd=2 ) 
  # }

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"comPlot",sep=""), type=saveType)
  }
}


