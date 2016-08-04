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
require(betareg)

flag_gem = 0
flag_bayesian = 0
flag_uncertainty = 0
flag_uncertainty2 = 1

###############################################################################
if (flag_gem) {
# copied from GEM_Exercise7.1.r
Irpinia.data <- read.csv("/Users/hyeuk/Project/GEM Physical Vulnerability Development Course documents/Session 7 - Uncertainty in deriving fragility functions/Worksheets & handouts/Data_Irpinia1980-Masonry.csv")
head(Irpinia.data)          # Take a quick look to see what variables are present
# remove zero total entry
row_to_remove <- Irpinia.data$Total == 0
Irpinia.data <- Irpinia.data[!row_to_remove,]
# not sure why removal of Nan makes difference.

#> names(Irpinia.data)
#[1] "DS0"   "DS1"   "DS2"   "DS3"   "DS4"   "DS5"   "Total" "IM"    "class"
Irpinia.data$Mod.Or.Severe <- (Irpinia.data$DS3 + Irpinia.data$DS4 + Irpinia.data$DS5)
Irpinia.data$Mod.Or.Severe.Rate <- Irpinia.data$Mod.Or.Severe/Irpinia.data$Total
Irpinia.data$logIMgal <- log(Irpinia.data$IM*100*9.81) # g -> gal
# noramlize IM
mean_im <- mean(Irpinia.data$logIMgal)
#sd_im <- sd(Irpinia.data$logIMgal)
Irpinia.data$logIMgal.centred <- Irpinia.data$logIMgal - mean_im
Irpinia.data$bin <- cut(Irpinia.data$logIMgal, breaks=seq(4.25,7.25,by=0.5), labels=seq(4.5,7.0,by=0.5))

# exploratory analysis
plot(Irpinia.data$logIMgal, Irpinia.data$Mod.Or.Severe.Rate)
plot(Irpinia.data$bin, Irpinia.data$Mod.Or.Severe.Rate) #boxplot
#table(Irpinia.data$bin) # showing frequency
#
#4.5   5 5.5   6 6.5   7 
# 20  36  52  36   8  12 

# the probit regression model
Irpinia.Model1A <- glm(Mod.Or.Severe.Rate ~ logIMgal.centred,
                      family=binomial(link="probit"),weights=Total,
                      data=Irpinia.data)
}
# logIM.limits <- extendrange(Irpinia.data$logIMgal)
# logIM.grid <- seq(logIM.limits[1],logIM.limits[2],length.out=500)

# #
# #   Extract the predicted curve and intervals from the fitted model, 
# #   on the probit scale.
# #   
# Model1A.predictions <- 
#   predict(Irpinia.Model1A,
#           newdata=data.frame(logIMgal=logIM.grid),
#           se.fit=TRUE)

# Model1A.predictions$pred <- pnorm(Model1A.predictions$fit)
# df.Model1A <- data.frame(logIM.grid, Model1A.predictions$pred)
# names(df.Model1A) <- c("logIM.grid","pred")

# #
# #   And plot (convert back to original PGA scale for compatibility with
# #   plots in Session 7.2)
# #
# x11()
# ggp0 <- ggplot(Irpinia.data)
# ggp0 <- ggp0 + geom_point(aes(x=logIMgal, y=Mod.Or.Severe.Rate))
# ggp0 <- ggp0 + geom_point(data = df.Model1A, aes(x=logIM.grid, y=pred)) 
#################################################

# dummy data set 1
set.seed(999)
n.reps <- 1
n.repeated <- 10
n <- 50
gain <- 1.0
threshold <- 0.0
sigmay <- 1.0
sx <- 0.5  # measurement error

# true covariate values
xtrue <- runif(n, -3, 3)
#y <- rnorm(n, alpha_ + beta_*xtrue, sigmay_)  # alpha=0, beta=1, sdy=1

# random subset to perform repeat covariate measurements
#which.repeated <- sample(n, n.repeated)
xerr <- xtrue + rnorm(n + (n.repeated * (n.reps - 1)), 0, sx)

logit_y <- 1/(1+exp(-gain*(xtrue-threshold)))

# beta_a <- mu*ph
ph <- 1.0
beta_a <- logit_y*ph
beta_b <- (1-logit_y)*ph
y2 <- rbeta(n, beta_a, beta_b)  # alpha=0, beta=1, sdy=1
y2[y2==1.0] = 1.0-1.0e-4

#glm_err <- glm(y2~xerr, family=beta, link=logit)
#print(glm_err)
# Call:
# lm(formula = y ~ x)

# Coefficients:
# (Intercept)            x  
#      1.4175       0.2184  

###############################################################################
# Bayesian approach
###############################################################################

# dataList = list(
#   x = Irpinia.data$logIMgal.centred, 
#   y = Irpinia.data$Mod.Or.Severe,
#   N = Irpinia.data$Total,
#   Ndata = length(Irpinia.data$Mod.Or.Severe)
# )

###############################################################################
if (flag_bayesian) {

dataList = list(
  x = xtrue, 
  y = y2,
  Ndata = length(xtrue)
)

dataList_error = list(
  x = xerr, 
  y = y2,
  Ndata = length(xerr)
)

model_string = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;// x
  vector<lower=0, upper=1>[Ndata] y ;// y 
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> ph; // sigma
}
model {
    vector<lower=0,upper=1>[Ndata] mu;
    vector<lower=0>[Ndata] beta_a;
    vector<lower=0>[Ndata] beta_b;

    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5);
    ph ~ cauchy(0, 2.5); // half-cauchy

    /*
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    ph ~ gamma(0.001, 0.001);
    */

    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*x[n]);
    }  

  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;
  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string)   

saveName = 'test'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

nChains = 4
iterSteps = 10000
burnInSteps =1000

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
for ( parName in c('a','b','ph') ) {
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
chain_s = mcmcMat[,"ph"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
summaryInfo = rbind( summaryInfo , 
                    "ph" = summarizePost( chain_s ) )

show(summaryInfo)

###############################################################################

saveName = 'test_error'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

nChains = 4
iterSteps = 10000
burnInSteps =1000

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList_error , 
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
for ( parName in c('a','b','ph') ) {
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
chain_s = mcmcMat[,"ph"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
summaryInfo = rbind( summaryInfo , 
                    "ph" = summarizePost( chain_s ) )

show(summaryInfo)

# betareg
#lm_xerr <- betareg(y2~xerr)
#print(lm_xerr)

}

###############################################################################
# adding Uncertainty
###############################################################################

if ( flag_uncertainty ) {

dataList = list(
  x = xerr, 
  y = y2,
  Ndata = length(xerr),
  sx = 0.5
)

model_string5 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;
  vector<lower=0, upper=1>[Ndata] y ;
  real<lower=0> sx; // measurement error
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> ph; // sigma
  vector[Ndata] xtrue;
}
model {
    vector<lower=0,upper=1>[Ndata] mu;
    vector<lower=0>[Ndata] beta_a;
    vector<lower=0>[Ndata] beta_b;
    /*
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    ph ~ gamma(0.001, 0.001);
    */
    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5);
    ph ~ cauchy(0, 2.5); // half-cauchy
    xtrue ~ normal(0, 5);
    x ~ normal(xtrue, sx);
    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*xtrue[n]);
    }  
  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;
  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

model_string3 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;// x
  vector[Ndata] y ;// y 
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> s; // sigma
  vector[Ndata] xtrue; 
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    s ~ gamma(0.001, 0.001);
    x ~ normal(xtrue, 1.0);
    y ~ normal(a + b*xtrue, s);
}
"# close quote for modelString

model_string2 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ; // normalized IM
  int<lower=0> y[Ndata] ; // no. of head 
  int<lower=1> N[Ndata] ; // no. of flip
}
parameters {
  real a; // parameter of linear predictor
  real<lower=0> b; // a+bx
  //real<lower=0> s; // sigma of observed im
  vector[Ndata] xtrue; // true x
}
/*
transformed parameters { // Matt trick
  vector[Ndata] z;
  z <- a + b * xtrue;
}
*/
model {
    a ~ normal(0, 1.0e+1);
    b ~ normal(0, 1.0e+1);
    //s ~ lognormal(log(0.35), 0.05); // qlnorm(c(0.025,0.975),log(0.35),0.2)
    //s ~ gamma(8.0, 1.0);
    //z <- operator*(b, xtrue)
    //y ~ binomial_logit(N,z);
    for (n in 1:Ndata) {
      xtrue[n] ~ normal(0.0,5.0);
      x[n] ~ normal(xtrue[n],0.35);
      y[n] ~ binomial(N[n],Phi_approx(a + b*xtrue[n]));
    }
}
"# close quote for modelString

#      //
#      //print("s=", s, "xtrue=", xtrue[n],"x=", x[n]);

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string5)   

saveName = 'test_uncertainty'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

nChains = 4
iterSteps = 20000
burnInSteps =5000

# mod1 <- stan(model_code = model_string2, data = dataList,
#     chains = nChains, iter = iterSteps,
#     thin = 1, warmup = burnInSteps)

# posteriors <- extract(mod1)

  # Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = 1 )

# #print(stanFit, digits_summary=3, pars=c('a','b'), probs=c(0.025, 0.5, 0.975))
# #traceplot(stanFit, pars=c('a','b','s'))

codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                  function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in c('a','b','ph') ) {
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
chain_s = mcmcMat[,"ph"]
#chain_sx = mcmcMat[,"sx"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
summaryInfo = rbind( summaryInfo , 
                    "ph" = summarizePost( chain_s ) )
#summaryInfo = rbind( summaryInfo , 
#                    "sx" = summarizePost( chain_sx ) )

show(summaryInfo)

#if ( !is.null(saveName) ) {
#  write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
#}

if (marPlot) { 
  # marginal distribution
  openGraph(width=6,height=3)
  layout( matrix( 1:3 , nrow=1) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(sigma)) #, main=paste("a=",parameters$a) )
#  histInfo = plotPost( chain_sx , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(c)) #, main=paste("a=",parameters$a) )
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
  pairs( cbind( chain_a , chain_b, chain_s )[plotIdx,] ,
         labels=c( expression(a) , expression(b), expression(sigma) ) , 
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

  y = dataList$y
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
  xLim= c( min(x) - xLimMult*xRang, max(x) + xLimMult*xRang )
  yLim= c( min(y) - yLimMult*yRang, max(y) + yLimMult*yRang )
  plot( x, y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
        xlab='IM' , ylab='Probability' , cex.lab=1.5 ,
        main=NULL) #paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
        #cex.main=1.33  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=50
  xComb = seq(xLim[1],xLim[2],length=501)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    lines( xComb, 1/(1+exp(-(chain_a[i] + chain_b[i]*xComb))) , col="skyblue" )
  }

  # for (i in 1:nSlice) {
  #   yComb = seq(yHDI_mat[i,1],yHDI_mat[i,2],length=75)
  #   lines( mmi_list[i] - 0*yComb , yComb , col="red" , lwd=2 ) 
  # }

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"comPlot",sep=""), type=saveType)
  }
}

}

if ( flag_uncertainty2 ) {
###############################################################################
# adding Uncertainty
###############################################################################

dataList = list(
  x = xerr, 
  y = y2,
  Ndata = length(xerr)
)

model_string6 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;
  vector<lower=0, upper=1>[Ndata] y ;
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> ph; // sigma
  real<lower=0> sx; // sigma of x
  vector[Ndata] xtrue;
}
model {
    vector<lower=0,upper=1>[Ndata] mu;
    vector<lower=0>[Ndata] beta_a;
    vector<lower=0>[Ndata] beta_b;
    /*
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    ph ~ gamma(0.001, 0.001);
    */
    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5);
    ph ~ gamma(2.0, 4.0); //
    //ph ~ cauchy(0, 2.5); // half-cauchy
    //sx ~ cauchy(0, 2.5); // half-cauchy
    sx ~ gamma(2.0, 4.0); 
    xtrue ~ normal(0, 3);
    x ~ normal(xtrue, sx);
    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*xtrue[n]);
    }  
  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;
  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

model_string4 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;// x
  vector[Ndata] y ;// y 
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> s; // sigma
  vector[Ndata] xtrue; 
  real<lower=0> sx; // sigma x
}
model {
    a ~ normal(0, 1.0e+2);
    b ~ normal(0, 1.0e+2);
    sx ~ gamma(10.0, 5.0);
    x ~ normal(xtrue, sx);
    y ~ normal(a + b*xtrue, s);
}
"# close quote for modelString

model_string2 = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ; // normalized IM
  int<lower=0> y[Ndata] ; // no. of head 
  int<lower=1> N[Ndata] ; // no. of flip
}
parameters {
  real a; // parameter of linear predictor
  real<lower=0> b; // a+bx
  //real<lower=0> s; // sigma of observed im
  vector[Ndata] xtrue; // true x
}
/*
transformed parameters { // Matt trick
  vector[Ndata] z;
  z <- a + b * xtrue;
}
*/
model {
    a ~ normal(0, 1.0e+1);
    b ~ normal(0, 1.0e+1);
    //s ~ lognormal(log(0.35), 0.05); // qlnorm(c(0.025,0.975),log(0.35),0.2)
    //s ~ gamma(8.0, 1.0);
    //z <- operator*(b, xtrue)
    //y ~ binomial_logit(N,z);
    for (n in 1:Ndata) {
      xtrue[n] ~ normal(0.0,5.0);
      x[n] ~ normal(xtrue[n],0.35);
      y[n] ~ binomial(N[n],Phi_approx(a + b*xtrue[n]));
    }
}
"# close quote for modelString

#      //
#      //print("s=", s, "xtrue=", xtrue[n],"x=", x[n]);

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string6)   

saveName = 'test_uncertainty6'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

nChains = 4
iterSteps = 50000
burnInSteps =15000

# mod1 <- stan(model_code = model_string2, data = dataList,
#     chains = nChains, iter = iterSteps,
#     thin = 1, warmup = burnInSteps)

# posteriors <- extract(mod1)

  # Get MC sample of posterior:
stanFit <- sampling( object=stanDso , 
                         data = dataList , 
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps , 
                         warmup = burnInSteps , 
                         #init = initsList , # optional
                         thin = 1 )

# #print(stanFit, digits_summary=3, pars=c('a','b'), probs=c(0.025, 0.5, 0.975))
# #traceplot(stanFit, pars=c('a','b','s'))

codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                  function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in c('a','b','ph','sx') ) {
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
chain_s = mcmcMat[,"ph"]
chain_sx = mcmcMat[,"sx"]

# summaryInfo
summaryInfo = rbind( summaryInfo , 
                   "a" = summarizePost( chain_a ) )
summaryInfo = rbind( summaryInfo , 
                   "b" = summarizePost( chain_b ) )
summaryInfo = rbind( summaryInfo , 
                    "s" = summarizePost( chain_s ) )
summaryInfo = rbind( summaryInfo , 
                    "sx" = summarizePost( chain_sx ) )

show(summaryInfo)

#if ( !is.null(saveName) ) {
#  write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
#}

if (marPlot) { 
  # marginal distribution
  openGraph(width=6,height=3)
  layout( matrix( 1:4 , nrow=1) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(sigma)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_sx , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(sigma_x)) #, main=paste("a=",parameters$a) )
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
  pairs( cbind( chain_a , chain_b, chain_s, chain_sx )[plotIdx,] ,
         labels=c( expression(a) , expression(b), expression(sigma), expression(sigma_x) ) , 
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

  y = dataList$y
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
  xLim= c( min(x) - xLimMult*xRang, max(x) + xLimMult*xRang )
  yLim= c( min(y) - yLimMult*yRang, max(y) + yLimMult*yRang )
  plot( x, y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
        xlab='IM' , ylab='Probability' , cex.lab=1.5 ,
        main=NULL) #paste( "Data w. Post. Pred. & ",postPredHDImass*100,"% HDI" ,sep="") , 
        #cex.main=1.33  )
  # Superimpose a smattering of believable regression lines:
  nPredCurves=50
  xComb = seq(xLim[1],xLim[2],length=501)
  for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    lines( xComb, 1/(1+exp(-(chain_a[i] + chain_b[i]*xComb))) , col="skyblue" )
  }

  # for (i in 1:nSlice) {
  #   yComb = seq(yHDI_mat[i,1],yHDI_mat[i,2],length=75)
  #   lines( mmi_list[i] - 0*yComb , yComb , col="red" , lwd=2 ) 
  # }

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"comPlot",sep=""), type=saveType)
  }
}

}