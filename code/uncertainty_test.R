# Modified from Stan-Ymet-Xmet-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
require(betareg)
source("DBDA2E-utilities.R")

flag_bayesian = 1
flag_uncertainty1 = 0
flag_uncertainty2 = 0

logit <- function(x) {
  return(log(x/(1-x))) #  
}

# inverse logit or logistic function
inv_logit <- function(x) {
  return(1/(1+exp(-x)))
}

# dummy data set 1
set.seed(999)
n.reps <- 1
n.repeated <- 10
n <- 50
gain <- 1.0
threshold <- 0.0
sx <- 0.5  # measurement error
ph <- 10.0

# true covariate values
xtrue <- runif(n, -3, 3)

# random subset to perform repeat covariate measurements
#which.repeated <- sample(n, n.repeated)
xerr <- xtrue + rnorm(n + (n.repeated * (n.reps - 1)), 0, sx)

mu <- inv_logit(gain*(xtrue-threshold))

# mu: mean
# ph: dispersion parameter 
# E(x) = mu, V(x) = mu(1-mu)/(1+ph)
beta_a <- mu*ph
beta_b <- (1-mu)*ph

# x_ <- seq(0,1, 0.001)
# for (i in 1:10) {
#   pdf_ <- dbeta(x_, beta_a[i], beta_b[i])
#   x11();plot(x_, pdf_)
# }

y <- rbeta(n, beta_a, beta_b)  # alpha=0, beta=1, sdy=1
y[y==1.0] = 1.0-1.0e-4

logit_y <- logit(y)

glm.y <- glm(logit_y~xtrue)

x.limits <- extendrange(xtrue)
x.grid <- seq(x.limits[1],x.limits[2],length.out=500)

glm.y.predictions <- predict(glm.y, newdata=data.frame(xtrue=x.grid), se.fit=TRUE)

x11()
plot(logit_y ~ xtrue)
points(glm.y.predictions$fit ~ x.grid)

###############################################################################
if (flag_bayesian) {

dataList = list(
  x = xtrue, 
  y = y,
  Ndata = length(xtrue)
)

# dataList_error = list(
#   x = xerr, 
#   y = y,
#   Ndata = length(xerr)
# )

model_string = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;// x
  vector<lower=0, upper=1>[Ndata] y ;// y 
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> ph; // dispersion parameter
}
model {
    vector<lower=0,upper=1>[Ndata] mu;
    vector<lower=0>[Ndata] beta_a;
    vector<lower=0>[Ndata] beta_b;

    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5);
    ph ~ cauchy(0, 2.5); // half-cauchy

    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*x[n]);
    }  

  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;

  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

model_string_normal = "
data {
  int<lower=1> Ndata ;
  vector[Ndata] x ;// x
  vector<lower=0, upper=1>[Ndata] y ;// y 
}
transformed data {
  vector[Ndata] logit_y;
  for (n in 1:Ndata) {
    logit_y[n] <- logit(y[n]);
  }
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real<lower=0> ph; // std
}
model {
    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5);
    ph ~ cauchy(0, 2.5); // half-cauchy

  logit_y ~ normal(a + b*x, ph);
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string_normal)   

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

# beta function with logit
# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string)   

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

comPlot = TRUE
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

}

# ph = 5.0

#           Mean     Median        Mode     ESS HDImass     HDIlow  HDIhigh
# a  -0.02514297 -0.0271015 -0.02895614 21760.8    0.95 -0.6354131 0.562617
# b   1.68139981  1.6821222  1.69598805 22525.9    0.95  1.3229607 2.024515
# ph  2.18735140  2.1679294  2.13866172 23657.9    0.95  1.7789941 2.641201
#    CompVal PcntGtCompVal ROPElow ROPEhigh PcntLtROPE PcntInROPE PcntGtROPE
# a       NA            NA      NA       NA         NA         NA         NA
# b       NA            NA      NA       NA         NA         NA         NA
# ph      NA            NA      NA       NA         NA         NA         NA

#           Mean      Median        Mode     ESS HDImass     HDIlow   HDIhigh
# a  -0.04904111 -0.04896401 -0.04733723 16128.1    0.95 -0.3167861 0.2165098
# b   1.00307827  1.00360975  1.01043296 15333.9    0.95  0.8018783 1.2012700
# ph  4.38108575  4.30770566  4.21962905 16198.6    0.95  2.6617481 6.1574643


# ph = 1.0

#         Mean    Median      Mode     ESS HDImass    HDIlow    HDIhigh CompVal
# a  -1.363662 -1.339287 -1.278378 24123.8    0.95 -3.542451  0.9035787      NA
# b   5.658388  5.666626  5.736183 22195.8    0.95  4.220801  7.0628186      NA
# ph  8.718020  8.638769  8.462283 21729.5    0.95  6.996567 10.4805991      NA

#          Mean     Median       Mode     ESS HDImass     HDIlow   HDIhigh
# a  -0.1235962 -0.1241187 -0.1414433 14243.9    0.95 -0.4383444 0.1994826
# b   1.1661290  1.1662286  1.1696557 13138.3    0.95  0.9273873 1.4137818
# ph  1.2935339  1.2718709  1.2415514 12908.1    0.95  0.7696670 1.8291748


###############################################################################

# saveName = 'test_error'; saveType='eps'
# marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

# # Get MC sample of posterior:
# stanFit <- sampling( object=stanDso , 
#                          data = dataList_error , 
#                          #pars = parameters , # optional
#                          chains = nChains ,
#                          iter = iterSteps , 
#                          warmup = burnInSteps , 
#                          #init = initsList , # optional
#                          thin = 1 )

# #print(stanFit, digits_summary=3, pars=c('a','b'), probs=c(0.025, 0.5, 0.975))
# #traceplot(stanFit, pars=c('a','b'))

# codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
#                                  function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

# #parameterNames = varnames(codaSamples) # get all parameter names
# for ( parName in c('a','b','ph') ) {
# diagMCMC( codaObject=codaSamples , parName=parName , 
#           saveName= saveName , saveType= 'png')
# }

# mcmcMat = as.matrix(codaSamples,chains=TRUE)

# #if ( !is.null(saveName) ) {
# #  saveRDS(mcmcMat, file=paste(saveName,'mcmcMat.RDS',sep=''))
# #}

# chainLength = nrow( mcmcMat )
# summaryInfo = NULL

# chain_a = mcmcMat[,"a"]
# chain_b = mcmcMat[,"b"]
# chain_s = mcmcMat[,"ph"]

# # summaryInfo
# summaryInfo = rbind( summaryInfo , 
#                    "a" = summarizePost( chain_a ) )
# summaryInfo = rbind( summaryInfo , 
#                    "b" = summarizePost( chain_b ) )
# summaryInfo = rbind( summaryInfo , 
#                     "ph" = summarizePost( chain_s ) )

# show(summaryInfo)

}

###############################################################################
# adding Uncertainty
###############################################################################

if ( flag_uncertainty1 ) {

dataList = list(
  x = xerr, 
  y = y,
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
    //xtrue ~ normal(0, 3);
    x ~ normal(xtrue, sx);
    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*xtrue[n]);
    }  
  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;
  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string5)   

saveName = 'test_uncertainty1'; saveType='eps'
marPlot = TRUE; pairsPlot = TRUE; comPlot= TRUE

nChains = 4
iterSteps = 20000
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

}

if ( flag_uncertainty2 ) {
###############################################################################
# adding Uncertainty
###############################################################################

dataList = list(
  x = xerr, 
  y = y,
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
    //xtrue ~ normal(0, 3);
    x ~ normal(xtrue, sx);
    for (n in 1:Ndata) {
      mu[n] <- inv_logit(a + b*xtrue[n]);
    }  
  beta_a <- mu * ph;
  beta_b <- (1-mu) * ph;
  y ~ beta(beta_a, beta_b);
}
"# close quote for modelString

# Translate to C++ and compile to DSO:
stanDso <- stan_model( model_code = model_string6)   

saveName = 'test_uncertainty2'; saveType='eps'
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

}