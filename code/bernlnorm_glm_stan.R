# regression bern-log-normal
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
require('qmap')
require(ggplot2)

logit <- function(p) {
  log(p/(1-p))
}

inv_logit <- function(x) {
  1/(1+exp(-x))
}

# generate synthetic data
gen_synthetic_data <- function(a, b, s, c, d, npts){

  # MMI
  xmin = 4
  xmax = 10  
  
  set.seed(123)
  
  x <- runif(npts, min = xmin, max = xmax)
  y <- numeric(length=npts)

  #mu_y <- inv_logit()
  #mu_lny <- log(mu_y) - 0.5*s**2.0 #mu = log(mux)-0.5*sigma**2
  #mu_lny <- inv_logit(a+b*x)
  
  mu_lny <- a+b*x
  sigma_lny <- s #1/k^2
  p <- inv_logit(c+d*x)

  #y <- rlnorm(npts, mu_lny, sigma_lny)
  for (i in 1:npts) {
    y[i] <- rbernlnorm(1, 1-p[i], mu_lny[i], sigma_lny)
  }

  #y <- rbernlnorm(npts, p, mu_lny, sigma_lny)

  dataList = list(
    x = x, 
    y = y,
    Ndata = npts
  )
  return(dataList)
}

a <- -10.0
b <- 0.5
c <- 10.0
d <- -2.0
s <- 0.5
npts <- 100

dataList <- gen_synthetic_data(a, b, s, c, d, npts)

model_string = "
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
  real<lower=0> s; // sigma lognormal constant cov
}
model {
    a ~ cauchy(0, 2.5);
    b ~ cauchy(0, 2.5); 
    c ~ cauchy(0, 2.5); 
    d ~ cauchy(0, 2.5); 
    s ~ cauchy(0, 2.5);  # half cauchy

    for (n in 1:Ndata) {
    (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
    if (y[n] > 0)
       y[n] ~ lognormal(a+b*x[n], s);
    }    

    /*
    y ~ lognormal(a + b*x, s);
    */

}
"# close quote for modelString

stanDso <- stan_model( model_code = model_string)

nChains = 4
iterSteps = 10000
burnInSteps =2000

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                         data = dataList ,
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps ,
                         warmup = burnInSteps ,
                         #init = initsList , # optional
                         thin = 1 )

# # increase npts
# dataList1 <- gen_synthetic_data(a, b, s, c, d, npts=200)

# # Get MC sample of posterior:
# stanFit1 <- sampling( object=stanDso ,
#                          data = dataList1 ,
#                          #pars = parameters , # optional
#                          chains = nChains ,
#                          iter = iterSteps ,
#                          warmup = burnInSteps ,
#                          #init = initsList , # optional
#                          thin = 1 )

#traceplot(stanFit1, pars=param_names)

#print(stanFit, digits_summary=3, pars=param_names, probs=c(0.025, 0.5, 0.975))
#print(stanFit1, digits_summary=3, pars=param_names, probs=c(0.025, 0.5, 0.975))

saveName = NULL

codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                 function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

parameterNames = varnames(codaSamples) # get all parameter names
parameterNames <- parameterNames[-length(parameterNames)] # remove the last element

for ( parName in parameterNames ) {
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

mmi_list <- seq(4.0, 10.0, 0.5)
nSlice = length(mmi_list)
nsamples = 1

# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
#for ( chainIdx in 1:chainLength ) {

estimate_fat_rate_HDI <- function(nsamples) {

  nfat = nsamples * chainLength  
  fat_rate <- matrix(NA, nrow = nSlice, ncol = nsamples * chainLength  )
  fat_rate_HDI <- matrix(NA, nrow=nSlice, ncol=2)

  for ( j in 1:nSlice ) {

    mmi <- mmi_list[j]

    for (k in 1:chainLength) {
      
      p1 <- inv_logit(chain_c[k] + chain_d[k]*mmi)
      y <- rbernlnorm(nsamples, 1-p1, chain_a[k] + chain_b[k]*mmi, chain_s[k])   
      fat_rate[j, seq(1+nsamples*(k-1), nsamples*k)] = y

    }

    fat_rate_HDI[j,] <- HDIofMCMC(fat_rate[j,])  

  }

  # comparison with data
  fat_rate_HDI <- as.data.frame(fat_rate_HDI)
  names(fat_rate_HDI) <- c("L", "U")
  fat_rate_HDI$mmi_bin <- mmi_list

  return(fat_rate_HDI)
}

fat_rate_HDI1 <- estimate_fat_rate_HDI(1)
fat_rate_HDI5 <- estimate_fat_rate_HDI(5)
#fat_rate_HDI10 <- estimate_fat_rate_HDI(10)
fat_rate_HDI100 <- estimate_fat_rate_HDI(100)

fig <- ggplot() +
 geom_errorbar(data = fat_rate_HDI5, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi_bin), width=0.2, colour='red')+
 geom_errorbar(data = fat_rate_HDI100, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi_bin), width=0.2, colour='black')

comPlot = FALSE

if (comPlot) {

  dat <- data.frame(mmi=dataList$x, rat=dataList$y)

  fig6 <- ggplot(data = dat, aes(x=mmi, y=rat)) +
   geom_point() +
   geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi_bin), width=0.2, colour='red')

   nPredCurves=50
   xComb = seq(4,10,length=501)
   yComb <- numeric(length=501)
   PostPred <- data.frame(yComb)

  # Data with superimposed regression lines and noise distributions:
   for ( i in floor(seq(from=1,to=chainLength,length=nPredCurves)) ) {
    
    PostPred$yComb <- exp(chain_a[i] + chain_b[i]*xComb)
    PostPred$xComb <- xComb

    fig6 <- fig6 + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
  }    

  PostPred$sengara <- 10^(0.622*xComb-8.033)

  fig6 <- fig6 + geom_line(data = PostPred, aes(x=xComb, y=sengara), linetype=2, size=0.5)

  fig6 <- fig6 + theme_bw(base_size=10) +
  xlab("MMI") + ylab("Fatality rate") 

}