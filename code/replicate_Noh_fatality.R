graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="replicate_Noh_fatality" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)

calculate_gamma_parameters <- function(mu0, sd0){
    tau = 1/sd0^2
    r <- mu0*tau
    s <- mu0^2*tau
    #newList <- list("float" = c(a0,b0), "names" = bar)    
return(list("r"=r,"s"=s))
}

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}

#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
# using zero trick
data {
   C <- 1000000    # this just has to be large enough to ensure all phi[i]'s < 1
                           # Note that JAGS will NOT WARN you if C is too small, but will sample from a truncated distribution
   for (j in 1:N) {
      ones[j] <- 1
      zeros[j] <- 0
   }
}
model {
      for (j in 1:N) {
         phi[j] <- dexp(lambda)T(a0,1) 
         ones[j] ~ dbern(phi[j])
         #phi[j] <-  log(sigma) + 0.5 * pow((x[j] - mu) / sigma, 2) + C
         #zeros[j] ~ dpois(phi[j])
   }
   #mu ~ dunif(-5.0,5.0)
   #sigma ~ dunif(0.0001, 5.0)
    lambda ~ dunif(0.0001,20.00)
    #p0 ~ dunif(0.0001,0.99)
    a0 ~ dunif(0.0001,0.4)
}
" # close quote for modelstring


#
# model {
#    for( i in 1:N ) {
#       #y_1[i] ~ dbern(p0)
#       #y_2[i] ~ dexp(lambda)T(a0,1)
#       #y[i] ~ dnorm((1-y_1[i])*y_2[i], 1000000000);
#       #y[i] ~ dsum(y_3[i],0)
#       y[i] ~ dexp(lambda)T(0.1,1)
#     }
#     #lambda ~ dgamma( 74.07 , 11.11 )
#     lambda ~ dunif(0.0001,20.00)
#     #p0 ~ dunif(0.0001,0.99)
#     #a0 ~ dunif(0.0001,0.4)
# }

      #y_branch[1] <- dbern(p0)
      #y_branch[2] <- dexp(lambda)
      #if_branch <- ifelse(y_branch[i,1]==0,1,2)
      #y[i] ~  y_branch[i,if_branch]
      #y[i] ~  y_branch[2]
      #y[i] ~ dgamma(lambda)

# model {
#    for( i in 1 : N ) {
#       # y[i] ~ dexp(lambda)
#    }
#     lambda ~ dgamma( 74.07 , 11.11 )
# }

# Z_branch[i, 1] ~ dcat(theta-gl[D[i], 1:K-gl])
# W_branch[i, 1] ~ dcat(phi-gl[z[i], 1:V])

# # else branch
# Z_branch[i, 2] ~ dcat(theta-loc[D[i], 1:K-loc])
# W_branch[i, 2] ~ dcat(phi-loc[z[i], 1:V])

# # decision here
# if_branch[i] <- 1 + step(-(X[i] - 0.5)) # 1 for "then" branch, 2 for "else" branch
# Z[i] ~ Z_branch[i, if_branch[i]]
# W[i] ~ W_branch[i, if_branch[i]]

#writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# synthetic data
#y = readRDS('syn_fatal.RDS')

npts = 100

#y <- rexp(npts,rate=1/0.107)
y <- rtrunc(npts,"exp",a=0.1,b=1,rate=1/0.107)
#q <- quantile(y,probs=0.3)
#y[y<q] = 0.0

dataList = list(
    y = y,
    N = length(y)
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

#initsList = list( p0 = 0.5, lambda = 1/mean(y)) 
#initsList = list( a0 = 0.5, lambda = 1/mean(y)) 
initsList = list( lambda = 1/mean(y)) 
#------------------------------------------------------------------------------
# RUN THE CHAINS

#parameters = c("lambda","p0")  # The parameter(s) to be monitored.
#parameters = c("a0","p0")  # The parameter(s) to be monitored.
parameters = c("lambda")  # The parameter(s) to be monitored.
#parameters = c("p0")  # The parameter(s) to be monitored.
adaptSteps = 5000              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( textConnection(modelstring) , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters, 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]] , ask=FALSE )
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples ) 
  show( effectiveChainLength )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

z0 = mcmcChain[, "lambda" ]
#z1 = mcmcChain[, "p0" ]
#z1 = mcmcChain[, "a0" ]

# specify preferred graph-file type:
graphType="eps" # or "jpg" or whatever you want

# Marginal distributions of mu and sigma:
openGraph(width=7,height=3.5)
layout(matrix(1:2,nrow=1))
plotPost( z0 , xlab=bquote(lambda))#,xlim=c(0.08,0.18) )
#plotPost( z1 , xlab=bquote(p[0]))#,xlim=c(0.08,0.18) )
#saveGraph(file=paste(fileNameRoot,"PostMuSigmaMarg_n",npts,sep=""),type=graphType)

# yp <- vector("numeric",length=length(z0))
# for (i in 1:length(z0)) {
#   #x_[i] <- (1-rbinom(1,1,z1[i]))*rexp(1,z0[i])
#   yp[i] <- (1-rbinom(1,1,z1[i]))*rtrunc(1,"exp",b=1,rate=z0[i])
# }

# X11()
# y.ordered = sort(y)
# yp.ordered = sort(yp)
# # plot the possible values of probability (0 to 1) against the ordered ozone data (sample quantiles of ozone)
# # notice the option type = 's' for plotting the step functions
# plot(y.ordered, (1:length(y))/length(y), type = 's', ylim = c(0, 1),col='blue')
# points(yp.ordered, (1:length(yp))/length(yp), type = 's', col='red',t='l')



# # Joint distribution of mu and sigma:
# nPts = length(muSample) ; nPtsForDisplay = min( nPts , 2000 )
# thinIdx = seq( 1 , nPts , nPts / nPtsForDisplay )
# openGraph(width=5,height=5)
# plot( muSample[thinIdx] , sigmaSample[thinIdx] , col="skyblue" ,
#       xlab=bquote(mu) , ylab=bquote(sigma) , cex.lab=1.5 )
# points( mean(muSample) , mean(sigmaSample) , pch="+" , cex=2 )
# text( mean(muSample) , mean(sigmaSample) ,
#       bquote( .(round(mean(muSample),1)) *"  "* .(round(mean(sigmaSample),1)) ),
#       adj=c(.5,-0.5) )
# saveGraph(file=paste(fileNameRoot,"PostMuSigmaJoint","_",modelStr,sep=""),type=graphType)

# # Posterior predictive check:
# openGraph(width=5,height=4)
# histInfo = hist( dataList$y , xlab="y (data)" , 
#                  main="Data with Posterior Pred. Distrib." , 
#                  breaks=20 , col="grey" , border="white" , prob=TRUE )
# yLim = range( histInfo$breaks )
# yComb = seq( yLim[1] , yLim[2] , length=501 )
# chainLength = length(muSample)
# chainIdxVec = floor( seq(1,chainLength,length=20) )

# if (modelStr == 'gamma') {
#   for ( i in chainIdxVec ) {
#     r0 = muSample[i]/sigmaSample[i]^2
#     s0 = muSample[i]^2/sigmaSample[i]^2
#     lines( yComb , dgamma( yComb , s0 , r0 ) ,
#            col="skyblue" ) }
# } else if (modelStr == 'beta') {
#   for ( i in chainIdxVec ) {
#     a0 = ((1-muSample[i])/sigmaSample[i]^2-1/muSample[i])*muSample[i]^2
#     b0 = a0*(1/muSample[i]-1)
#     lines( yComb , dbeta( yComb , a0 , b0 ) ,
#            col="skyblue" )}
# } else if (modelStr == 'tnorm') {
#   for ( i in chainIdxVec ) {
#     lines( yComb , dnorm( yComb , mean=muSample[i] , sd=sigmaSample[i] ) ,
#            col="skyblue" )}
# }
# saveGraph(file=paste(fileNameRoot,"PostPredict","_",modelStr,sep=""),type=graphType)



#------------------------------------------------------------------------------