# extension of bernlnorm_stan.R
#-------------------------------------------------------------------------------
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)

set.seed(123)

# environment
if (Sys.info()['sysname']=="Darwin") {
  comPath = '/Users/hyeuk/Projects/fatality/'
  codePath = paste(comPath,'code/',sep="")
  plotPath = paste(comPath,'plot/',sep="")
  dataPath = paste(comPath,'data/',sep="")
  saveType='eps'
} else {
  comPath = '/home/547/hxr547/Projects/fatality/'
  codePath = paste(comPath,'code/',sep="")
  plotPath = paste(comPath,'plot/',sep="")
  dataPath = paste(comPath,'data/',sep="")
  saveType='pdf'
}

# load utilities
source(paste(codePath,"utilities_graphic.R",sep=""))
source(paste(codePath,"utilities.R",sep=""))

# load stan models
source(paste(codePath,"stan_models.R",sep=""))

data_str = 'case1'
#link_str = 'log_add'
#link_str = 'binomial'
link_str = 'hi'
saveName = paste(data_str,'_berngamma_',link_str,sep="")
marPlot = TRUE;

# Translate to C++ and compile to DSO:
if (link_str=='log') {
  stanDso <- stan_model( model_code = model_berngamma_log)
} else if (link_str=='log_add') {
  stanDso <- stan_model( model_code = model_berngamma_log_add)
} else if (link_str=='binomial') {
  stanDso <- stan_model( model_code = model_binomialgamma_log_add)
} else if (link_str=='logit') {
  stanDso <- stan_model( model_code = model_berngamma_logit)
} else if (link_str=='hi') {
  stanDso <- stan_model( model_code = model_bernlnorm_add_hi)
}

# read fatality data
if (data_str=='wald') {
  dat <- read.csv(paste(dataPath,'DATA_WALD_COR_ROUND_12_Feb_2013.csv',sep=""), header=0)
} else if (data_str=='worden') {
  dat <- read.csv(paste(dataPath,'DATA_WORDEN_COR_ROUND_12_Feb_2013.csv',sep=""), header=0)
} else if (data_str=='case1') {
  dat <- read.csv(paste(dataPath,'case1.csv',sep=""), header=0)
} else if (data_str=='case2') {
  dat <- read.csv(paste(dataPath,'case2.csv',sep=""), header=0)
} else if (data_str=='case3') {
  dat <- read.csv(paste(dataPath,'case3.csv',sep=""), header=0)
}
names(dat) <- c("pop","fat","mmi","mmi_bin","id")
dat$rat <- dat$fat/dat$pop

# read expo cat data
dummy <- read_expo_cat(paste(dataPath,'EXPO_CAT_2007_12.csv', sep=""))
pop <- dummy[[1]]
expo_cat <- dummy[[2]]
mmi_list <- dummy[[3]]

dat$mmi_id <- sapply(dat$mmi,function(x) which(mmi_list==x))

if (link_str=='binomial') {
    dataList = list(
      x = dat$mmi,
      y = dat$rat,
      pop = dat$pop,
      fat = dat$fat,
      Ndata = length(dat$mmi),
      Nnew = length(mmi_list),
      xnew = mmi_list
    )
} else if (link_str=='hi') {
    dataList = list(
      x = dat$mmi,
      y = dat$rat,
      Ndata = length(dat$mmi),
      N_mmi = length(mmi_list),
      mmi_bin = dat$mmi_id,
      Nnew = length(mmi_list),
      xnew = mmi_list
    )
} else {
    dataList = list(
      x = dat$mmi,
      y = dat$rat,
      Ndata = length(dat$mmi),
      Nnew = length(mmi_list),
      xnew = mmi_list
    )
}

nChains = 4
burnInSteps = 1000
iterSteps = burnInSteps + 8000

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                     data = dataList ,
                     #pars = parameters , # optional
                     chains = nChains ,
                     iter = iterSteps ,
                     warmup = burnInSteps ,
                     #init = initsList , # optional
                     thin = 1 )

#stanSamples <- extract(stanFit)
#paramNames <- names(stanSamples)
#paramNames <- paramNames[-length(paramNames)] # remove the last

# diagplot #2
codaSamples = mcmc.list( lapply( 1:ncol(stanFit) ,
                                 function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
#paramNames = varnames(codaSamples)
#paramNames <- paramNames[-length(paramNames)] # remove the last
paramNames <- c('a','b','c','d','s')

for ( parName in paramNames ) {
  diagMCMC( codaObject=codaSamples, parName=parName ,
          saveName=paste(plotPath,saveName,'_',sep=""), saveType=saveType)
}

# traceplot
if (Sys.info()['sysname']=="Linux") {
  png(file=paste(plotPath,saveName,"_traceplot.png",sep=""), height=5, width=7, units="in", res=72)
  rstan::traceplot(stanFit, pars=paramNames, nrow=2, ncol=3)
  dev.off()
} else {
  openGraph(width=7, height=5)
  rstan::traceplot(stanFit, pars=paramNames, nrow=2, ncol=3)
  saveGraph( file=paste(plotPath,saveName,"_traceplot",sep=""), type=saveType)
}

mcmcMat = as.matrix(codaSamples,chains=TRUE)
chainLength = nrow( mcmcMat )

summaryInfo = NULL
chain <- data.frame(dummy = numeric(chainLength))
for ( parName in paramNames ) {
  chain[[parName]] <- mcmcMat[,parName]
  summaryInfo = rbind( summaryInfo, summarizePost( chain[[parName]] ))
}
chain$dummy <- NULL
row.names(summaryInfo) <- paramNames

if (marPlot) {

  if (Sys.info()['sysname']=="Linux") {
    png(file=paste(plotPath,saveName,"_PostMarg.png",sep=""), height=3*2, width=3*3, units="in", res=600)
  } else {
    openGraph(width=3*3, height=3*2)
  }

  # layout( matrix( 1:6 , nrow=2, byrow=TRUE) )

  # for ( parName in paramNames ) {
  #   histInfo = plotPost( chain[[parName]] , cex.lab = 1.75 , showCurve=FALSE, xlab=parName) #, main=paste("a=",parameters$a) )
  # }

  layout( matrix( 1:6 , nrow=2) )

  histInfo = plotPost( chain[['a']] , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain[['c']] , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(c)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain[['b']] , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain[['d']] , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(d)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain[['s']] , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(alpha)) #, main=paste("a=",parameters$a) )

  # saveGraph( file=paste(plotpath,"figure5.eps",sep=""), type=saveType)

  if ( !is.null(saveName) ) {
    if (Sys.info()['sysname']=="Linux") {
      dev.off()
    } else {
      saveGraph( file=paste(plotPath, saveName,"_PostMarg",sep=""), type=saveType)
    }
  }
}

# extract predicted fatality rate by MMI
sel_col <- vector("character", length(mmi_list))
for (i in seq(1, length(mmi_list))) {
  sel_col[[i]] <- paste('ynew[', i, ']', sep='')
}
fat_rate_by_mmi <- mcmcMat[, sel_col]

#nsamples <- 5
#if (link_str=='log') {
#  temp <- estimate_fat_rate_HDI_berngamma_log(nsamples)
#} else if (link_str=='logit') {
#  temp <- estimate_fat_rate_HDI_berngamma_logit(nsamples)
#}
#fat_rate <- temp[[1]]
#fat_rate_HDI <- temp[[2]]

# remove any fatality rate >= 1.0
#ff <- which(fat_rate >= 1.0, arr.ind=TRUE)
#if (nrow(ff) > 0) {
#  fat_rate <- fat_rate[,-ff[,2]]
#}

fat_by_event <- pop %*% t(fat_rate_by_mmi) # nevents x (nsamples*chainLength)

# HDI of estimated fatality
dummy = compute_HDI_prob(fat_by_event)
fatHDI = dummy[[1]]
prob_mag = dummy[[2]]

count_order <- count_order_match(prob_mag, expo_cat)
show(count_order)

if ( !is.null(saveName) ) {
  saveRDS(mcmcMat, file=paste(dataPath,saveName,'_mcmcMat.RDS',sep=""))
  write.csv(summaryInfo , file=paste(dataPath,saveName,"_SummaryInfo.csv",sep=""))
  #saveRDS(fat_rate, file=paste(dataPath,saveName,'_fat_rate.RDS',sep=''))

  #df.fat_rate <- as.data.frame(t(fat_rate))
  #names(df.fat_rate) <- mmi_list
  # remove .5 for inasafe
  #inasafe_mmi_list <- seq(4.0, 10.0)
  #df.fat_rate_inasafe <- df.fat_rate[as.character(inasafe_mmi_list)]
  #write.csv(df.fat_rate_inasafe, file=paste(dataPath,saveName,'_fat_rate_inasafe.csv',sep=""), row.names=FALSE)

  #saveRDS(fat_rate_HDI, file=paste(dataPath,saveName,'_fat_rate_HDI.RDS',sep=''))
  write.csv(fat_rate_by_mmi, file=paste(dataPath,saveName,'_fat_rate_by_mmi.csv',sep=""), row.names=FALSE)
  saveRDS(fat_by_event, paste(dataPath,saveName, '_fat_by_event.RDS',sep=""))
  saveRDS(fatHDI, paste(dataPath, saveName, '_fatHDI.RDS',sep=""))
  saveRDS(prob_mag, paste(dataPath, saveName, '_prob_mag.RDS',sep=""))
}

