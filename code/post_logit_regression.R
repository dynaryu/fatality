# postprocessing of logit_regression_raijin.R

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
# CTRL + L clear console
#require(rstan)
#source("DBDA2E-utilities.R")
#require('qmap')
library(MASS)
require(ggplot2)
library(scales)
library(reshape2)

# (0,1) -> (-inf, inf)
logit <- function(p) {
  log(p/(1-p))
}

# [-inf, inf] -> (0, 1)
inv_logit <- function(x) {
  1/(1+exp(-x))
}

count_order_match2 <- function(prob_mag, expo_cat) {

  nevents <- length(expo_cat$mag_obs)
  match_idx <- numeric(length=nevents)

  for (i in 1:nevents) {
    match_idx[i] <- sum(prob_mag[i,] >= prob_mag[i,][expo_cat$mag_obs[i]],na.rm=TRUE)
  }

  mytable <- table(match_idx, expo_cat$mag_obs)
  return(mytable)
}

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
  #eqName <- sel_dat$eqName

  # 0-10(1), 10-100(2), 100-1000(3), 1000-10000(4)
  mag_observed <- floor(log10(observed)) + 1
  mag_observed[is.infinite(mag_observed)] = 1 # 0 range

  #mag_observed <- floor(log10(observed)) + 2
  #mag_observed[is.infinite(mag_observed)] = 1 # 0-1 range

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

# Model based on case1
#mcmcMat_file = 'wald_berngamma_n15kmcmcMat.RDS'
#mcmcMat_file = 'case3_berngamma_logit_fat_rate.RDS'
marPlot = FALSE #TRUE
#saveName ='wald_berngamma_n15k'
saveName = 'case1_berngamma'
saveType = 'eps'

dummy <- read_expo_cat(fileName = '/Users/hyeuk/Projects/fatality/data/EXPO_CAT_2007_12.csv')
pop <- dummy[[1]]
expo_cat <- dummy[[2]]

datapath = '/Users/hyeuk/Projects/fatality/data/'
plotpath = '/Users/hyeuk/Projects/fatality/plot/'
mcmcMat_file = paste(datapath, 'case1_berngamma_log_add_mcmcMat.RDS', sep="")

mcmcMat <- readRDS(mcmcMat_file)
# chainLength = nrow( mcmcMat )
# summaryInfo = NULL
#
chain_a = mcmcMat[,"a"]
chain_b = mcmcMat[,"b"]
chain_c = mcmcMat[,"c"]
chain_d = mcmcMat[,"d"]
chain_s = mcmcMat[,"s"]


# figure 5
if (marPlot) {
  # marginal distribution
  #openGraph(width=9,height=5)
  dev.new(width=9, height=5)
  layout( matrix( 1:6 , nrow=2) )

  histInfo = plotPost( chain_a , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(a)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_c , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(c)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_b , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(b)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_d , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(d)) #, main=paste("a=",parameters$a) )
  histInfo = plotPost( chain_s , cex.lab = 1.75 , showCurve=FALSE, xlab=bquote(alpha)) #, main=paste("a=",parameters$a) )

  saveGraph( file=paste(plotpath,"figure5_marginal.eps",sep=""), type=saveType)
}

###############################################################################
# figure 6
# comparison with data
#dat <- readRDS('DATA_WALD_COR_ROUND_12_Feb_2013.RDS')
#fat_rate_HDI <- readRDS('wald_berngamma_n15kfat_rate_HDI.RDS')
dat <- read.csv(paste(datapath, 'case1.csv', sep=""), header = FALSE)
names(dat) <- c('pop', 'fat', 'mmi', 'mmi_round', 'id')
dat$rat <- dat$fat/dat$pop

fat_rate_by_mmi <- read.csv(paste(datapath, 'case3_berngamma_log_add_fat_rate_by_mmi.csv', sep=""))
mmi_list = seq(4, 10, 0.5)

# zero
prob_zero_by_mmi <- apply(fat_rate_by_mmi==0, FUN=sum, MARGIN = 2)/dim(fat_rate_by_mmi)[1]

# fat_rate
fat_rate_HDI <- matrix(NA, ncol=3, nrow=dim(fat_rate_by_mmi)[2])

for (i in 1:length(mmi_list)) {
  x <- fat_rate_by_mmi[[paste('ynew.',i,'.',sep="")]]
  #x <- x[x>0]
  dummy <- HDIofMCMC(x)
  fat_rate_HDI[[i, 1]] <- mmi_list[i]
  fat_rate_HDI[[i, 2]] <- dummy[1]
  fat_rate_HDI[[i, 3]] <- dummy[2]
}

#fat_rate_HDI <- readRDS(paste(datapath, '', sep=""))
fat_rate_HDI <- as.data.frame(fat_rate_HDI)
names(fat_rate_HDI) <- c("mmi","L","U")

# mmi_list = seq(4,10,0.5)
# fat_rate_HDI$mmi_bin <- mmi_list
#######
# logscale
######
fig6 <- ggplot(data = dat, aes(x=mmi, y=rat)) +
  geom_point() +
geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.1, colour='red') +
scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^0)) +
scale_x_continuous(limits=c(3.5,10.5), breaks=seq(4.0, 10.0, 1.0))
#geom_point(data = dat, aes(x=mmi, y=rat)) +
# geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.2, colour='red')
# scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
 # geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.2, colour='red')
# geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L), width=0.2, colour='red')

 #nPredCurves=50
 yComb <- numeric(length=length(mmi_list))
 PostPred <- data.frame(yComb)
 PostPred$xComb <- mmi_list

 #idx <- sample(seq(1, dim(fat_rate_by_mmi)[1]), 50)
 # PostPred <- rbind(as.double(mmi_list), as.double(fat_rate_by_mmi[1,]))
 # PostPred <- as.data.frame(t(PostPred))
 # fig6 <- fig6 + geom_line(data = PostPred, aes(x=V1, y=V2), colour='skyblue', size=0.1)

# # Data with superimposed regression lines and noise distributions:
 for ( i in idx ) {

  PostPred$yComb <- exp(mcmcMat[i,"a"] + mcmcMat[i,"b"]*mmi_list)
  #PostPred$xComb <- xComb

  #PostPred$yComb <- as.double(fat_rate_by_mmi[i, ])
  #PostPred$xComb <- mmi_list
  fig6 <- fig6 + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
  # fig6 <- fig6 + geom_point(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)

    }

PostPred$sengara <- 10^(0.622*PostPred$xComb-8.033)

fig6 <- fig6 + geom_line(data = PostPred, aes(x=xComb, y=sengara), linetype=2, size=0.5)

fig6 <- fig6 + theme_bw(base_size=10) +
xlab("MMI") + ylab("Fatality rate")

ggsave(paste(plotpath,'fig6_new',saveName, '.', saveType, sep= ''),  width = 16, height = 8, unit="cm",  fig6)


#######
# figure6 in linear scale
#######
fig6b <- ggplot(data = dat, aes(x=mmi, y=rat)) +
  geom_point() +
  geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.2, colour='red') +
  scale_x_continuous(limits=c(4.0,10.0), breaks=seq(4.0, 10.0, 1.0))

  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1))
#geom_point(data = dat, aes(x=mmi, y=rat)) +
# geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.2, colour='red')
# scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
# geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L, y=0.5*(L+U), x = mmi), width=0.2, colour='red')
# geom_errorbar(data = fat_rate_HDI, aes(ymax = U, ymin = L), width=0.2, colour='red')

#nPredCurves=50
yComb <- numeric(length=length(mmi_list))
PostPred <- data.frame(yComb)
PostPred$xComb <- mmi_list

#idx <- sample(seq(1, dim(fat_rate_by_mmi)[1]), 50)
# PostPred <- rbind(as.double(mmi_list), as.double(fat_rate_by_mmi[1,]))
# PostPred <- as.data.frame(t(PostPred))
# fig6 <- fig6 + geom_line(data = PostPred, aes(x=V1, y=V2), colour='skyblue', size=0.1)

# # Data with superimposed regression lines and noise distributions:
for ( i in idx ) {

  PostPred$yComb <- exp(mcmcMat[i,"a"] + mcmcMat[i,"b"]*mmi_list)
  #PostPred$xComb <- xComb

  #PostPred$yComb <- as.double(fat_rate_by_mmi[i, ])
  #PostPred$xComb <- mmi_list
  fig6b <- fig6b + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
  # fig6 <- fig6 + geom_point(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)

}

PostPred$sengara <- 10^(0.622*PostPred$xComb-8.033)

fig6b <- fig6b + geom_line(data = PostPred, aes(x=xComb, y=sengara), linetype=2, size=0.5)

fig6b <- fig6b + theme_bw(base_size=10) +
  xlab("MMI") + ylab("Fatality rate")

ggsave(paste(plotpath,'fig6b_new',saveName, '.', saveType, sep= ''),  width = 16, height = 8, unit="cm",  fig6b)


# plot mu, and p, beta

# # Data with superimposed regression lines and noise distributions:
yComb <- numeric(length=length(mmi_list))
PostPred <- data.frame(yComb)
PostPred$xComb <- mmi_list

########
# mean
########

fig7 <- ggplot()
for ( i in idx ) {
  PostPred$yComb <- exp(mcmcMat[i, "a"] + mcmcMat[i, "b"]*PostPred$xComb)
  fig7 <- fig7 + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
}

#####
# sigma
#####

fig9 <- ggplot()
for ( i in idx ) {
  PostPred$yComb <- exp(mcmcMat[i, "a"] + mcmcMat[i, "b"]*PostPred$xComb)/sqrt(mcmcMat[i, "s"])
  fig9 <- fig9 + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
}

#####
# cov
#####

fig10 <- ggplot()
for ( i in idx ) {
  PostPred$yComb <- 1/sqrt(mcmcMat[i, "s"])
  fig10 <- fig10 + geom_line(data = PostPred, aes(x=xComb, y=yComb), colour='skyblue', size=0.1)
}


#######
# zero prob
#######
yComb <- numeric(length=length(mmi_list))
PostPred <- data.frame(yComb)
PostPred$xComb <- mmi_list
PostPred$yComb <- prob_zero_by_mmi
fig8 <- ggplot(data=PostPred, aes(x=xComb, y=yComb)) +
  geom_line() +
  scale_x_continuous(limits=c(4.0,10.0), breaks=seq(4.0, 10.0, 1.0)) +
  xlab("MMI") +
  ylab("Probability of zero fatality rate") +
  theme_bw(base_size=10)
ggsave(paste(plotpath,'figure6b_zero_prob.eps',sep=""),  width = 8, height = 8, unit="cm",  fig8)

####
# fig 7a

plot_prob_one_event <- function( prob_mag, k, expo_cat, saveName, saveType ) {

  # names_str <- c('0-10^0','10^0-10^1','10^1-10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6','10^7-10^8')
  names_str <- c('0-10^1','10^1-10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6','10^7-10^8')

  label_str <- c(expression(0-10^1), expression(10^1-10^2),
    expression(10^2-10^3), expression(10^3-10^4), expression(10^4-10^5),
    expression(10^5-10^6), expression(10^6-10^7))

  prob_mag_event <- prob_mag[k,]

  nmax <- length(prob_mag_event)

  event1 <- matrix(NA, ncol=nmax, nrow=1)
  event1[1,] <- prob_mag_event
  event1 <-data.frame(event1)

  names(event1) <- names_str[1:nmax]
  event1$eqID <- expo_cat$eqID[k]
  event1.long<-melt(event1,id.vars=c("eqID"))

  event1_ <- matrix(NA, ncol=3, nrow=1)
  event1_[1,1] <- expo_cat$obs[k]
  event1_[1,2] <- expo_cat$mag_obs[k]
  event1_[1,3] <- 0*prob_mag_event[event1_[1,2]]
  event1_ <- data.frame(event1_)
  names(event1_) <- c('obs','x','y')

  myPlot <- ggplot(event1.long,aes(variable,value))+
   geom_bar(stat="identity", width=0.5, fill='grey')+
   coord_flip() +
   scale_x_discrete(labels= label_str) +
   geom_point(data=event1_,aes(x,y), col='red', size=4) +
   #facet_wrap(~eqID,ncol=1) +
   ylab("Probability") +
   xlab("Fatality ranges") +
   theme_bw(base_size=10)

  if ( !is.null(saveName) ) {

    ggsave(paste(plotpath, 'fig7b_', saveName, '.', saveType, sep= ''),  width = 8, height = 8, unit="cm",  myPlot)

  }
  return(myPlot)
}

# figure 7b
# prob_mag <- readRDS('wald_berngamma_n15kprob_mag.RDS')
prob_mag <- readRDS(paste(datapath, 'case3_berngamma_log_add_prob_mag.RDS', sep=""))

# one event
# histogram of one event (200006041628, Bengkulu, Indonesia) INDEX: 42
k = 42
saveName <- 'case3'
myPlot3 <- plot_prob_one_event(prob_mag, k, expo_cat, saveName, saveType)

# figure 7a
#fat <- readRDS('wald_berngamma_n15kfat_by_event_42.RDS')
#fatHDI <- readRDS('wald_berngamma_n15kfatHDI.RDS')
fat <- readRDS(paste(datapath, 'case3_berngamma_log_add_fat_by_event.RDS', sep=""))
fatHDI <- readRDS(paste(datapath, 'case3_berngamma_log_add_fatHDI.RDS', sep=""))
fat_event <- as.data.frame(fat[k,])
names(fat_event) <- c('fat')

event1_ <- matrix(NA, ncol=3, nrow=1)
event1_[1,1] <- expo_cat$obs[k]
event1_[1,2] <- expo_cat$mag_obs[k]
event1_[1,3] <- 0
event1_ <- data.frame(event1_)
names(event1_) <- c('obs', 'x', 'y')

#fatHDI <- HDIofMCMC(fat_event$fat)
fatHDI_ <- fatHDI[k,]
fatHDI_ <- as.data.frame(fatHDI_)
#prob_mag <- readRDS('./prob_mag_no_corr_n100_worden.RDS')

#fatHDI <- fatHDI[k,]
#fatHDI <- as.data.frame(fatHDI)

#show(is.atomic(estimate_fatality))
#names(fatHDI) <- "x"
fatHDI_$y <- c(0.0,0.0)
#fatHDI$obs <- expo_cat$obs[k]
#fatHDI$eqID <- expo_cat$eqID[k]

#show(is.atomic(estimate_fatality))

#fatHDI$obs[expo_cat$obs==0] = NA

myPlot1 <- ggplot(fat_event,aes(x=fat)) +
    #geom_histogram() +
      geom_histogram(aes(y=..count../sum(..count..)), fill='grey') +
          #geom_histogram(aes(y = ..density..), fill='cyan') +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x)),
                   limits = c(10^-1, 10^4)) +
    geom_point(data=event1_,aes(obs,y), col='red', size=3) + #, shape='x'
    geom_line(data = fatHDI_, aes(x=fatHDI_, y=y), size=0.5) +
    geom_point(data = fatHDI_, aes(x=fatHDI_, y=y), shape='|', size=2) +
    ylab("Normalized frequency") +
    xlab("Fatality") +
    theme_bw(base_size=10)

ggsave(paste(plotpath,'figure7a_',saveName,'.', saveType, sep=''),  width = 8, height = 8, unit="cm",  myPlot1)

########
## example showing fatality estimation
########

dat <- read.csv(paste(datapath, 'case1.csv', sep=""), header = FALSE)
names(dat) <- c('pop', 'fat', 'mmi', 'mmi_round', 'id')

#######
# zero prob
#######


