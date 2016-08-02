graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#require(rstan)
#source("DBDA2E-utilities.R")
library(MASS)
require(ggplot2)
library(scales)
library(reshape2)

process_csvdata <- function(csvfileName) {
  
  dat <- read.csv(csvfileName, header=0)
  names(dat) <- c("pop","fat","mmi","mmi_bin","id")
  
  #dat$fat[dat$fat==0] = 1.0e-3 # very small number instead of 0 fatality
  dat$rat <- dat$fat/dat$pop
  ndat <- dat[dat$fat > 0, ]
  
  #dat.df <- data.frame(dat$mmi, dat$rat)
  #names(dat.df) <- c("MMI","rate")
  
  #saveRDS()
  return(ndat)
}

#library(ggplot2)

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

plotCI <- function(fatHDI, observed, eqID, saveName, saveType) {

  fatHDI <- as.data.frame(fatHDI)

  #show(is.atomic(estimate_fatality))

  names(fatHDI) <- c("L","U")
  fatHDI$obs <- observed
  fatHDI$eqID <- eqID

  #show(is.atomic(estimate_fatality))

  fatHDI$obs[observed==0] = NA
  #estimate_fatality$x <- seq(1, length(observed))

  #show(is.atomic(estimate_fatality))

  #show(names(estimate_fatality))

  # 5%-95% HDI plot

  #ts <- as.character(estimate_fatality$eqID)

  #show(ts)

  #openGraph(6, 8)


  myPlot <- ggplot(fatHDI, aes(y = as.character(eqID), x = obs), xlab = 'eqID') +
    geom_point(size = 4, col='red') +
    geom_errorbarh(aes(xmax = U, xmin = L)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x)),
                   limits = c(10^-1, 10^4)) +
    ylab("Earthquake ID") +
    xlab("Number of fatality") +
    theme_bw() +
    theme(axis.text.y=element_text(angle=0, hjust=1))
  myPlot <- myPlot + annotation_logticks(sides = "tb")                # Default: log ticks on bottom and left
    #scale_x_log10(limits = c(0,1e8)) +

  #return(estimate)

  if ( !is.null(saveName) ) {
    #saveGraph( file=paste(saveName,"CI",sep=""), type=saveType)

    #postscript(file = paste(saveName, 'CI.eps', sep= ''))
    #print(myPlot)
    #dev.off()

    ggsave(paste('CI_',saveName, '.', saveType, sep= ''), width = 18, height = 30, unit="cm",  myPlot)
    #dev.off()
  }
  #return(myPlot)
}

# one event (histogram) comparison approach I vs. II
# plot CI
plot_prob <- function( prob_mag, expo_cat, saveName, saveType, bar_width=0.3, pt_size=2.5, base_size = 10, no_col = 6 ) {

  names_str <- c('0-10^0','10^0-10^1','10^1-10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6','10^7-10^8')

  label_str <- c(expression(0-10^0), expression(10^0-10^1), expression(10^1-10^2),
    expression(10^2-10^3), expression(10^3-10^4), expression(10^4-10^5),
    expression(10^5-10^6), expression(10^6-10^7))

  nevents = nrow(prob_mag)
  nmax = ncol(prob_mag)

  output <- data.frame(prob_mag)
  names(output) <- names_str[1:nmax]

  output$eqID <- expo_cat$eqID

  df.long <- melt(output,id.vars=c("eqID"))

  event1_ <- matrix(NA, ncol=3, nrow=nevents)
  event1_[,1] <- expo_cat$eqID
  event1_[,2] <- expo_cat$mag_obs
  event1_[,3] <- 0*event1_[,2]
  #for (i in 1:nevents) {
  #  event1_[i, 3] <- prob_mag[i, event1_[i, 2]]
  #}
  event1_ <- data.frame(event1_)
  names(event1_) <- c('eqID','x','y')
  #event1_.long <- melt(event1_, id.vars=c("eqID"))

  myPlot <- ggplot(df.long,aes(variable,value))+
     geom_bar(stat="identity", width=bar_width)+
     coord_flip() +
     scale_x_discrete(labels= label_str) +
     geom_point(data=event1_,aes(x,y), col='red', size=pt_size) + #, shape='x'
     facet_wrap(~eqID,ncol=no_col) +
     ylab("Probability") +
     xlab("Fatality ranges") +
     theme_bw(base_size=base_size)

  if ( !is.null(saveName) ) {

    ggsave(paste('pb_events_',saveName, '.', saveType, sep= ''), width = 18, height = 30, unit="cm",  myPlot)

  }
}

# plot PB
plot_prob_one_event <- function( prob_mag, k, expo_cat, saveName, saveType ) {

  names_str <- c('0-10^0','10^0-10^1','10^1-10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6','10^7-10^8')

  label_str <- c(expression(0-10^0), expression(10^0-10^1), expression(10^1-10^2),
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

    ggsave(paste('pb_',saveName, '.', saveType, sep= ''),  width = 8, height = 8, unit="cm",  myPlot)

  }
  return(myPlot)
}

estimate_fatality_PAGER <- function( parameters , pop, mmi_list ) {

  # comparison with Kishor's model
  #params_PAGER = c(14.05, 0.17, 1.74) # USGS OpenFile report 2009-1136
  #params_PAGER = c(14.05, 0.17, 2.15) # EQ spectra (Jaiswal and Wald, 2010)
  #params_PAGER = c(13.249,  0.151, 1.641) # http://pubs.usgs.gov/of/2009/1136/pdf/PAGER%20Implementation%20of%20Empirical%20model.xls(Jaiswal, 2009)
  #estimate_PAGER <- estimate_fatality_PAGER(params_PAGER, pop, mmi_list)
  #plotCI(estimate_PAGER, observed)

  #med_kishor = 13.249; beta_kishor = 0.151;  # v2.0
  #x_trev, y_trev = 1.03, 10.75 # Trevo Allen

  #14.05 0.17 2.15


  # Posterior prediction:
  nevents = nrow(pop)
  nmmi = ncol(pop)

  fat_rate <- plnorm(mmi_list, meanlog=log(parameters[1]), parameters[2])

  # Generate posterior predicted y values.
  # This gets only one y value, at each x, for each step in the chain.

  fat_by_event = pop %*% fat_rate

  # 5% ~ 95% of estimated fatality
  fatHDI <- matrix( NA, nrow = nevents, ncol = 3)
  fatHDI[,1] <- exp(parameters[3]*qnorm(0.05)+log(fat_by_event))
  fatHDI[,2] <- exp(parameters[3]*qnorm(0.95)+log(fat_by_event))
  fatHDI[,3] <- fat_by_event

  order_mag <- ceiling(log10(max(fatHDI[,2])))

  prob_mag <- matrix(NA, nrow = nevents, ncol = (order_mag + 1))

  prob_mag[, 1] = plnorm( 1, meanlog = log(fatHDI[,3]), sdlog = parameters[3] + 0*fatHDI[,3])

  for (j in 2:(order_mag+1)) {
    b2 <- plnorm( 10^(j-1), meanlog = log(fatHDI[,3]), sdlog = parameters[3] + 0*fatHDI[,3])
    b1 <- plnorm( 10^(j-2), meanlog = log(fatHDI[,3]), sdlog = parameters[3] + 0*fatHDI[,3])
    prob_mag[,j] <- b2 - b1
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

#-------------------------------------------------------------------------------
# Main procedure
#-------------------------------------------------------------------------------

# read expo cat for reference
dummy <- read_expo_cat()
pop <- dummy[[1]]
expo_cat <- dummy[[2]]
mmi_list <- dummy[[3]]

# PAGER CAT
params_PAGER = c(13.249,  0.151, 1.641)
dummy <- estimate_fatality_PAGER( params_PAGER , pop, mmi_list )
fatHDI_PAGER <- dummy[[1]]
prob_mag_PAGER <- dummy[[2]]
saveName = 'PAGER'
saveType = 'eps'
plotCI(fatHDI_PAGER[,1:2], expo_cat$obs, expo_cat$eqID, saveName, saveType)
plot_prob( prob_mag, expo_cat, saveName, saveType)

# # # no correlaiton n100
saveName = 'no_corr_n100_worden'
saveType = 'eps'
prob_mag <- readRDS('./prob_mag_no_corr_n100_worden.RDS')
plot_prob( prob_mag, expo_cat, saveName, saveType)

# CI
fatHDI <- readRDS('./fatHDI_no_corr_n100_worden.RDS')
plotCI(fatHDI, expo_cat$obs, expo_cat$eqID, saveName, saveType)

# one event
# histogram of one event (200006041628, Bengkulu, Indonesia) INDEX: 42
k = 42
myPlot3 <- plot_prob_one_event(prob_mag, k, expo_cat, saveName, saveType)

fat <- readRDS('fat_by_event_no_corr_n100_worden_event42.RDS')
fat_event <- as.data.frame(fat)

event1_ <- matrix(NA, ncol=3, nrow=1)
event1_[1,1] <- expo_cat$obs[k]
event1_[1,2] <- expo_cat$mag_obs[k]
event1_[1,3] <- 0
event1_ <- data.frame(event1_)
names(event1_) <- c('obs','x','y')

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
    geom_histogram(aes(y=..count../sum(..count..)), fill='grey') +
          #geom_histogram(aes(y = ..density..), fill='cyan') +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x)),
                   limits = c(10^-1, 10^4)) +
    geom_point(data=event1_,aes(obs,y), col='red', size=1) + #, shape='x'
    geom_line(data = fatHDI_, aes(x=fatHDI_, y=y), line='black', size=1) +
    geom_point(data = fatHDI_, aes(x=fatHDI_, y=y), shape='|', size=2) +
    ylab("Normalized frequency") +
    xlab("Fatality") +
    theme_bw(base_size=10)

ggsave('figure5.eps',  width = 8, height = 8, unit="cm",  myPlot1)

# myPlot2 <- ggplot(fat_event,aes(x=fat)) +
#     geom_histogram(aes(y=..count../sum(..count..)), fill='grey') +
#           #geom_histogram(aes(y = ..density..), fill='cyan') +
#     #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#     #               labels = trans_format("log10", math_format(10^.x)),
#     #               limits = c(10^-1, 10^4)) +
#     geom_point(data=event1_,aes(obs,y), col='red', size=6) + #, shape='x'
#     geom_line(data = fatHDI_, aes(x=fatHDI_, y=y), line='black', size=1) +
#     geom_point(data = fatHDI_, aes(x=fatHDI_, y=y), shape='|', size=8) +
#     ylab("Normalized frequency") +
#     xlab("Fatality") +
#     theme_bw(base_size=10)

# ggsave('figure5b.eps',  width = 8, height = 8, unit="cm",  myPlot2)

# saveName = 'no_corr_n100_wald'
# saveType = 'eps'
# prob_mag <- readRDS('./prob_mag_no_corr_n100_wald.RDS')
# plot_prob( prob_mag, expo_cat, saveName, saveType)

# # # CI
# fatHDI <- readRDS('./fatHDI_no_corr_n100_wald.RDS')
# plotCI(fatHDI, expo_cat$obs, expo_cat$eqID, saveName, saveType)

# # no correlaiton n100
# saveName = 'no_corr_n100'
# saveType = 'eps'
# prob_mag_no_corr_n100 <- readRDS('/Volumes/RYU@GA/nci/prob_mag_no_corr_n100.RDS')
# plot_prob( prob_mag_no_corr_n100, expo_cat, saveName, saveType)

# # CI
# fatHDI_no_corr_n100 <- readRDS('/Volumes/RYU@GA/nci/fatHDI_no_corr_n100.RDS')
# plotCI(fatHDI_no_corr_n100, expo_cat$obs, expo_cat$eqID, saveName, saveType)

# # correlation
# saveName = 'corr_n100'
# saveType = 'eps'
# prob_mag_corr_n100 <- readRDS('/Volumes/RYU@GA/nci/prob_mag_corr_n100.RDS')
# plot_prob( prob_mag_corr_n100, expo_cat, saveName, saveType)

# # CI
# fatHDI_corr_n100 <- readRDS('/Volumes/RYU@GA/nci/fatHDI_corr_n100.RDS')
# fatHDI_corr_n100[fatHDI_corr_n100[,1]< 1.0e-1,1] = 0.1
# plotCI(fatHDI_corr_n100, expo_cat$obs, expo_cat$eqID, 'xxx', saveType)

# plot of fatality data
#csvfileName_worden = '/Users/hyeuk/Projects/fatality/data/DATA_WORDEN_COR_ROUND_12_Feb_2013.csv'
datapath = '/Users/hyeuk/Projects/fatality/data/'
dat1 <- process_csvdata(paste(datapath,'case1.csv', sep=""))
dat2 <- process_csvdata(paste(datapath,'case2.csv', sep=""))

# exploratory analysis
# dat <- dat_worden

#sd_ <- tapply(dat$rat, dat$mmi_bin, sd)

#mean_ <- tapply(dat$rat, dat$mmi_bin, mean)

#X11()
#qplot(names(sd_), log(mean_)) + xlab("MMI") + ylab("log mean of fatality rate")
  #dev.copy2eps(file='log_mean_fatality.eps')

# myPlot7 <- ggplot(dat1, aes(x=mmi,y=rat)) +
#  geom_point() +
#  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
#  #geom_smooth() +
#  #scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
#  #                  labels = trans_format("log10", math_format(10^.x)),
#  #                  limits = c(10^-12, 10^-3)) +
#  xlab("MMI") +
#  ylab("Fatality rate") +
#  theme_bw(base_size=10)
# 
# ggsave('figure4b.eps',  width = 8, height = 8, unit="cm",  myPlot7)


myPlot8 <- ggplot(dat1, aes(x=factor(mmi_bin),y=rat)) +
 geom_boxplot() +
 geom_point(position = position_jitter(width = 0.3), size=0.5, colour='red')+
 scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
 xlab("MMI") +
 ylab("Fatality rate") +
 theme_bw(base_size=10)
ggsave('figure4a.eps',  width = 12, height = 8, unit="cm",  myPlot8)

# mf_labeller <- function(var, value){
#     value <- as.character(value)
#     if (var=="id") {
#         value[value=="1"]   <- "Nias"
#         value[value=="2"]   <- "Yogyakarta"
#         value[value=="3"]   <- "West Sumatra"
#         value[value=="4"]   <- "West Java"
#         #value[value=="5"]   <- "All events"
#     }
#     return(value)
# }
# 
# myPlot9 <- ggplot(dat_total, aes(x=mmi_bin,y=rat)) +
#  geom_boxplot() +
#  facet_grid(~id, labeller = mf_labeller) +
#  xlab("MMI") +
#  ylab("Fatality rate") +
#  theme_bw(base_size=10)
# ggsave('figure9.eps',  width = 8, height = 8, unit="cm",  myPlot8)

# aggregate
dat1$GMICE <- 'Finite Fault'
dat2$GMICE <- 'Point Source'

dat1_all <- dat1
dat1$id <- 5

dat2_all <- dat2
dat2_all$id <- 5

dat_total <- rbind(dat1, dat2, dat1_all, dat2_all)
dat_total$id <- factor(dat_total$id)

# mf_labeller <- function(var, value){
#     value <- as.character(value)
#     if (var=="id") {
#         value[value=="1"]   <- "Nias"
#         value[value=="2"]   <- "Yogyakarta"
#         value[value=="3"]   <- "West Sumatra"
#         value[value=="4"]   <- "West Java"
#         value[value=="5"]   <- "All events"
#     }
#     return(value)
# }

label_metrics <- function(x) {
    x[x == "1"]   <- "Nias"
    x[x == "2"]   <- "Yogyakarta"
    x[x == "3"]   <- "West Sumatra"
    x[x == "4"]   <- "West Java"
    x[x == "5"]   <- "All events"
    x
  }

param_labeller <- ggplot2::as_labeller(label_metrics)

#myPlot <- ggplot(dat_total,aes(x=mmi, y=rat, shape=id))+

myPlot2 <- ggplot(data=dat_total, aes_string(x='mmi', y='rat')) +
 facet_grid(GMICE~id, labeller = param_labeller) +
 geom_point(size=0.5) +      # Use hollow circles
 scale_x_continuous(limits=c(3.5,8.5)) +
 scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
  #facet_wrap(~GMICE,ncol=2) +
 #scale_shape_manual(values=c(15,16,17,18)) + # Use a hollow circle and triangle
 #scale_shape_discrete(name="EQ event",
 #                         breaks=c("1", "2", "4", "3"),
 #                        labels=c("Nias", "Yogyakarta", "West Java","West Sumatra"), solid=FALSE) +
 #guides(color=guide_legend(override.aes=list(fill=NA))) +
 #guides(colour = guide_legend(nrow = 2)) +
 xlab("MMI") +
 ylab("Fatality rate") +
 theme_bw(base_size=10)

ggsave('figure3.eps',  width = 12, height = 8, unit="cm",  myPlot2)

 #theme(legend.position="right", legend.key = element_rect(colour = NA),
 # legend.title = element_text(size = 6),
 # legend.text=element_text(size=6))







# myPlot <- ggplot(dat_total,aes(x=mmi, y=rat, shape=id))+
#  geom_point(size=1) +      # Use hollow circles
#  scale_x_continuous(limits=c(3.5,9.0)) +
#  facet_wrap(~GMICE,ncol=2) +
#  #scale_shape_manual(values=c(15,16,17,18)) + # Use a hollow circle and triangle
#  scale_shape_discrete(name="EQ event",
#                          breaks=c("1", "2", "4", "3"),
#                          labels=c("Nias", "Yogyakarta", "West Java","West Sumatra"), solid=FALSE) +
#  #guides(color=guide_legend(override.aes=list(fill=NA))) +
#  #guides(colour = guide_legend(nrow = 2)) +
#  xlab("MMI") +
#  ylab("Fatality rate") +
#  theme_bw(base_size=10) +
#  theme(legend.position="right", legend.key = element_rect(colour = NA),
#   legend.title = element_text(size = 6),
#   legend.text=element_text(size=6))
# 
# ggsave('figure6.eps',  width = 16, height = 8, unit="cm",  myPlot)


 # geom_bar(stat="identity", width=0.5, fill='grey')+
 # coord_flip() +
 # scale_x_discrete(labels= label_str) +
 # geom_point(data=event1_,aes(x,y), col='red', size=4) +
 # #facet_wrap(~eqID,ncol=1) +
 # ylab("Probability") +
 # xlab("Fatality ranges") +
 # theme_bw(base_size=10)






