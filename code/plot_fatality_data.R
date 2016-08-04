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


#-------------------------------------------------------------------------------
# Main procedure
#-------------------------------------------------------------------------------

# plot of fatality data
#csvfileName_worden = '/Users/hyeuk/Projects/fatality/data/DATA_WORDEN_COR_ROUND_12_Feb_2013.csv'
datapath = '/Users/hyeuk/Projects/fatality/data/'
plotpath = '/Users/hyeuk/Projects/fatality/plot/'
dat1 <- process_csvdata(paste(datapath,'case1.csv', sep=""))
dat3 <- process_csvdata(paste(datapath,'case3.csv', sep=""))

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
ggsave(paste(plotpath,'figure4a_case1_log.eps', sep=""), width = 8, height = 8, unit="cm",  myPlot8)

myPlot8a <- ggplot(dat1, aes(x=factor(mmi_bin),y=rat)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.3), size=0.5, colour='red')+
  xlab("MMI") +
  ylab("Fatality rate") +
  theme_bw(base_size=10)
ggsave(paste(plotpath,'figure4a_case1_linear.eps', sep=""), width = 8, height = 8, unit="cm",  myPlot8a)


myPlot8b <- ggplot(dat3, aes(x=factor(mmi_bin),y=rat)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.3), size=0.5, colour='red')+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),limits=c(10^-7, 10^-1)) +
  xlab("MMI") +
  ylab("Fatality rate") +
  theme_bw(base_size=10)
ggsave(paste(plotpath,'figure4a_case3.eps', sep=""), width = 12, height = 8, unit="cm",  myPlot8b)

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
dat1$GMICE <- 'BMKG'
dat3$GMICE <- 'EXPO-CAT'

dat1_all <- dat1
dat1_all$id <- 5

dat3_all <- dat3
dat3_all$id <- 5

dat_total <- rbind(dat1, dat3, dat1_all, dat3_all)
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
  scale_x_continuous(limits=c(3.5,9.5), breaks=seq(4.0, 9.0, 1.0)) +
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

ggsave(paste(plotpath,'figure3_by_setting.eps', sep=""),  width = 16, height = 8, unit="cm",  myPlot2)

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






