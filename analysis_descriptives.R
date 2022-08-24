# o.ryan@uu.nl 17 August

# This file loads the raw data and creates the descriptive statistics 
  #  - numbers of observations of different types
  #  - response, mean and SD distributions
  #  - high risk response rates

library(RColorBrewer)
library(scales)
library(plyr)
library(dplyr)
source("aux_functions.R")

#load data
data <- readRDS(file = paste0("Data/SI_EMA_01_2022.RDS"))

# create empty highRisk variable for later use
data$HighRisk <- 0

#split data into ema and burst data frames
ema <- dplyr::filter(data, Type=="EMA")
burst <- dplyr::filter(data, Type=="Burst")

# put data into person-specific list
ids <- unique(data$ID)
ldata <- lapply(ids,function(i){
  data[data$ID == i,]
})

for(i in 1:length(ids)){
  ldata[[i]]$ID <- i
}

# ----------------------------------------------------------------------
# ------------------------- Basic Descriptives -------------------------
# ----------------------------------------------------------------------

# number of observations
nobs <- unlist(lapply(ldata,nrow))
sum(nobs)
range(nobs)
# IQR number of observations
quantile(nobs, c(.25,.75))

# How many days on average?
ndays <- unlist(lapply(ldata, function(l) max(l$RunTime)/24))

mean(ndays); range(ndays); quantile(ndays, c(.25,.75))
  hist(ndays)

# ----------------------------------------------------------------------
# ------------------ Descriptives Visualization  -----------------------
# ----------------------------------------------------------------------

# Prepare color palettes for plotting

colvec <- brewer.pal(4, "Spectral")
col2 <- brewer.pal(11,"Spectral")[c(1,11)]



# -------- First, describe the frequency of responses collapsed across all participants -----
# desire and intent EMA
d1 <- density(ema$SI_DesireKill, na.rm = TRUE, bw= .4)
d1$y <- c(0,d1$y[d1$x >0 & d1$x < 10],0)
d1$x <- c(0,d1$x[d1$x >0 & d1$x < 10],10)

d2 <- density(ema$SI_Intent, na.rm = TRUE, bw= .4)
d2$y <- c(0,d2$y[d2$x >0 & d2$x < 10],0)
d2$x <- c(0,d2$x[d2$x >0 & d2$x < 10],10)


# desire and intent Burst
d1s <- density(burst$SI_DesireKill, na.rm = TRUE, bw= .4)
d1s$y <- c(0,d1s$y[d1s$x >0 & d1s$x < 10],0)
d1s$x <- c(0,d1s$x[d1s$x >0 & d1s$x < 10],10)

d2s <- density(burst$SI_Intent, na.rm = TRUE, bw= .4)
d2s$y <- c(0,d2s$y[d2s$x >0 & d2s$x < 10],0)
d2s$x <- c(0,d2s$x[d2s$x >0 & d2s$x < 10],10)


# ----- Second, calculate all person-specific means and SDs of Desire and Intent -----

#means
ema_data_m <- plyr::ddply(ema, .(ID), summarise, mean_d = mean(SI_DesireKill, na.rm = TRUE), 
                          mean_i = mean(SI_Intent, na.rm = TRUE))
burst_data_m <- plyr::ddply(burst, .(ID), summarise, mean_d = mean(SI_DesireKill, na.rm = TRUE), 
                            mean_i = mean(SI_Intent, na.rm = TRUE))

#SDs
ema_data_sd <- plyr::ddply(ema, .(ID), summarise, sd_d = sd(SI_DesireKill, na.rm = TRUE), 
                           sd_i = sd(SI_Intent, na.rm = TRUE))
burst_data_sd <- plyr::ddply(burst, .(ID), summarise, sd_d = sd(SI_DesireKill, na.rm = TRUE), 
                             sd_i = sd(SI_Intent, na.rm = TRUE))

# Create figure showing densities of responses and person-specific means + SDs
pdf("figures/descriptives_v1.pdf")
par(cex.axis=1, cex.lab=1, mfrow = c(2,2))
cex_r = 1

# make plots
#  ----------- densities ------------
plot.new()
  plot.window(xlim = c(-0.1,10.1), ylim = c(0,.7))
  axis(1)
  axis(2)
    lines(d1, xlim = c(0,10), col = col2[1])
    polygon(d1, col = alpha(col2[1],.5), lty = 0)
    lines(d2, xlim = c(0,10), col = col2[2])
    polygon(d2, col = alpha(col2[2],.5), lty = 0)

  title(xlab = "Response", ylab= "Density")
  legend(
    "topright",
    c("Desire","Intent"),
    fill = alpha(col2,.5),
    border = col2,
    xpd = TRUE,
    bty = "n"
  )

# second density
plot.new()
  plot.window(xlim = c(-0.1,10.1), ylim = c(0,.7))
  axis(1)
  axis(2)
    lines(d1s, xlim = c(0,10), col = col2[1])
    polygon(d1s, col = alpha(col2[1],.5), lty = 0)
    lines(d2s, xlim = c(0,10), col = col2[2])
    polygon(d2s, col = alpha(col2[2],.5), lty = 0)

  title(xlab = "Response", ylab= "Density")
  legend(
    "topright",
    c("Desire","Intent"),
    fill = alpha(col2,.5),
    border = col2,
    xpd = TRUE,
    bty = "n"
  )

# ----------- Mean + SD scatterplots ------------

# mean vs mean
plot.new()
  plot.window(c(0,10), c(0,10))
  axis(1); axis(2)
  title(xlab = "Person-Specific Mean EMA", ylab = "Person-Specific Mean Burst")
  abline(lm(burst_data_m$mean_d ~ ema_data_m$mean_d), col = col2[1], lwd =1)
    points(ema_data_m$mean_d, burst_data_m$mean_d, col = alpha(col2[1],.6), pch = 19)
    text(x = 9.25, y = 1, paste0("r = ",round(cor(ema_data_m$mean_d, burst_data_m$mean_d),2)), 
         col = col2[1], cex = cex_r)

  abline(lm(burst_data_m$mean_i ~ ema_data_m$mean_i), col = col2[2], lwd = 1)
    points(ema_data_m$mean_i, burst_data_m$mean_i, col = alpha(col2[2],.6), pch = 19)
    text(x = 9.25, y = 1.75, paste0("r = ",round(cor(ema_data_m$mean_i, burst_data_m$mean_i),2)), 
         col = col2[2], cex = cex_r)


# sd vs sd
plot.new()
  plot.window(c(0,5), c(0,5))
  axis(1); axis(2)
  title(xlab = "Person-Specific SD EMA", ylab = "Person-Specific SD Burst")
  abline(lm(burst_data_sd$sd_d ~ ema_data_sd$sd_d), col = col2[1], lwd =1)
      points(ema_data_sd$sd_d, burst_data_sd$sd_d, col = alpha(col2[1],.6), pch = 19)
      text(x = 4.75 -.125, y = .5, paste0("r = ",round(cor(ema_data_sd$sd_d, burst_data_sd$sd_d, "pairwise.complete.obs"),2)), 
           col = col2[1], cex = cex_r)

  abline(lm(burst_data_sd$sd_i ~ ema_data_sd$sd_i), col = col2[2], lwd =1)
      points(ema_data_sd$sd_i, burst_data_sd$sd_i, col = alpha(col2[2],.6), pch = 19)
      text(x = 4.75-.125, y = .5 +.375, paste0("r = ",round(cor(ema_data_sd$sd_i, burst_data_sd$sd_i, "pairwise.complete.obs"),2)), 
           col = col2[2], cex = cex_r)

dev.off()

# ----------------------------------------------------------------------
# ---------------------- High-Risk Observations  -----------------------
# ----------------------------------------------------------------------

# code observations as high risk
ldata <- lapply(ldata,function(matrix){
  hr <- which(matrix[,"SI_Intent"] > 7)
  if(length(hr)>0){
    matrix[hr,"HighRisk"] <- 1
  }
  matrix
})

# does a subject produce a "high risk" observation?
tfh_i <- unlist(lapply(ldata,function(mat){
  any(mat[,"HighRisk"] == 1)
}))

# how many high risk observations per person?
nh_i <- unlist(lapply(ldata,function(mat){
  sum(mat[,"HighRisk"] == 1)
}))


sum(tfh_i)# 31 participants
sum(nh_i) # 1213 observations

# ----------------------------------------------------------------
# Identify when a burst measurement captures a high-risk observation
# which is preceded by a "non high-risk" EMA measurement 
# ----------------------------------------------------------------

# Get list with burst information per person
bldata <- get_burst(ldata, include_end = FALSE)

# hr_u is a list with indicators for whether the burst episodes are "uniquely" informative 
# about high risk observations
hr_u <- lapply(bldata,function(l){
  store <- rep(0,length(unique(l$burst_episode))) 
  for(b in unique(l$burst_episode)){
    sel <- l[l[,"burst_episode"]==b,]
    if(any(sel[sel[,"Type2"]!="Start","HighRisk"]==1)){
      if(sel[sel[,"Type2"]=="Start","HighRisk"]==0){
        store[b] <- 1
      }
    }
  }
  store
  
})

# sum "unique" high risk burst observations per person, combine
n_unique <- do.call("rbind",lapply(hr_u,sum))

# 74 observations burst episodes give a high risk observation not found in the preceding EMA measurement
sum(n_unique)

# 18 people have a high risk observation captured by a burst episode (not captured by the preceding EMA)
length(which(n_unique!=0))

# How many people have a high risk observation captured by a burst
# where OTHERWISE we would never capture any high risk observation?
EMA <- filter(data, Type == "EMA")
EMA_risk <- filter(EMA, SI_Intent > 7)
dim(EMA_risk)
nrow(EMA)

nrow(EMA_risk)/nrow(EMA)

# burst
burst <- filter(data, Type == "Burst")
burst_risk <- filter(burst, SI_Intent > 7)
dim(burst_risk)
nrow(burst)

nrow(burst_risk)/nrow(burst)

#544 high risk responses 
ids_ema_risk <- unique(EMA_risk$ID)

ids_bu_risk <- ids[which(n_unique!=0)]
length(ids_ema_risk)

# 6 people have bursts pick up high risk info 
# where otherwise no high risk info would be found
sum(!(ids_bu_risk %in% ids_ema_risk))






