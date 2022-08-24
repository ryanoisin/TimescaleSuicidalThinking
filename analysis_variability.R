# o.ryan@uu.nl 17 August

# Reproduces the "Variability as a function of time-interval" analysis

source("aux_functions.R")
library(dplyr)
library(collapse)

# load data
data <- readRDS(file = paste0("Data/SI_EMA_01_2022.RDS"))
ids <- unique(data$ID)

# put data into person-specific list
ldata <- lapply(ids,function(i){
  data[data$ID == i,]
})

# overwrite 
for(i in 1:length(ids)){
  ldata[[i]]$ID <- i
}


nobs <- unlist(lapply(ldata,nrow))

# ------------------------------------------------------------
# ----- Quantify variability using the mode (descriptives) ---
# ------------------------------------------------------------

# note: here p_mode is defined as the proportion of observations
# which are equal to the mode

get_pmode <- function(l){
  c(sum(l$SI_DesireKill == fmode(l$SI_DesireKill), na.rm = TRUE)/nrow(l),
    sum(l$SI_Intent== fmode(l$SI_Intent), na.rm = TRUE)/nrow(l))
}

# obtian the modal value itself
get_mode <- function(l){
  c(fmode(l$SI_DesireKill),fmode(l$SI_Intent))
}

mode <- t(sapply(ldata, get_mode))
# how often do values appear as the mode?
table(mode[,1]); table(mode[,2])
# cross-table: co-occurence of mode values
table(mode[,1],mode[,2])

# calcualte proportion of values equal to the mode
pmode <- t(sapply(ldata, get_pmode))
# # desire
 hist(pmode[,1]); mean(pmode[,1]); quantile(pmode[,1], c(.25, .5, .75))
# # intent
 hist(pmode[,2]); mean(pmode[,2]); quantile(pmode[,2], c(.25, .5, .75))

# how many individuals show no variation on each variable?
sum(pmode[,1] == 1)
sum(pmode[,2] == 1)

# correlation
cor(pmode)

# ------------------------------------------------------------
# ----------------- Analyse variability: Part I --------------
# ------------------------------------------------------------

# First, look for ways to cluster participants together

# Find those individuals who have < 5 percent observations not equal to the mode
tf_pmode <- apply(pmode, 1, function(row) any(row > .95))
sum(tf_pmode)

# variable specific
d_pmode <- apply(pmode, 1, function(row) row[1] > .95)
i_pmode <- apply(pmode, 1, function(row) row[2] > .95)

sum(d_pmode)
sum(i_pmode)

# here we see our three groups, lower bottom-right cell is empty
table(d_pmode, i_pmode)

# ---------- examine low variance people for DESIRE -----#

# extract individuals with low variance on desire
lowd <- ldata[d_pmode]

# for(i in 1:length(lowd)){
#   tmp <- lowd[[i]]
#   plot.new()
#   plot.window(xlim = c(0,max(tmp$RunTime)), ylim = c(0,10))
#   axis(1); axis(2)
#   title(main = paste0("ID = ", i))
#   lines(tmp[,"SI_DesireKill"], x= tmp$RunTime, col = "blue", type = "b")
#   lines(tmp[,"SI_Intent"],  x= tmp$RunTime, col = "red", type = "b")
# }
# # 

# Detect episodes of elevated desire
  # see aux_functions.R for details of how this is computed
elist <- lapply(lowd, function(l){
  detectelev(l$SI_DesireKill, l$RunTime, l$Type, l$Time)$elev_stats
})

elist
sum(unlist(lapply(elist,nrow))!=0)
# 9 participants show episodes

# extract episode durations
durations <- unlist(lapply(elist, function(l) l[,"duration"]))
durations <- na.omit(durations)
# hist(durations)

# Summarize episode durations
mean(durations)
quantile(durations, c(0.25, .5, .75))

# sum(durations < 3)/length(durations)

# ---------- examine low variance people for INTENT -----#

lowi <- ldata[i_pmode]
# for(i in 1:length(lowi)){
#   tmp <- lowi[[i]]
#   plot.new()
#   plot.window(xlim = c(0,max(tmp$RunTime)), ylim = c(0,10))
#   axis(1); axis(2)
#   title(main = paste0("ID = ", i))
#   lines(tmp[,"SI_DesireKill"], x= tmp$RunTime, col = "blue", type = "b")
#   lines(tmp[,"SI_Intent"],  x= tmp$RunTime, col = "red", type = "b")
# }
# 
# lowi_e <- lapply(lowi, function(l){
#   l$ind_label <- detectelev(l$SI_Intent, l$RunTime, l$Type, l$Time)$ind_label
#   return(l)
# })

ilist <- lapply(lowi, function(l){
  detectelev(l$SI_Intent, l$RunTime, l$Type, l$Time)$elev_stats
})
ilist
ilist_s <- ilist[unlist(lapply(ilist,nrow))!=0]
ilist_s <- lapply(ilist_s, function(l){
  # some participants have one episode, transpose to correct format
  if(ncol(l) ==1) as.data.frame(t(l)) else l
})

length(ilist_s)
# 26 parcticipants

i_elev_df <- do.call("rbind", ilist_s)
nrow(i_elev_df)
# 67 episodes


durations <- unlist(lapply(ilist_s, function(l) l$duration))
length(durations)
durations <- na.omit(durations)
hist(durations)
mean(durations); quantile(durations, c(.25,.5, .75))
median(durations)
# sum(durations < 3)/length(durations)

# ------------------------------------------------------------
# ----------------- Analyse variability: Part II -------------
# ------------------------------------------------------------

# we now examine the remaining participants, who exhibit relatively
# higher variability in desire and/or intent

# ----------------- Desire -------------

# additional data cleaning step: omit those with very short time-series
tsl <- unlist(lapply(ldata, nrow))
dropl <- which(tsl < 14)

# omit those who have very low variability in desire 
dropd <- which(d_pmode)
drop <- unique(c(dropd, dropl))
length(drop); length(ids)

ldata_d <- ldata[-drop]
# manually dropping one strange participant who only has short and very long gaps
ldata_d <- ldata_d[-70]
length(ldata_d)
# 85 participants

# idenfity whether Desire has changed in value from one occasion to the next
change_d <- lapply(ldata_d, function(l){
  dt <- l$RunTime[-1] - l$RunTime[-nrow(l)]
  c_d <- l$SI_DesireKill[-1] - l$SI_DesireKill[-nrow(l)] != 0
  data.frame(ID = l$ID[-1], dt = dt, c_d = c_d)
})

# check number of observations in each bin
unlist(lapply(change_d, function(l) length(which(l$dt < .75)))) 
unlist(lapply(change_d, function(l) length(which(l$dt >= .75 & l$dt <= 3)))) 
unlist(unlist(lapply(change_d, function(l) length(which(l$dt > 3)))))

# check that all individuals have at least two observations in each bin
which(unlist(lapply(change_d, function(l) length(which(l$dt >= .75 & l$dt <= 3)))) < 2)


# proportion of observation pairs that show variability in each bin
pp1_d <- unlist(lapply(change_d, function(l){
  tmp  <- l[l$dt < .75,]
  sum(tmp$c_d, na.rm = T)/nrow(tmp)
}))

pp2_d <- unlist(lapply(change_d, function(l){
  tmp  <- l[l$dt >= .75 & l$dt <= 3,]
  sum(tmp$c_d, na.rm = T)/nrow(tmp)
}))

pp3_d <- unlist(lapply(change_d, function(l){
  tmp  <- l[l$dt > 3,]
  sum(tmp$c_d, na.rm = T)/nrow(tmp)
}))


# collect in matrix
ppmat <- cbind(pp1_d,pp2_d,pp3_d)

# descriptives
mean(pp1_d); mean(pp2_d); mean(pp3_d)



# -- Plot for Desire -------------
set.seed(123)
# jitter slightly on x-axis for visibility
xpp1 <- jitter(rep(1, length(pp1_d)))
xpp2 <- jitter(rep(2, length(pp2_d)))
xpp3 <- jitter(rep(3, length(pp3_d)))




library(scales)
pdf("figures/ts_desire.pdf",4.5,4.5)
plot.new()
plot.window(xlim = c(0.9,3.1), ylim = c(0,1))
axis(2); axis(1, at = c(1,2,3), labels = c("<45 min", "45 min - 3hrs", "> 3 hrs"))
title(ylab = "Proportion of change", xlab = "Timescale",
      # main = "Desire")
      main = "")
points(y = pp1_d, x =xpp1, col = alpha("red",0.3), pch = 18)
points(y = pp2_d, x =xpp2, col = alpha("red",0.3), pch = 18)
points(y = pp3_d, x =xpp3, col = alpha("red",0.3), pch = 18)


# see how people change

for(i in 1:nrow(ppmat)){
  lines(y = ppmat[i,], x = c(1,2,3), type = "l", col = alpha("grey",0.5))
}

points(y = mean(pp1_d), 1, col = alpha("black",0.5), pch = 18, cex = 2)
points(y = mean(pp2_d), 2, col = alpha("black",0.5), pch = 18, cex = 2)
points(y = mean(pp3_d), 3, col = alpha("black",0.5), pch = 18, cex = 2)

dev.off()


# -----------------------------------------------------------------------------
# ----------------- Intent -------------
# repeat analysis for intent


# omit those who have very low variability in desire 

dropi <- which(i_pmode)
drop_intent <- unique(c(dropi, dropl))
length(drop_intent); length(ids)

ldata_i <- ldata[-drop_intent]
length(ldata_i)
# 59 participants

# idenfity whether Desire has changed in value from one occasion to the next
change_i <- lapply(ldata_i, function(l){
  dt <- l$RunTime[-1] - l$RunTime[-nrow(l)]
  c_i <- l$SI_Intent[-1] - l$SI_Intent[-nrow(l)] != 0
  data.frame(ID = l$ID[-1], dt = dt, c_i = c_i)
})

# check number of observations in each bin
unlist(lapply(change_i, function(l) length(which(l$dt < .75)))) 
unlist(lapply(change_i, function(l) length(which(l$dt >= .75 & l$dt <= 3)))) 
unlist(unlist(lapply(change_i, function(l) length(which(l$dt > 3)))))

# check that all individuals have at least two observations in each bin 
drop2 <- which(unlist(lapply(change_i, function(l) length(which(l$dt >= .75 & l$dt <= 3)))) < 2)
# drop 1 individual who only has very long or very short gaps
change_i <- change_i[-drop2]

# proportion of observation pairs that show variability in each bin
pp1_i <- unlist(lapply(change_i, function(l){
  tmp  <- l[l$dt < .75,]
  sum(tmp$c_i, na.rm = T)/nrow(tmp)
}))

pp2_i <- unlist(lapply(change_i, function(l){
  tmp  <- l[l$dt >= .75 & l$dt <= 3,]
  sum(tmp$c_i, na.rm = T)/nrow(tmp)
}))

pp3_i <- unlist(lapply(change_i, function(l){
  tmp  <- l[l$dt > 3,]
  sum(tmp$c_i, na.rm = T)/nrow(tmp)
}))


# collect in matrix
ppmat <- cbind(pp1_i,pp2_i,pp3_i)

# descriptives
mean(pp1_i); mean(pp2_i); mean(pp3_i)



# -- Plot for Intent-------------
set.seed(123)
# jitter slightly on x-axis for visibility
xpp1 <- jitter(rep(1, length(pp1_i)))
xpp2 <- jitter(rep(2, length(pp2_i)))
xpp3 <- jitter(rep(3, length(pp3_i)))




library(scales)
pdf("figures/ts_intent.pdf",4.5,4.5)
plot.new()
plot.window(xlim = c(0.9,3.1), ylim = c(0,1))
axis(2); axis(1, at = c(1,2,3), labels = c("<45 min", "45 min - 3hrs", "> 3 hrs"))
title(ylab = "Proportion of change", xlab = "Timescale",
      # main = "Desire")
      main = "")
points(y = pp1_i, x =xpp1, col = alpha("red",0.3), pch = 18)
points(y = pp2_i, x =xpp2, col = alpha("red",0.3), pch = 18)
points(y = pp3_i, x =xpp3, col = alpha("red",0.3), pch = 18)


# see how people change

for(i in 1:nrow(ppmat)){
  lines(y = ppmat[i,], x = c(1,2,3), type = "l", col = alpha("grey",0.5))
}

points(y = mean(pp1_i), 1, col = alpha("black",0.5), pch = 18, cex = 2)
points(y = mean(pp2_i), 2, col = alpha("black",0.5), pch = 18, cex = 2)
points(y = mean(pp3_i), 3, col = alpha("black",0.5), pch = 18, cex = 2)

dev.off()
