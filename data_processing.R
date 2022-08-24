# o.ryan@uu.nl

# this file loads the complete dataset and;
    # - processes the data to be appropriate for ctsem
    # - creates subsets for the ctsem and msm subset analysis
    # - outputs RDS data files


# ----------------------------------------------------------------------
# --------- Source helper functions, Load Data -------------------------
# ----------------------------------------------------------------------

source("aux_functions.R")

library(collapse)

data <- readRDS(file = "Data/SI_EMA_01_2022.RDS")

nobs <- table(data$ID)
ids <- unique(data$ID)
nid <- length(ids)

# how many "burst observations" per person?
nb <- sapply(ids,function(i)
  sum(data$Type[data$ID==i] == "Burst"))


# ----------------------------------------------------------------------
# --------------- Data cleaning and selection  -------------------------
# ----------------------------------------------------------------------


# Save each individual persons new dataset in a list
  # create subsets while doing this
data_list_out <- sapply(ids, function(i){
  # Select an individual
  tmp <- data[data$ID == i,]
  # Calculate consecutive time-differences
  tmp$Difftime <- c(0,tmp[-1,"RunTime"] - tmp[-nrow(tmp),"RunTime"])
  tmp$Diffafter <- c(tmp$Difftime[-1],9999)
  # Now I have, for each row, how long ago the previous observation was, and how long ago the next observation is
  # I want to label as "short" all observations which either follow or are followed by 
  # an observation with a spacing of less than 1.5 hrs
  tmp$Type2 <- "Long"
  tmp[tmp$Difftime < 1.5 | tmp$Diffafter < 1.5, "Type2" ] <- "Short"
  # But we also want to include 
  # tmp[,-9]
  tmp
  
}, simplify = FALSE)

# Select cases with sufficient data length

# drop individuals with short time series
drop_ind1 <- which(sapply(data_list_out, nrow) < 10)

alldata <- do.call("rbind", data_list_out)[,c("ID","SI_DesireKill", "SI_Intent","RunTime")]
data_drop_low_n <- do.call("rbind", data_list_out[-drop_ind1])[,c("ID","SI_DesireKill", "SI_Intent","RunTime")]

# saveRDS(alldata,"Data/data_full.RDS")
# saveRDS(data_drop_low_n,"Data/data_full_low_n_drop.RDS")

# check the mode
get_pmode <- function(l){
  c(sum(l$SI_DesireKill != fmode(l$SI_DesireKill), na.rm = TRUE)/nrow(l),
    sum(l$SI_Intent!= fmode(l$SI_Intent), na.rm = TRUE)/nrow(l))
}

# Now, let's look for people who have very little variance on either variable
# since we have a discrete scale we can look at the mode

pmode <- sapply(data_list_out,get_pmode)
# colnames(pmode) <- c("Desire_ema","Intent_ema","Desire_short","Intent_short")

# Find those individuals who have < 5 percent observations not equal to the mode
tf_pmode <- apply(pmode, 2, function(row) any(row < .05))
drop_ind2 <- which(tf_pmode)

drop_ind <- unique(c(drop_ind1, drop_ind2))
length(drop_ind); length(ids)

# remove problematic cases
data_list_clean <- data_list_out[-drop_ind]
length(data_list_clean)


data_list_std <- lapply(data_list_clean, function(tmp){
  tmp <- tmp[,c("ID","SI_DesireKill", "SI_Intent","RunTime")]
  means <- apply(tmp[,c("SI_DesireKill","SI_Intent")],2,mean)
  tmp[,c("SI_DesireKill","SI_Intent")] <- apply(tmp[,c("SI_DesireKill","SI_Intent")],2,scale)
  tmp
})


# save both versions: the cleaned unstandardized and standardized

data_all_list <- list(data_clean = 
                        do.call("rbind",data_list_clean)[,c("ID","SI_DesireKill", "SI_Intent","RunTime")],
                      data_std = do.call("rbind",data_list_std))

saveRDS(data_all_list, "/ct_analysis_alldata/data_full_ctsem.RDS")


# ----------------------------------------------------------------------
# ------------ Create files for subset analysis ------------------------
# ----------------------------------------------------------------------

# Create subsets for "short" and "long" spacing
data_short_list <- lapply(data_list_out,function(df) {
  tmp <- df[df$Type2=="Short",]
  # Make new runtime and difftime variables
  tmp$RunTime <- round(as.numeric(difftime(tmp$Time,tmp$Time[1], units = "hours")),3)
  tmp$Difftime <- c(0,tmp[-1,"RunTime"] - tmp[-nrow(tmp),"RunTime"])
  tmp
})

# Note that there will be some overlap here, because EMA measurements which directly 
# precede a burst are counted in "short"
# I.e., all "long" observations are ema, but not all "ema" measurements are "long"
data_ema_list <- lapply(data_list_out,function(df) {
  tmp <- df[df$Type=="EMA",]
  # Make new runtime and difftime variables
  tmp$RunTime <- round(as.numeric(difftime(tmp$Time,tmp$Time[1], units = "hours")),3)
  tmp$Difftime <- 0
  tmp[-1,]$Difftime <- tmp[-1,"RunTime"] - tmp[-nrow(tmp),"RunTime"]
  tmp
})

# # now let's get some descriptives for each person across both datasets
# hist(sapply(data_ema_list, nrow))
# hist(sapply(data_short_list,nrow))

# Drop those people who have less than 10 observations in either dataset
drop_l1 <- which(sapply(data_ema_list, nrow) < 10)
drop_l2 <- which(sapply(data_short_list, nrow) < 10)
drop_l <- unique(c(drop_l1,drop_l2))


# Now, let's look for people who have very little variance on either variable
# since we have a discrete scale we can look at the mode
# Note - there are NAs! Set na.rm = FALSE above to see where

pmode <- t(rbind(sapply(data_ema_list,get_pmode),
      sapply(data_short_list, get_pmode)))
colnames(pmode) <- c("Desire_ema","Intent_ema","Desire_short","Intent_short")

# Find those individuals who have < 5 percent observations not equal to the mode
tf_pmode <- apply(pmode, 1, function(row) any(row < .05))
drop_m <- which(tf_pmode)

drop <- unique(c(drop_l, drop_m))
length(drop); length(ids)

# create new dataframes/lists
data_1 <- data[!(data$ID %in% ids[drop]),] # dropped 55
data_ema_list1 <- data_ema_list[-drop]
data_short_list1 <- data_short_list[-drop]

# Check how many burst observations we have per person in the new dataset
nobs_1 <- table(data_1$ID)
ids_1 <- unique(data_1$ID)
nid_1 <- length(ids_1)

nb_1 <- sapply(ids_1,function(i)
  sum(data_1$Type[data_1$ID==i] == "Burst"))

# hist(nobs_1)
# hist(nb_1/nobs_1)

# check the relative size of the short dataset with the ema/long dataset
nobs_short <- do.call("c",lapply(data_short_list1,nrow))
nobs_ema <- do.call("c",lapply(data_ema_list1,nrow))

# hist(nobs_short/nobs_ema, main = "Ratio of Short timediff observations to EMA observations")

# Let's check the time-intervals which are present in all datasets
data_short_df <- do.call("rbind",data_short_list1)
data_ema_df <- do.call("rbind",data_ema_list1)


# quantile(data_short_df$Difftime); quantile(data_short_df$Difftime, c(.83))
#  hist(data_short_df$Difftime)
# 
# quantile(data_ema_df$Difftime); quantile(data_ema_df$Difftime, c(.05))
# hist(data_ema_df$Difftime)

saveRDS(data_ema_df, file="Data/data_ema_50.RDS")
saveRDS(data_short_df, file="Data/data_short_50.RDS")




