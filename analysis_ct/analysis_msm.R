# code to run the markov switching model analyses in the main text

#load packages
library(msm)
library(tidyverse)
library(DataCombine)
source("aux_functions.R")
library(qgraph)


# function to re-code continuous variables into states
get_states <- function(vector,th1 = 4 ,th2 = 7){
  out <- rep(NA,length(vector))
  out[vector < th1 | vector ==th1] <- 1
  out[vector > th1 & vector < th2 | vector ==th2 ] <- 2
  out[vector > th2] <- 3
  out
}


#### load data for full data analysis
data_list <- readRDS(file = "data/data_full_ctsem.RDS")
data <- data_list$data_clean

dim(data)

# -------------------------------------------------------------------------------------------
# ------------------------------------------ EMA --------------------------------------------
# -------------------------------------------------------------------------------------------

# Desire States

data$state_desire <- get_states(data$SI_DesireKill)
data$state_intent <- get_states(data$SI_Intent)


## -------------------------- Markov Chain Modeling -------------------------
## ------- Multi-State Model: Desire
state_tab <- statetable.msm(state_desire, ID, data = data)     ## transition frequencies
state_tab        ## these are at a discrete 1-observation level
sum(state_tab)

trans_mat <- prop.table(state_tab, 1)    ## transition matrix (conditional relative frequencies)
trans_mat %>% round(3)
rowSums(trans_mat)

## --- fit model
## define the structure of the transition intensity matrix Q
# We chose this because low probabilities of transition from state 1 to 3 and back
# Q <- rbind(c(0, 0.1, 0.001), c(.1, 0, .1), c(.001, .1, 0))
Q <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
rownames(Q) <- colnames(Q) <- c("low", "mild",  "severe")
Q       ## which instantaneous (!) transitions are allowed; main diagonals should be 0.

desire_msm <- msm(state_desire ~ RunTime, subject = ID, data = data, qmatrix = Q, control=list(fnscale=4000))
desire_msm                    ## transition intensities (rates)
qmatrix.msm(desire_msm , ci = "none")  %>% round(3)      ## same in matrix format

#sojurn time
## average sojourn time for each state (in hours)
sojourn.msm(desire_msm)
# given a state, what is the probability of the next state?
  # note: low -> severe and severe -> low transitions not possible due to model specification
  # (and seen rarely in the data, as can be seen in state_tab and trans_mat)
pnext.msm(desire_msm)
# interpretation: when mild desire, 2-times more likely to transition to severe than to low

pmat_desire <- pmatrix.msm(desire_msm, t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_desire
# note - pmat shows an exponential decay over t
# all.equal(pmat$estimates %*% pmat$estimates, pmatrix.msm(SI_msm, t = 2, ci = "normal")$estimates)
pmat_desire$estimates %>% round(2) ## rows sum up to 1
# rowSums(pmat_desire$estimates)

## ------- Multi-State Model: Intent --------
state_tab2 <- statetable.msm(state_intent, ID, data = data)     ## transition frequencies
state_tab2        ## these are at a discrete 1-observation level
sum(state_tab2)

trans_mat2 <- prop.table(state_tab2, 1) ## transition matrix (conditional relative frequencies)
trans_mat2 %>% round(3)
rowSums(trans_mat2)

# fit model
intent_msm <- msm(state_intent~ RunTime, subject = ID, data = data, qmatrix = Q, control=list(fnscale=4000))
intent_msm

sojourn.msm(intent_msm)
pnext.msm(intent_msm)
# interpretation: when mild intent, slightly more likely to transition to low than severe

pmat_intent <- pmatrix.msm(intent_msm , t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_intent

# --------------------------------------------------------------------------------------------------
# -------------------------------------  Visualize model estimates ---------------------------------
# --------------------------------------------------------------------------------------------------

# ------------- desire ---------------------------
colmate <- matrix(c("#3288BD"), 3,3,byrow = T)

laymat = matrix(c(1,-1,
                  -1,0,
                  1,1),3,2,byrow =T)
ltymat = matrix(c(1,1,1,
                  2,1,1,
                  2,2,1),3,3,byrow =T)

# pdf(file = "./ctsem_new/figures/msmnet_short_desire.pdf")
pdf(file = "figures/msmnet_desire_full.pdf")
plotnet(pmat_desire$estimates, lower = pmat_desire$L, upper = pmat_desire$U, laymat = laymat,
        digits = 2, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times
sojourn.msm(desire_msm)$estimates
# calculate angles for clock visualization
 360 - (sojourn.msm(desire_msm)$estimates/12)*360
 360 - ((sojourn.msm(desire_msm)$estimates[1]-12)/12)*360


# ------------- intent---------------------------
# pdf(file = "./ctsem_new/figures/msmnet_short_intent.pdf")
pdf(file = "figures/msmnet_intent_full.pdf")
plotnet(pmat_intent$estimates, lower = pmat_intent$L, upper = pmat_intent$U, laymat = laymat,
        digits = 2, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times

sojourn.msm(intent_msm)$estimates
 360 - (sojourn.msm(intent_msm)$estimates/12)*360
 360 - ((sojourn.msm(intent_msm)$estimates[1]-12)/12)*360


 # ---------- Save estimates ---------
qboot_desire  <- qmatrix.msm(desire_msm , ci = "bootstrap", cores = 2)        ## same in matrix format
 round(qboot_desire,3)
 round(sojourn.msm(desire_msm) ,3)
 round(pnext.msm(desire_msm)[2,], 3)[c(1,3),]

 qboot_intent <- qmatrix.msm(intent_msm , ci = "bootstrap", cores = 2)
 round(sojourn.msm(intent_msm) ,3)
 round(pnext.msm(intent_msm)[2,], 3)[c(1,3),]

