# code to run the markov switching model analyses in the main text

#load packages
library(msm)
library(tidyverse)
library(DataCombine)


# function to re-code continuous variables into states
get_states <- function(vector,th1 = 4 ,th2 = 7){
  out <- rep(NA,length(vector))
  out[vector < th1 | vector ==th1] <- 1
  out[vector > th1 & vector < th2 | vector ==th2 ] <- 2
  out[vector > th2] <- 3
  out
}


#### load data for full data analysis
#TBD

#### load data for subset analysis

ema <- readRDS("Data/data_ema_50.RDS")
short <-  readRDS("Data/data_short_50.RDS")
### create SI states

length(unique(ema$ID)); length(unique(short$ID))

dim(ema); dim(short)

# -------------------------------------------------------------------------------------------
# ------------------------------------------ EMA --------------------------------------------
# -------------------------------------------------------------------------------------------

# Desire States

ema$state_desire <- get_states(ema$SI_DesireKill)
ema$state_intent <- get_states(ema$SI_Intent)


## -------------------------- Markov Chain Modeling -------------------------
## ------- Multi-State Model: Desire
state_tab <- statetable.msm(state_desire, ID, data = ema)     ## transition frequencies
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

ema_desire_msm <- msm(state_desire ~ RunTime, subject = ID, data = ema, qmatrix = Q, control=list(fnscale=4000))
ema_desire_msm                    ## transition intensities (rates)
qmatrix.msm(ema_desire_msm , ci = "none")  %>% round(3)      ## same in matrix format

#sojurn time
## average sojourn time for each state (in hours)
sojourn.msm(ema_desire_msm )
pnext.msm(ema_desire_msm )


pmat_ema_desire <- pmatrix.msm(ema_desire_msm, t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_ema_desire
# note - pmat shows an exponential decay over t
# all.equal(pmat$estimates %*% pmat$estimates, pmatrix.msm(SI_msm, t = 2, ci = "normal")$estimates)
pmat$estimates %>% round(2) ## rows sum up to 1
rowSums(pmat$estimates)

## ------- Multi-State Model: Intent
state_tab2 <- statetable.msm(state_intent, ID, data = ema)     ## transition frequencies
state_tab2        ## these are at a discrete 1-observation level
sum(state_tab2)   

trans_mat2 <- prop.table(state_tab2, 1)    ## transition matrix (conditional relative frequencies)
trans_mat2 %>% round(3)
rowSums(trans_mat2)

# fit model
ema_intent_msm <- msm(state_intent~ RunTime, subject = ID, data = ema, qmatrix = Q, control=list(fnscale=4000))
ema_intent_msm

sojourn.msm(ema_intent_msm)
pnext.msm(ema_intent_msm)

pmat_ema_intent <- pmatrix.msm(ema_intent_msm , t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_ema_intent



# -------------------------------------------------------------------------------------------
# ------------------------------------------ Short --------------------------------------------
# -------------------------------------------------------------------------------------------
short$state_desire <- get_states(short$SI_DesireKill)
short$state_intent <- get_states(short$SI_Intent)


# ------------- Desire -----------------
short_desire_msm <- msm(state_desire ~ RunTime, subject = ID, data = short, qmatrix = Q, control=list(fnscale=4000))
short_desire_msm                    ## transition intensities (rates)
qmatrix.msm(short_desire_msm , ci = "none")  %>% round(3)      ## same in matrix format

#sojurn time
## average sojourn time for each state (in hours)
sojourn.msm(short_desire_msm )
pnext.msm(short_desire_msm )


pmat_short_desire <- pmatrix.msm(short_desire_msm, t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_short_desire

# ------------- Intent -----------------

short_intent_msm <- msm(state_intent~ RunTime, subject = ID, data = short, qmatrix = Q, control=list(fnscale=4000))
short_intent_msm                    ## transition intensities (rates)
qmatrix.msm(short_intent_msm , ci = "none")  %>% round(3)      ## same in matrix format

#sojurn time
## average sojourn time for each state (in hours)
sojourn.msm(short_intent_msm )
pnext.msm(short_intent_msm )


pmat_short_intent <- pmatrix.msm(short_intent_msm, t = 1, ci = "normal")   ## t = 1 ... time interval to estimate transition probabilities
pmat_short_desire

# --------------------------------------------------------------------------------------------------
# Visualizations

# source("./ctsem_new/plotnet.R")
source("aux_functions.R")
library(qgraph)

# desire 
colmate <- matrix(c("#3288BD"), 3,3,byrow = T)

laymat = matrix(c(1,-1,
                  -1,0,
                  1,1),3,2,byrow =T)
ltymat = matrix(c(1,1,1,
                  2,1,1,
                  2,2,1),3,3,byrow =T)

# pdf(file = "./ctsem_new/figures/msmnet_short_desire.pdf")
pdf(file = "figures/msmnet_short_desire.pdf")
plotnet(pmat_short_desire$estimates, lower = pmat_short_desire$L, upper = pmat_short_desire$U, laymat = laymat,
        digits = 2, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times
sojourn.msm(short_desire_msm)$estimates
360 - (sojourn.msm(short_desire_msm)$estimates/12)*360

# pdf(file = "./ctsem_new/figures/msmnet_ema_desire.pdf")
pdf(file = "figures/msmnet_ema_desire.pdf")
plotnet(pmat_ema_desire$estimates, lower = pmat_ema_desire$L, upper = pmat_ema_desire$U, laymat = laymat,
        digits = 3, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times
sojourn.msm(ema_desire_msm)$estimates

360 - (c(sojourn.msm(ema_desire_msm)$estimates[1] - 36, sojourn.msm(ema_desire_msm)$estimates[2:3])/12)*360


# intent#

# pdf(file = "./ctsem_new/figures/msmnet_short_intent.pdf")
pdf(file = "figures/msmnet_short_intent.pdf")
plotnet(pmat_short_intent$estimates, lower = pmat_short_intent$L, upper = pmat_short_intent$U, laymat = laymat,
        digits = 2, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times

sojourn.msm(short_intent_msm)$estimates
360 - (sojourn.msm(short_intent_msm)$estimates/12)*360

# pdf(file = "./ctsem_new/figures/msmnet_ema_intent.pdf")
pdf(file = "figures/msmnet_ema_intent.pdf")
plotnet(pmat_ema_intent$estimates, lower = pmat_ema_intent$L, upper = pmat_ema_intent$U, laymat =laymat,
        digits = 3, labels = c("Low","Mild","Severe"), label.scale.equal = TRUE, curve = 2, edge.color = colmate,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .8, asize = 5,
        maximum = 1,
        lty = ltymat)
dev.off()

# sojourn times
sojourn.msm(ema_intent_msm)$estimates

360 - (c(sojourn.msm(ema_intent_msm)$estimates[1] - 48, sojourn.msm(ema_intent_msm)$estimates[2:3])/12)*360


