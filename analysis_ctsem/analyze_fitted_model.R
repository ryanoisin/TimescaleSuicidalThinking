library(rstan)
library(ctsem)

#  obj <- readRDS("shortnew_MCMC.RDS")
# # 
# # sumobj <- summary(obj)
# # saveRDS(sumobj, file="sumobj_short.RDS")
# sumobj <- readRDS("sumobj_short.RDS")
obj <- readRDS("ct_analysis_alldata/full_ml.RDS")

sumobj <- summary(obj)
plot(obj)

plotobj <- ctStanDiscretePars(obj, plot = T)
sumobj
# objhmc <- readRDS("shortnew_HMC.RDS")
# summary(objhmc)
# plot(objhmc)
# 
# objml <- readRDS("shortnew_ML.RDS")
# summary(objml)
# 
# plot(objml)

# they all look the same! MCMC maybe looks a bit more convincingly converged

# Individual specific parameters from MCMC?
indpars <- ctStanSubjectPars(obj, pointest = TRUE)
# [1,2] is the small one == d_l_i_l
# d = 1, i = 2
indlist <- list()
for(i in 1:49){
indlist[[i]] <- matrix(c(indpars[1,i,"drift_D_l"], indpars[1,i,"drift_D_l_I_l"],
         indpars[1,i,"drift_I_l_D_l"], indpars[1,i,"drift_I_l"]),2,2,byrow = T,
         dimnames = list(c("D","I"), c("D","I")))
}
# 
# saveRDS(list(sumobj = obj, indlist = indlist), file = "mcmcres_short.RDS")
# 

## Analyse/Plot
sumobj$popmeans

drift_short <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","Mean"],
                      2,2,byrow = T)
drift_short_l <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","2.5%"],
                        2,2,byrow = T)
drift_short_u <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","97.5%"],
                        2,2,byrow = T)

source("../ctsem_new/plotnet.R")
library(qgraph)

# convert to "1-hour" parameters ? 

library(RColorBrewer)
brewer.pal(2,name = "Spectral")
colmat <- matrix(c("#D53E4F","#3288BD",
         "grey", "#D53E4F"), 2,2,byrow = T)

pdf("figures/ctpopnet_short.pdf")
plotnet(t(drift_short), t(drift_short_l), t(drift_short_u), edge.color = colmat, 
        digits = 3, labels = c("Desire","Intent"), label.scale.equal = TRUE, curve = 2,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .9, asize = 5,
        maximum = 2.5)
dev.off()

# plotnet(t(indlist[[5]]))


# ----------------- load EMA results -------------------------

objema <- readRDS("emanew_MCMC.RDS")

sumobje <- summary(objema)
plot(objema)


tmp <- ctStanDiscretePars(objema, plot = T, ggcode = T)


plotobje <- ctStanDiscretePars(objema, plot = T)
pdf("figures/popexp_ema.pdf")
plotobje
# indpars <- ctStanSubjectPars(fit = objema, pointest = FALSE)
# stored here (function doesn't work properly)
# save(out, file = "indpars_ema.RDS")

drift_shorte <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","Mean"],
                      2,2,byrow = T)
drift_short_le <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","2.5%"],
                        2,2,byrow = T)
drift_short_ue <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","97.5%"],
                        2,2,byrow = T)

source("plotnet.R")
library(qgraph)

# convert to "1-hour" parameters ? 
colmate <- matrix(c("#D53E4F","#3288BD",
                    "#D53E4F", "#D53E4F"), 2,2,byrow = T)
pdf("figures/ctpopnet_ema.pdf")
plotnet(t(drift_shorte), t(drift_short_le), t(drift_short_ue), edge.color = colmate, 
        digits = 3, labels = c("Desire","Intent"), label.scale.equal = TRUE, curve = 2,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .9, asize = 5,
        maximum = 2.5)
dev.off()

# 
# 
# rstan::stan_dens(obj$stanfit$stanfit, pars = c("pop_DRIFT"))
# rstan::stan_rhat(stanfit$stanfit)
# # rstan::stan_diag(stanfit$stanfit)
# 
# # plot trace of all possible population parametrs
# substrings = c("pop_") 
# out <- c()
# for (subsi in substrings) {
#   out <- c(out, 
#            stanfit$stanfit@model_pars[grep(paste0("^", subsi),
#                                            stanfit$stanfit@model_pars)])
#  }
# rstan::stan_trace(stanfit$stanfit, pars=out)
# 
# # plot trace only of free parmeters
# out2 <- c("pop_DRIFT","pop_DIFFUSION", "pop_T0VAR[1,2]")
# 
# rstan::stan_trace(stanfit$stanfit, pars = out2)
# 
# pairs(stanfit$stanfit, pars = c("lp__", "pop_DRIFT"))
