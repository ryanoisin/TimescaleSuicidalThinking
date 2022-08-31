
# file to do initial analysis of ctsem fitted model objects

library(rstan)
library(ctsem)
library(qgraph)
library(RColorBrewer)
source("aux_functions.R")
# I use plotnet() below from aux_functions

# ----------------- load "short" results results -------------------------

obj <- readRDS("files/shortnew_MCMC.RDS")
# Create summary object (can take some time)
# sumobj <- summary(obj)
# saveRDS(sumobj, file="files/sumobj_short.RDS")
sumobj <- readRDS("files/sumobj_short.RDS")
# plot(obj)
# can be used to inspect model parameters and trace plots

# plotobj <- ctStanDiscretePars(obj, plot = T)
# sumobj


# Individual specific parameters from MCMC?
# indpars <- ctStanSubjectPars(obj, pointest = TRUE)
# # [1,2] is the small one == d_l_i_l
# # d = 1, i = 2
# indlist <- list()
# for(i in 1:49){
#   indlist[[i]] <- matrix(c(indpars[1,i,"drift_D_l"], indpars[1,i,"drift_D_l_I_l"],
#                            indpars[1,i,"drift_I_l_D_l"], indpars[1,i,"drift_I_l"]),2,2,byrow = T,
#                          dimnames = list(c("D","I"), c("D","I")))
# }


## Analyse/Plot
sumobj$popmeans

drift_short <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","Mean"],
                      2,2,byrow = T)
drift_short_l <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","2.5%"],
                        2,2,byrow = T)
drift_short_u <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","97.5%"],
                        2,2,byrow = T)



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



# ----------------- load EMA results -------------------------

objema <- readRDS("files/emanew_MCMC.RDS")

# sumobje <- summary(objema)
#  saveRDS(sumobje, file="files/sumobj_ema.RDS")
sumobje <- readRDS("files/sumobj_ema.RDS")

# plot(objema)
# tmp <- ctStanDiscretePars(objema, plot = T, ggcode = T)


# plotobje <- ctStanDiscretePars(objema, plot = T)

# indpars <- ctStanSubjectPars(fit = objema, pointest = FALSE)
# stored here (function doesn't work properly)
# save(out, file = "indpars_ema.RDS")

drift_shorte <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","Mean"],
                       2,2,byrow = T)
drift_short_le <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","2.5%"],
                         2,2,byrow = T)
drift_short_ue <- matrix(sumobje$parmatrices[sumobje$parmatrices$matrix == "DRIFT","97.5%"],
                         2,2,byrow = T)
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

