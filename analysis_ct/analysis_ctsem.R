# file to do analysis of ctsem fitted model object

library(rstan)
library(ctsem)
library(qgraph)
library(RColorBrewer)
library(ctnet)
source("aux_functions.R")
# I use plotnet() and IRFfun() below from aux_functions.R

# --------------------------------------------------------
# ----------------- load results -------------------------
# --------------------------------------------------------

obj <- readRDS("files/full_MCMC_std.RDS")
# check trace plots and posterior densities
# plot(obj)
# Create summary object (can take some time)
# sumobj <- summary(obj)
# Save summary object
# saveRDS(sumobj, file="files/sumobj_full.RDS")
sumobj <- readRDS("files/sumobj_full.RDS")


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

# --------------------------------------------------------
# ---- inspect and plot drift matrix parameters ----------
# --------------------------------------------------------

sumobj$popmeans

drift <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","Mean"],
                      2,2,byrow = T)
drift_l <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","2.5%"],
                        2,2,byrow = T)
drift_u <- matrix(sumobj$parmatrices[sumobj$parmatrices$matrix == "DRIFT","97.5%"],
                        2,2,byrow = T)


colmat <- matrix(c("#D53E4F","#3288BD",
                   "grey", "#D53E4F"), 2,2,byrow = T)

pdf("figures/ctpopnet_full.pdf")
plotnet(t(drift), t(drift_l), t(drift_u), edge.color = colmat,
        digits = 3, labels = c("Desire","Intent"), label.scale.equal = TRUE, curve = 2,
        edge.label.cex = 1.1,
        vsize = 12, rescale = F, layscale = .9, asize = 5,
        maximum = 2.5)
dev.off()

# --------------------------------------------------------
# --------- Create lagged effects and IRF plots ----------
# --------------------------------------------------------

# set up colors
colvec <- brewer.pal(4, "Spectral")
colirfs <- brewer.pal(11,"Spectral")[c(1,11)]

# extract posterior samples of the drift matrix
post <- ctsem::ctExtract(obj)
drift_post <- post$pop_DRIFT
dts <- seq(0,24,.1)

# use the ctnet package function getCIs (can take some time)
# here we obtain CIs for the "lagged effects"  at different time-intervals
# phidt_CI <- sapply(dts,function(dt){
#   getCIs(drift_post,simplify=TRUE, FUN=expm::expm, const = dt)
# }, simplify = "array")

# save this object because it takes some time to run
# saveRDS(phidt_CI,"files/phidt_CI_full.RDS")
 phidt_CI <- readRDS("files/phidt_CI_full.RDS")

# here we use the same procedure to obtain CIs for the IRF at different times
# IRFfun is defined in aux_functions.R
# IRF_CI <- sapply(dts,function(dt){
  # getCIs(drift_post,simplify=TRUE, FUN=IRFfun, const = dt, start = c(1,0))
# }, simplify = "array")

# save this object because it takes some time to run
 # saveRDS(IRF_CI,"files/IRF_CI_full.RDS")
 IRF_CI <- readRDS("files/IRF_CI_full.RDS")


 # --------- Make lagged effects plot ----------

pdf("figures/popexp_full.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plotPhi(CI_obj = phidt_CI, dts = dts,  index = "all", leg = TRUE, colvec = colvec,poly =TRUE,
 drawaxis = F, main = "", xlab, cex.lab = 1.5)
axis(1,at = c(0,seq(4,24,4)))
axis(2)
dev.off()

# ------------- Make IRF plot ----------------

# set-up
ind <- c(1,2)
ylim = c(min(IRF_CI[, ind, ]), max(IRF_CI[, ind, ]))
colvec2 <-brewer.pal(11,"Spectral")[c(1,11)]

pdf("figures/popIRF_full.pdf")
poly = TRUE # optional: plot CIs using polygon function

# set up plotting window
par(cex.axis=1.5, cex.lab=1.5)
plot.new()
plot.window(xlim = c(dts[1], max(dts)), ylim = c(0,1))

# plot the IRF and it's CIs
for (i in 1:length(ind)) {
  lines(dts, IRF_CI[2, ind[i], ], col = colvec2[i])
  lines(dts, IRF_CI[1, ind[i], ], col = colvec2[i], lty = 2)
  lines(dts, IRF_CI[3, ind[i], ], col = colvec2[i], lty = 2)
  if(isTRUE(poly)){
    polygon(x = c(dts, rev(dts)),
            y = c( IRF_CI[1, ind[i], ], rev( IRF_CI[3, ind[i], ])),
            col = alpha(colvec2[i], .1), lty = 0)}
}

# add points
abline(h = 0, col = "grey")
points(dts[1], IRF_CI[2,1,1], pch = 18, cex = 3, col = colvec2[1])
points(dts[1], IRF_CI[2,2,1], pch = 18, cex = 3, col = colvec2[2])

axis(1,at = c(0,seq(4,24,4)))
axis(2)
title(main = "",
      xlab = "Time",
      ylab = "Process Value")
legend("topright", col = c(colvec2,"gray"), lty = c(1,1,NA), pch = c(NA,NA,18),
       legend = c("Desire","Intent","Impulse Value"),xpd = T, bty = "n",
       cex = 1.5, pt.cex = 3)
dev.off()


