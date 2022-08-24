# file to plot "curves" and IRFs

library(ctnet)
library(RColorBrewer)
colvec <- brewer.pal(4, "Spectral")
colirfs <- brewer.pal(11,"Spectral")[c(1,11)]
# IRF function
IRFfun <- function(drift, start){
  expm::expm(drift)%*%start
}

 # obj <- readRDS("shortnew_MCMC.RDS")
obj <- readRDS("ct_analysis_alldata/full_ml.RDS")

post_s <- ctsem::ctExtract(obj)
s_drift <- post_s$pop_DRIFT
dts <- seq(0,24,.1)

# use the ctnet package function getCIs
phidt_CI <- sapply(dts,function(dt){
  getCIs(s_drift,simplify=TRUE, FUN=expm::expm, const = dt)
}, simplify = "array")

# use the 

pdf("figures/popexp_short.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plotPhi(CI_obj = phidt_CI, dts = dts,  index = "all", leg = TRUE, colvec = colvec)
        # drawaxis = F, main = "", xlab, cex.lab = 1.5)
axis(1,at = c(0,seq(4,24,4)))
axis(2)
dev.off()
# IRFs

IRF_s <- sapply(dts,function(dt){
  getCIs(s_drift,simplify=TRUE, FUN=IRFfun, const = dt, start = c(1,0))
}, simplify = "array")

# see if this works
# plot(IRF_s[1,1,])
# plotPhi(CI_obj = IRF_s, dts = dts, index = "all", leg = FALSE, colvec = colvec[c(1,4)], 
#         drawaxis = F, main = "", xlab= "", ylab = "", ind = c(1,2))

# manually
ind <- c(1,2)
ylim = c(min(IRF_s[, ind, ]), max(IRF_s[, ind, ]))
# colvec2 <- colvec[2:3]
 colvec2 <-brewer.pal(11,"Spectral")[c(1,11)]
 
pdf("figures/popIRF_short.pdf")
poly = TRUE
par(cex.axis=1.5, cex.lab=1.5)
plot.new()
plot.window(xlim = c(dts[1], max(dts)), ylim = c(0,1))
abline(h = 0, col = "grey")

for (i in 1:length(ind)) {
  lines(dts, IRF_s[2, ind[i], ], col = colvec2[i])
  lines(dts, IRF_s[1, ind[i], ], col = colvec2[i], lty = 2)
  lines(dts, IRF_s[3, ind[i], ], col = colvec2[i], lty = 2)
  if(isTRUE(poly)){
    polygon(x = c(dts, rev(dts)), 
            y = c( IRF_s[1, ind[i], ], rev( IRF_s[3, ind[i], ])),
            col = alpha(colvec2[i], .1), lty = 0)}
}
# add points
points(dts[1], IRF_s[2,1,1], pch = 18, cex = 3, col = colvec2[1])
points(dts[1], IRF_s[2,2,1], pch = 18, cex = 3, col = colvec2[2])
axis(1,at = c(0,seq(4,24,4)))
axis(2)
title(main = "",
      xlab = "Time",
      ylab = "Process Value")
legend("topright", col = c(colvec2,"gray"), lty = c(1,1,NA), pch = c(NA,NA,18),
       legend = c("Desire","Intent","Impulse Value"),xpd = T, bty = "n",
       cex = 1.5, pt.cex = 3)
dev.off()

# -------------------- Repeat for EMA --------------------------------------


 obje <- readRDS("emanew_MCMC.RDS")


post_e <- ctsem::ctExtract(obje)
e_drift <- post_e$pop_DRIFT
dts <- seq(0,24,.1)

phidt_CIe <- sapply(dts,function(dt){
  getCIs(e_drift,simplify=TRUE, FUN=expm::expm, const = dt)
}, simplify = "array")

pdf("figures/popexp_ema.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plotPhi(CI_obj = phidt_CIe, dts = dts,  index = "all", leg = TRUE, colvec = colvec, 
        drawaxis = F, main = "", xlab, cex.lab = 1.5)
axis(1,at = c(0,seq(4,24,4)))
axis(2)
dev.off()
# IRFs

IRF_e <- sapply(dts,function(dt){
  getCIs(e_drift,simplify=TRUE, FUN=IRFfun, const = dt, start = c(1,0))
}, simplify = "array")

# see if this works
# plot(IRF_s[1,1,])
# plotPhi(CI_obj = IRF_s, dts = dts, index = "all", leg = FALSE, colvec = colvec[c(1,4)], 
#         drawaxis = F, main = "", xlab= "", ylab = "", ind = c(1,2))

# manually
ind <- c(1,2)
ylim = c(min(IRF_e[, ind, ]), max(IRF_e[, ind, ]))

pdf("figures/popIRF_ema.pdf")
par(cex.axis=1.5, cex.lab=1.5)
poly = TRUE
plot.new()
plot.window(xlim = c(dts[1], max(dts)), ylim = c(0,1))
abline(h = 0, col = "grey")

for (i in 1:length(ind)) {
  lines(dts, IRF_e[2, ind[i], ], col = colvec2[i])
  lines(dts, IRF_e[1, ind[i], ], col = colvec2[i], lty = 2)
  lines(dts, IRF_e[3, ind[i], ], col = colvec2[i], lty = 2)
  if(isTRUE(poly)){
    polygon(x = c(dts, rev(dts)), 
            y = c( IRF_e[1, ind[i], ], rev( IRF_e[3, ind[i], ])),
            col = alpha(colvec2[i], .1), lty = 0)}
}
# add points
points(dts[1], IRF_e[2,1,1], pch = 18, cex = 3, col = colvec2[1])
points(dts[1], IRF_e[2,2,1], pch = 18, cex = 3, col = colvec2[2])
axis(1,at = c(0,seq(4,24,4)))
axis(2)
title(main = "",
      xlab = "Time",
      ylab = "Process Value")
legend("topright", col = c(colvec2,"gray"), lty = c(1,1,NA), pch = c(NA,NA,18),
       legend = c("Desire","Intent","Impulse Value"),xpd = T, bty = "n",
       cex = 1.5, pt.cex = 3)
dev.off()
# 
# saveRDS(list(s_drift = s_drift, phidt_CIs = phidt_CI, IRF_s = IRF_s,
#              e_drift = e_drift, phidt_CIe = phidt_CIe, IRF_e = IRF_e), file= "pop_post_helpers.RDS")

