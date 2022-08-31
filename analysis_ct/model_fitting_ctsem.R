# code to fit ctsem model to full dataset
  # Note: due to long run-time this model was fit on a different machine than the one
  # used to run the analyses. For details see files/sessionInfo_ctsem.txt

library(ctsem)
library(rstan)
set.seed(1234)

# ----------------------------------------------------------------------
# --------- Load and Prepare Data  -------------------------------------
# ----------------------------------------------------------------------

data_list <- readRDS(file = "data/data_full_ctsem.RDS")
data <- data_list$data_std

colnames(data)[4] <- "TI"

# ---------- CT-VAR Model Estimation -------------

# make variable labels
lab <- c("D","I")
p <- 2


# Specify simpler CT-VAR model
ctmodel <- ctModel(type='stanct',
                   manifestNames=c("SI_DesireKill","SI_Intent"),
                   latentNames= paste0(lab,"_l"),
                   LAMBDA = diag(nrow=p),
                   DRIFT = "auto",
                   # MANIFESTMEANS = matrix(data=0, nrow=p, ncol=1),
                   # MANIFESTVAR=diag(0,p),
                   CINT = matrix(data = 0, nrow = p, ncol = 1),
                   DIFFUSION ="auto",
                   # T0MEANS=matrix(data = 0, nrow = p, ncol = 1),
                   # T0VAR=matrix(data = c(1,"t0cov_est","t0cov_est",1), nrow = p, ncol = p),
                   time = "TI",
                   id = "ID")

# let drift and diffusion parameters vary across individuals -- drift only is heavily influenced by unmodeled diffusion heterogeneity
ctmodel$pars$indvarying[7:14] <- TRUE

# For bayesian MCMC (slow, more reliable/stable), optimize = FALSE, noprriors = FALSE
chains = 3
cores =3
iter = 5000
ctfit <- ctStanFit(data, ctmodel,
                   nlcontrol = list(maxtimestep = 100),
                   optimize = FALSE,
                   nopriors = FALSE,
                   chains = chains, cores = cores,plot=10,
                   iter = iter, fit = TRUE, verbose = 0)
saveRDS(ctfit, "full_MCMC_std.RDS")

sink("sessionInfo_ctsem.txt")
sessionInfo()
sink()
