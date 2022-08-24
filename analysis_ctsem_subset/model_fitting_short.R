# code to fit the CT-VAR(1) model to the "short" dataset

library(ctsem)
library(rstan)
set.seed(1234)

# ----------------------------------------------------------------------
# --------- Load and Prepare Data  -------------------------------------
# ----------------------------------------------------------------------

data_ema_gonly <- readRDS(file="Data/data_short_50.RDS")

ids <- unique(data_ema_gonly$ID)

data_ema_gonly <- data_ema_gonly[,c("ID","SI_DesireKill", "SI_Intent","RunTime")]
colnames(data_ema_gonly)[4] <- "TI"

for(i in 1:length(ids)){
  
  ema <- data_ema_gonly[data_ema_gonly$ID == ids[i],]
  ema[,c("SI_DesireKill","SI_Intent")] <- apply(ema[,c("SI_DesireKill","SI_Intent")],2, scale)
  ema[,"ID"] <- i
  data_ema_gonly[data_ema_gonly$ID == ids[i],] <- ema
}


# ---------- CT-VAR Model Estimation -------------

# make variable labels
lab <- c("D","I")
p <- 2

# Specify simple CT-VAR model
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

# For bayesian MCMC (slow, more reliable/stable)
chains = 3
cores =3
iter = 5000
ctfit_ema <- ctStanFit(data_ema_gonly, ctmodel,
                       nlcontrol = list(maxtimestep = 100),
                       optimize = FALSE,
                       nopriors = FALSE,
                       chains = chains, cores = cores,plot=10,
                       iter = iter, fit = TRUE, verbose = 1)
saveRDS(ctfit_ema, "files/shortnew_MCMC.RDS")
