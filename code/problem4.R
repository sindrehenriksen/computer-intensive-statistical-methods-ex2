#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# Solution to excercise 4
## ---- ess
# effective sample size
library(tidyverse)
library(coda)
library(kableExtra)
load("data/input.Rdata")
load("data/mcmc_df.Rdata")
ess <- effectiveSize(MCMC_list[MCMC_list$is_burnin==F, 2:9])
ess
## ---- break
save(file = "../code/data/ess4.Rdata", ess)
## ---- essPrint
load("../code/data/ess4.Rdata")
cat("ESS[v1] =", ess[1])
cat("ESS[v2] =", ess[2])
cat("ESS[v3] =", ess[3])
cat("ESS[u1] =", ess[4])
cat("ESS[u2] =", ess[5])
cat("ESS[u3] =", ess[6])
cat("ESS[kappa_u] =", ess[7])
cat("ESS[kappa_v] =", ess[8])


## ---- relESS
load("../code/data/samples.Rdata")
time <- as.numeric(samples$run_time[1])
relESS <- mean(c(ess[7],ess[8]))/time
relESS
## ---- break
save(file = "../code/data/ress.Rdata",relESS)
## ---- rESSprint
load("../code/data/ress.Rdata")
cat("Relative effective sample size:", relESS)
