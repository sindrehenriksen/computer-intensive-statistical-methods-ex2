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
ess <- effectiveSize(data.frame(MCMC_list$kappa_u[MCMC_list$is_burnin==F],
                                MCMC_list$kappa_v[MCMC_list$is_burnin==F]))
ess
## ---- break
save(file = "data/ess4.Rdata", ess)
## ---- essPrint
load("data/ess4.Rdata")
cat("ESS[kappa_u] =", ess[1])
cat("ESS[kappa_v] =", ess[2])
