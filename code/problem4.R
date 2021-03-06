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
kable(ess,caption = "\\label{tab:ess}Effective sample size
      of the parameters",booktabs = T) 
## ---- break
save(file = "../code/data/ess4.Rdata", ess)
## ---- essPrint
load("../code/data/ess4.Rdata")
kable(ess,caption = "\\label{tab:ess}Effective sample size
      of the parameters",booktabs = T) 

## ---- relESS
load("../code/data/samples.Rdata")
time <- as.numeric(samples$run_time[1])
relESS <- mean(c(ess[7],ess[8]))/time
cat("Relative effective sample size:", relESS)
## ---- break
save(file = "../code/data/ress.Rdata",relESS)
## ---- rESSprint
load("../code/data/ress.Rdata")
cat("Relative effective sample size:", relESS)
