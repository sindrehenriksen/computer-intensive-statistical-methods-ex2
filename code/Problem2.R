#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# Use log densities

# Libraries ------------------------------------
library(spam)
library(fields, warn.conflicts = FALSE)
source("./data/ex2_additionalFiles/dmvnorm.R")
source("./data/ex2_additionalFiles/ess.R")
str(Oral)
attach(Oral)
library(colorspace)
col <- diverge_hcl(8) # blue - red
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)


# loading data
load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")
M <- 50000

## ---- 21


## ---- 22


## ---- 23


## ---- 24
d_eta <- function(z,y,u,kappa_u,kappa_v){
  
}
## ---- 25