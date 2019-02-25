#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
rm(list = ls())
# Libraries ------------------------------------
library(spam)
library(fields, warn.conflicts = FALSE)
library(INLA)
g <- system.file("demodata/germany.graph", package="INLA")
source("./data/ex2_additionalFiles/dmvnorm.R")
source("./data/ex2_additionalFiles/ess.R")
str(Oral)
attach(Oral)
library(colorspace)
col <- diverge_hcl(8) # blue - red
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)


# loading data
load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")
smoking = read.table("./data/ex2_additionalFiles/smoking.dat")