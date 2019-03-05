#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# Solution to excercise 4
# effective sample size
library(tidyverse)
library(spam)
library(ggpubr)
load("data/input.Rdata")
load("data/samples.Rdata")
ess <- effectiveSize(data.frame(samples$kappa_u,samples$kappa_v))
