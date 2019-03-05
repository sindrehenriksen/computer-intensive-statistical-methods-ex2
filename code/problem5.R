#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(tidyverse)
library(spam)
library(ggpubr)
library(colorspace)
load("data/input.Rdata")
load("data/samples.Rdata")
burnin = 5000
M = length(samples$u[,1])
u_med_exp <- function(u){
  temp_median = apply(u,2, median, na.rm = TRUE)
  posterior = exp(temp_median)
}
col <- diverge_hcl(8) # blue - red
pdf(file = "../figures/germany_data.pdf")
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)
dev.off()
pdf(file = "../figures/germany_post.pdf")
germany.plot(u_med_exp(samples$u[seq(burnin,M),]), col = col, legend=TRUE)
dev.off()
