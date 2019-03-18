#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- 5
library(spam)
library(colorspace)
load("data/input.Rdata")
load("data/samples.Rdata")
burnin = 5000
M = length(samples$u[,1])
u_med_exp <- function(u){
  temp_median = apply(u,2, median, na.rm = TRUE)
  posterior = exp(temp_median)
}
# color scheme
col <- diverge_hcl(8) # blue - red
pdf(file = "../figures/germany_5a.pdf", width=10, height=8)
par(mfrow = c(1,2))
# standardised mortality rates
germany.plot(Oral$Y/Oral$E, 
             col=col, 
             legend=TRUE,
             main="Standardised mortaility rates",
             cex.main=1)
# spatial structured effects
germany.plot(u_med_exp(samples$u[seq(burnin,M),]), 
             col = col, 
             legend=TRUE,
             main="Spatial structured effects",
             cex.main=1)

## ---- break
dev.off()
