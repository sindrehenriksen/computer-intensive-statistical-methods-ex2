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


set.seed(123)
# loading data
R = load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")


input <- list(
  y = Oral$Y, 
  E = Oral$E, 
  n = length(Oral$Y), 
  alpha = 1, 
  beta = 0.01,
  M = 50000,
  R = R
  )

get_b <- function(input,z){
  return(y+E%*%exp(z))
}

get_c <- function(input,z){
  return(input$E%*%exp(z))
}


## ---- 21
# Draw samples from the full condition of kappa_u
r_kappa_u <- function(input){
  shape = (input$n-1)/2 + input$alpha
  rate = 1/2 * (t(u) %*%R ) %*% u + input$beta
  return(rgamma(n = 1, shape = shape, rate = rate))
}

## ---- 22
# Draw samples from the full condition of kappa_v
r_kappa_v <- function(input){
  shape = n/2 + input$alpha
  rate = 1/2 * (t(eta-u)) %*% (eta - u) + input$beta
  return(rgamma(n = 1 ,shape = shape, rate = rate))
}
## ---- 23
# draw samples from the full conditional of u
r_u <- functuion(input,kappa_v,kappa_v){
  inverted.calc <- solve(diag.spam(k_v, input$n) + kappa_u*R)
  return(rmvnorm.canonical(
    n = 1, 
    b = inverted.calc*k_v%*%eta, 
    Q = inverted.calc)
  )
}

## ---- 24
# draw samples from the proposal density of eta
r_eta_prop <- function(input,z,u,kappa_u,kappa_v){
  c_mat = get_c(input,z)
  b_mat = get_b(input,z)
  return(rmvnorm.canonical(
    n = 1,
    b = kappa_v%*%u + b_mat,
    Q = diag.spam(kappa_v, n) + diag.spam(c_mat))
  )
}
## ---- 25
# acceptance of proposal eta




