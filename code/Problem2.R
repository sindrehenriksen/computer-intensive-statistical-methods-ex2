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
library(here)
here()
col <- diverge_hcl(8) # blue - red
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)


set.seed(123)
# loading data
load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")


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
  return(input$y+E*exp(z$sample))
}

get_c <- function(input,z){
  return(input$E*exp(z$sample))
}

d_eta <- function(input,eta,kappa_v,u){
  return(-1/2*t(eta)%*%(diag.spam(kappa_v,input$n)%*%eta) + 
           t(eta)%*%(kappa_v*u) + t(eta)%*%input$y-t(exp(eta))%*%E)
}

acceptance_prob <- function(input,eta_prop,eta,kappa_v,u){
  return(min(1,exp(d_eta(input,eta_prop,kappa_v,u) - d_eta(input,eta,kappa_v,u) +
               - eta_prop$prob + eta$porb)))
}


## ---- 21
# Draw samples from the full condition of kappa_u
r_kappa_u <- function(input,u){
  shape = (input$n-1)/2 + input$alpha
  rate = 0.5 * t(u) %*% input$R %*% u + input$beta
  return(rgamma(n = 1, shape = shape, rate = rate))
}

## ---- 22
# Draw samples from the full condition of kappa_v
r_kappa_v <- function(input,eta,u){
  shape = input$n/2 + input$alpha
  rate = 1/2 * (t(eta$sample-u)) %*% (eta$sample - u) + input$beta
  return(rgamma(n = 1 ,shape = shape, rate = rate))
}
## ---- 23
# draw samples from the full conditional of u
r_u <- function(input,kappa_u,kappa_v){
  computation<- diag.spam(kappa_v, input$n) + kappa_u*input$R
  u = rmvnorm.canonical(n = 1,  b = kappa_v*computation%*%eta$sample, Q = computation)
  return(c(u))
}

## ---- 24
# draw samples from the proposal density of eta
r_eta_prop <- function(input,z,u,kappa_v){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  print(b_vec)
  b = c(kappa_v*u + b_vec)
  Q = diag.spam(kappa_v, input$n) + diag.spam(c_vec)
  print(b[1])
  print(Q)
  sample = rmvnorm.canonical(n = 1, b = b, Q = Q)
  prob = dmvnorm.canonical(x = sample, b = b, Q = Q, log = TRUE)
  return(list(sample=sample,prob=prob))
}
## ---- 25
# acceptance of proposal eta
u = c(rep_len(0.0,input$n))
eta <- list(sample = c(rep_len(0.0,input$n)), prob = 1)
i = 0
while (i < input$M){
  kappa_u = r_kappa_u(input,u)
  kappa_v = r_kappa_v(input,eta,u)
  u = r_u(input,kappa_u,kappa_v)
  eta_prop = r_eta_prop(input,eta,u,kappa_v)
  accept_prob <- acceptance_prob(input,eta_prop,eta,kappa_v,u)
  if(ruinf(1) < accept_prob){
    eta = eta_prop
    i = i + 1
  }
}


