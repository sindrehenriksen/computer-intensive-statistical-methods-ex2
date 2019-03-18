#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- libs2
library(spam)
library(tidyverse)
library(fields, warn.conflicts = FALSE)
source("./data/ex2_additionalFiles/dmvnorm.R")
library(colorspace)
library(ggpubr)

## ---- data2
load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")
str(Oral)
attach(Oral)

# a list of all the input variables to make code more readable
input <- list(
  y = Oral$Y,
  E = Oral$E,
  n = length(Oral$Y),
  alpha = 1,
  beta = 0.01,
  R = R
  )

## ---- functions
## ---- get_b
# b values
get_b <- function(input,z){
  return(input$y+input$E*exp(z)*(z-1))
}

## ---- get_c
# c values
get_c <- function(input,z){
  return(input$E*exp(z))
}

## ---- r_kappa_u
# Draw samples from the full condition of kappa_u
r_kappa_u <- function(input,u){
  shape = (input$n-1)/2 + input$alpha
  rate = 1/2* t(u) %*% input$R %*% u + input$beta
  return(rgamma(n = 1, shape = shape, rate = rate))
}

## ---- r_kappa_v
# Draw samples from the full condition of kappa_v
r_kappa_v <- function(input,eta,u){
  shape = input$n/2 + input$alpha
  rate = 1/2 * t(eta - u) %*% (eta - u) + input$beta
  return(rgamma(n = 1 ,shape = shape, rate = rate))
}

## ---- r_u
# draw samples from the full conditional of u
r_u <- function(input,kappa_u,kappa_v,eta){
  Q = diag.spam(kappa_v, input$n) + kappa_u*input$R
  b = kappa_v * eta
  return(c(rmvnorm.canonical(n = 1,  b = b, Q = Q)))
}

## ---- r_eta_prop
# draw samples from the proposal density of eta
r_eta_prop <- function(input,z,u,kappa_v){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  b = kappa_v*u + b_vec
  Q = diag.spam(kappa_v, input$n) + diag.spam(c_vec)
  return(c(rmvnorm.canonical(n = 1, b = b, Q = Q)))
}

## ---- d_eta_p
# finding the p-value of eta from the
# full condition p(eta|...)
d_eta_p <- function(input,eta,kappa_v,u){
  return(-1/2*t(eta)%*%diag.spam(kappa_v,input$n)%*%eta +
           t(eta)%*%(kappa_v*u) +
           t(eta)%*%input$y -
           t(exp(eta))%*%input$E)
}

## ---- d_eta_q
# finding the p-value of eta from the
# proposal density q(eta|...)
d_eta_q <- function(input,eta,z,kappa_v,u){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  b = kappa_v*u + b_vec
  Q = diag.spam(kappa_v,input$n) + diag.spam(c_vec)
  return(dmvnorm.canonical(x = eta, b = b, Q = Q, log = TRUE))
}

## ---- acceptance_prob
# calculating the acceptance probability
acceptance_prob <- function(input,eta_prop,eta,kappa_v,u){
  return(min(1,exp(
      d_eta_p(input, eta_prop, kappa_v, u) +
      d_eta_q(input, eta, eta_prop, kappa_v, u) -
      d_eta_p(input, eta, kappa_v, u) -
      d_eta_q(input, eta_prop, eta, kappa_v, u))))
}

## ---- myMCMC
# running MCMC simulation
# return a list of eta, u, kappa_u, kappa_v, and the acceptance
M <- 70000
set.seed(123)
myMCMC <- function(input, M){
  # for keeping track of completition
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  # chosen u
  u = runif(input$n)
  # eta from proposal density
  eta = r_eta_prop(input,u,u,kappa_v =180)
  # storing all etas,vs, kappas and the acceptance prob
  eta_samples = matrix(NA,nrow=M,ncol=input$n)
  u_samples = matrix(NA,nrow=M,ncol=input$n)
  kappa_u_samples = vector()
  kappa_v_samples = vector()
  accept_vec = vector()
  for (i in seq(1,M)){
    # for keeping track of completition
    setTxtProgressBar(pb, i)
    # drawing kappas
    kappa_u = r_kappa_u(input,u)
    kappa_v = r_kappa_v(input,eta,u)
    # drawin u
    u = r_u(input,kappa_u,kappa_v,eta)
    # drawing eta
    eta_prop = r_eta_prop(input,eta,u,kappa_v)
    # cacluating acceptance
    accept_prob = acceptance_prob(input,eta_prop,eta,kappa_v,u)
    if(runif(1) < accept_prob){
      eta = eta_prop
    }
    eta_samples[i,] = eta
    u_samples[i,] = u
    accept_vec = c(accept_vec, accept_prob)
    kappa_u_samples[i] = kappa_u
    kappa_v_samples[i] = kappa_v
  }
  return(list(
    eta = eta_samples,
    u = u_samples,
    accept = accept_vec,
    kappa_u = kappa_u_samples,
    kappa_v = kappa_v_samples
  ))
}

## ---- res2
# runs the MCMC and stores the time it took
run_time <- system.time(samples <- myMCMC(input, M))
# adding runtime to the list from MCMC
samples$run_time = run_time
# saving for use in other code
save(samples,file = "data/samples.Rdata")
save(input,file ="data/input.Rdata")
cat("Average acceptance rate:", mean(samples$accept))
cat("Run time of MCMC:", as.numeric(run_time[1]))
## ---- break
runaccept <- data.frame(run_time = as.numeric(samples$run_time[1]), 
                        accept = mean(samples$accept))
save(file = "../code/data/runaccept.Rdata",runaccept)
## ---- timeAccept
load("../code/data/runaccept.Rdata")
cat("Average acceptance rate:", runaccept$accept)
cat("Run time of MCMC:", runaccept$run_time)
