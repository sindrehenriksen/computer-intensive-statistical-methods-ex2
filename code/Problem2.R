#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- libs2
library(spam)
library(tidyverse)
library(fields, warn.conflicts = FALSE)
source("./data/ex2_additionalFiles/dmvnorm.R")
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
# b values 
get_b <- function(input,z){
  return(input$y+input$E*exp(z)*(z-1))
}

# c values
get_c <- function(input,z){
  return(input$E*exp(z))
}


# Draw samples from the full condition of kappa_u
r_kappa_u <- function(input,u){
  shape = (input$n-1)/2 + input$alpha
  rate = 1/2* t(u) %*% input$R %*% u + input$beta
  return(rgamma(n = 1, shape = shape, rate = rate))
}


# Draw samples from the full condition of kappa_v
r_kappa_v <- function(input,eta,u){
  shape = input$n/2 + input$alpha
  rate = 1/2 * t(eta - u) %*% (eta - u) + input$beta
  return(rgamma(n = 1 ,shape = shape, rate = rate))
}


# draw samples from the full conditional of u
r_u <- function(input,kappa_u,kappa_v,eta){
  Q = diag.spam(kappa_v, input$n) + kappa_u*input$R
  b = kappa_v * eta
  return(c(rmvnorm.canonical(n = 1,  b = b, Q = Q)))
}


# draw samples from the proposal density of eta
r_eta_prop <- function(input,z,u,kappa_v){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  b = kappa_v*u + b_vec
  Q = diag.spam(kappa_v, input$n) + diag.spam(c_vec)
  return(c(rmvnorm.canonical(n = 1, b = b, Q = Q)))
}


# finding the probability of eta from the
# full condition p(eta|...)
d_eta_p <- function(input,eta,kappa_v,u){
  return(-1/2*t(eta)%*%diag.spam(kappa_v,input$n)%*%eta + 
           t(eta)%*%(kappa_v*u) + 
           t(eta)%*%input$y - 
           t(exp(eta))%*%input$E)
}

# finding the probability of eta from the
# estimated full condition q(eta|...)
d_eta_q <- function(input,eta,z,kappa_v,u){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  b = kappa_v*u + b_vec
  Q = diag.spam(kappa_v,input$n) + diag.spam(c_vec)
  return(dmvnorm.canonical(x = eta, b = b, Q = Q, log = TRUE))
}

# calculating the acceptance probability
acceptance_prob <- function(input,eta_prop,eta,kappa_v,u){
  return(min(1,exp(
      d_eta_p(input, eta_prop, kappa_v, u) + 
      d_eta_q(input, eta, eta_prop, kappa_v, u) - 
      d_eta_p(input, eta, kappa_v, u) - 
      d_eta_q(input, eta_prop, eta, kappa_v, u))))
}


# running a MCMC
M <- 70000
myMCMC <- function(input, M){
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  # choosing kappa from the prior
  kappa_u = rgamma(n = 1, shape = input$alpha, rate = input$beta)
  kappa_v = rgamma(n = 1, shape = input$alpha, rate = input$beta)
  # choosing u around the mean
  u = c(rep_len(0.0, input$n))
  eta <- r_eta_prop(input,u,u,kappa_v)
  eta_samples <- matrix(NA,nrow=M,ncol=input$n)
  u_samples <- matrix(NA,nrow=M,ncol=input$n)
  kappa_u_samples <- vector()
  kappa_v_samples <- vector()
  accept_vec <- vector()
  for (i in seq(1,M)){
    setTxtProgressBar(pb, i)
    kappa_u = r_kappa_u(input,u)
    kappa_v = r_kappa_v(input,eta,u)
    u = r_u(input,kappa_u,kappa_v,eta)
    eta_prop = r_eta_prop(input,eta,u,kappa_v)
    accept_vec <- c(accept_vec, acceptance_prob(input,eta_prop,eta,kappa_v,u))
    if(runif(1) < accept_vec[i]){
      eta = eta_prop
    }
    eta_samples[i,] = eta
    u_samples[i,] = u
    kappa_u_samples = c(kappa_u_samples,kappa_u)
    kappa_v_samples = c(kappa_v_samples,kappa_v)
  }
  return(list(
    eta = eta_samples,
    u = u_samples,
    accept_prob <- accept_vec,
    kappa_u = kappa_u_samples,
    kappa_v = kappa_v_samples
  ))
}

## ---- save2
run_time <- system.time(samples <- myMCMC(input, M))
samples$run_time = run_time
save(samples,file = "data/samples.Rdata")
save(input,file ="data/input.Rdata")
