rm(list = ls())

# Use log densities

# Libraries ------------------------------------
library(spam)
library(ggplot2)
library(tidyverse)
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
  R = R
  )

get_b <- function(input,z){
  return(input$y+input$E*exp(z)*(z-1))
}

get_c <- function(input,z){
  return(input$E*exp(z))
}

d_eta <- function(input,eta,kappa_v,u){
  return(-1/2*t(eta)%*%diag.spam(kappa_v,input$n)%*%eta + 
           t(eta)%*%(kappa_v*u) + t(eta)%*%input$y-t(exp(eta))%*%input$E)
}

acceptance_prob <- function(input,eta_prop,eta,kappa_v,u){
  accept <- exp(d_eta(input,eta_prop$sample,kappa_v,u) - 
                  d_eta(input,eta$sample,kappa_v,u) -
                  eta_prop$prob + 
                  eta$prob)
  if (accept > 1){
    return(1)
  }
  return(accept)
}


## ---- 21
# Draw samples from the full condition of kappa_u
r_kappa_u <- function(input,u){
  shape = (input$n-1)/2 + input$alpha
  rate = 1/2* t(u) %*% input$R %*% u + input$beta
  return(rgamma(n = 1, shape = shape, rate = rate))
}

## ---- 22
# Draw samples from the full condition of kappa_v
r_kappa_v <- function(input,eta,u){
  shape = input$n/2 + input$alpha
  rate = 1/2 * (t(eta - u)) %*% (eta - u) + input$beta
  return(rgamma(n = 1 ,shape = shape, rate = rate))
}

## ---- 23
# draw samples from the full conditional of u
r_u <- function(input,kappa_u,kappa_v,eta){
  Q = diag.spam(kappa_v, input$n) + kappa_u*input$R
  b = kappa_v * eta
  u = c(rmvnorm.canonical(n = 1,  b = b, Q = Q))
  return(u)
}

## ---- 24
# draw samples from the proposal density of eta
r_eta_prop <- function(input,z,u,kappa_v){
  c_vec = get_c(input,z)
  b_vec = get_b(input,z)
  b = kappa_v*u + b_vec
  Q = diag.spam(kappa_v, input$n) + diag.spam(c_vec)
  sample = c(rmvnorm.canonical(n = 1, b = b, Q = Q))
  prob = dmvnorm.canonical(x = sample, b = b, Q = Q, log = TRUE)
  return(list(sample=sample,prob=prob))
}

## ---- 25
# acceptance of proposal eta
M = 50000
myMCMC <- function(M){
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  u = c(rep_len(0.0,input$n))
  eta <- r_eta_prop(input, u, u, kappa_v = 0.0001)
  i = 0
  eta_samples <- matrix(NA,nrow=M,ncol=input$n)
  u_samples <- matrix(NA,nrow=M,ncol=input$n)
  kappa_u_samples <- vector()
  kappa_v_samples <- vector()
  steps <- seq(1,M)
  count = 0
  for (i in steps){
    setTxtProgressBar(pb, i)
    kappa_u = r_kappa_u(input,u)
    kappa_v = r_kappa_v(input,eta$sample,u)
    u = r_u(input,kappa_u,kappa_v,eta$sample)
    eta_prop = r_eta_prop(input,eta$sample,u,kappa_v)
    accept_prob <- acceptance_prob(input,eta_prop,eta,kappa_v,u)
    if(runif(1) < accept_prob){
      eta = eta_prop
      count = count+1
    }
    eta_samples[i,] = eta$sample
    u_samples[i,] = u
    kappa_u_samples = c(kappa_u_samples,kappa_u)
    kappa_v_samples = c(kappa_v_samples,kappa_v)
  }
  vs <- eta_samples - u_samples
  samples <- tibble(
    steps = steps,
    eta = eta_samples[,1],
    v = vs[,1],
    u1 = u_samples[,1],
    u2 = u_samples[,2],
    kappa_u = kappa_u_samples,
    kappa_v = kappa_v_samples
  )
  return(samples)
}
samples <- myMCMC(M)
ggplot(samples,aes(x = steps)) +
  geom_line(aes(y = v))
ggplot(samples,aes(x = steps)) + 
  geom_line(aes(y = kappa_v))
ggplot(samples,aes(x = steps)) + 
  geom_line(aes(y = kappa_u))
ggplot(samples,aes(x = steps)) + 
  geom_line(aes(y = u1))
ggplot(samples,aes(x = steps)) + 
  geom_line(aes(y = u2))
ggplot(samples,aes(x = steps)) +
  geom_line(aes(y = eta))

cat("Accepted iterations: ",samples$count)

acf(samples$kappa_u)
acf(samples$kappa_v)
acf(samples$v)
acf(samples$u2)
acf(samples$u1)
acf(samples$eta)
