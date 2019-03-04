#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())


library(spam)
library(tidyverse)
library(fields, warn.conflicts = FALSE)
source("./data/ex2_additionalFiles/dmvnorm.R")
library(colorspace)
library(ggpubr)
library(gridExtra)
## ---- 2

load("./data/ex2_additionalFiles/tma4300_ex2_Rmatrix.Rdata")
attach(Oral)
col <- diverge_hcl(8) # blue - red
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)

# a list of all the input variables to make code more readable
input <- list(
  y = Oral$Y, 
  E = Oral$E, 
  n = length(Oral$Y), 
  alpha = 1, 
  beta = 0.01,
  R = R
  )

# getting the b-values
get_b <- function(input,z){
  return(input$y+input$E*exp(z)*(z-1))
}

# getting the c values
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
  sample = c(rmvnorm.canonical(n = 1, b = b, Q = Q))
  prob = dmvnorm.canonical(x = sample, b = b, Q = Q,log = TRUE)
  return(list(sample=sample,prob=prob))
}


# finding the probability of eta
d_eta <- function(input,eta,kappa_v,u){
  return(-1/2*t(eta)%*%diag.spam(kappa_v,input$n)%*%eta + 
           t(eta)%*%(kappa_v*u) + 
           t(eta)%*%input$y - 
           t(exp(eta))%*%input$E)
}


# calculating the acceptance probability
acceptance_prob <- function(input,eta_prop,eta,kappa_v,u){
  return(min(1,exp(
    d_eta(input,eta_prop$sample,kappa_v,u) + 
      r_eta_prop(input,eta_prop$sample,u,kappa_v)$prob - 
      eta_prop$prob - 
      d_eta(input,eta$sample,kappa_v,u))))
}


# running a MCMC
M <- 70000
burnin <- 10000
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
for (i in seq(1,M)){
  setTxtProgressBar(pb, i)
  kappa_u = r_kappa_u(input,u)
  kappa_v = r_kappa_v(input,eta$sample,u)
  u = r_u(input,kappa_u,kappa_v,eta$sample)
  eta_prop = r_eta_prop(input,eta$sample,u,kappa_v)
  accept_prob <- acceptance_prob(input,eta_prop,eta,kappa_v,u)
  if(runif(1) < accept_prob){
    eta = eta_prop
  }
  eta_samples[i,] = eta$sample
  u_samples[i,] = u
  kappa_u_samples = c(kappa_u_samples,kappa_u)
  kappa_v_samples = c(kappa_v_samples,kappa_v)
}

v <- data.frame(
  steps = seq(1:M),
  v1 = eta_samples[,1] - u_samples[,1],
  v2 = eta_samples[,242] - u_samples[,242],
  v3 = eta_samples[,493] - u_samples[,493]
)

u <- data.frame(
  steps = seq(1:M),
  u1 = u_samples[,1],
  u2 = u_samples[,242],
  u3 = u_samples[,493]
)

kappa <- data.frame(
  steps = seq(1,M),
  kappa_u = kappa_u_samples,
  kappa_v = kappa_v_samples
)

# Plotting v samples
grid.arrange(
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = v1, colour = "slateblue")),
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = v2, colour = "deeppink")),
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = v3, colour = "grey41")),
  nrow = 3,
  ncol = 1
  )

# Plotting u samples
ggarrange(
  ggplot(u, aes(x = steps)) + 
    geom_line(aes(y = u1)),
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = u2)),
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = u3))
)


# Plotting Kappas
ggarrange(
  ggplot(kappa, aes(x = steps)) + 
    geom_line(aes(y = kappa_u)),
  ggplot(v, aes(x = steps)) + 
    geom_line(aes(y = kappa_v))
)

# test MCMC
eta <- data.frame(
  steps = seq(burnin:M),
  sam = eta_samples[burnin:M,1]
)
ggplot(eta, aes(x = steps)) + 
  geom_line(aes(y = sam))

acf(eta$sam)
# --- 
  
ggsave("../figures/posteriorSamps.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)