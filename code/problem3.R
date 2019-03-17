setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# this is the solution to excercise 3
# it contains mostly plotting of data from excercise 2
## ---- libs3
library(tidyverse)
library(spam)
library(ggpubr)
set.seed(123)
load("data/input.Rdata")
load("data/samples.Rdata")

## ---- dfCreate
r_cols <- sort(sample(1:input$n,3,replace = F))
M <- length(samples$eta[,1])
burnin = 5000
steps <- seq(1,M)
burnin_step = seq(burnin,M)

MCMC_list <- data.frame(
  steps = steps,
  v1 = samples$eta[steps,r_cols[1]] - samples$u[steps,r_cols[1]],
  v2 = samples$eta[steps,r_cols[2]] - samples$u[steps,r_cols[2]],
  v3 = samples$eta[steps,r_cols[3]] - samples$u[steps,r_cols[3]],
  u1 = samples$u[steps,r_cols[1]],
  u2 = samples$u[steps,r_cols[2]],
  u3 = samples$u[steps,r_cols[3]],
  kappa_u = samples$kappa_u[steps],
  kappa_v = samples$kappa_v[steps],
  is_burnin = c(rep(TRUE,burnin),rep(FALSE,M-burnin)))
save(MCMC_list,file = "data/mcmc_df.Rdata")

## ---- 3a
# trace plots with burn in removed
fig_3a <- ggarrange(
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = v1, color = is_burnin)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = v2, color = is_burnin)) +
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = v3)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = u1)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = u2)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = u3)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps,y = kappa_u)) + 
    geom_line(color = "grey24") + rremove("x.text") + rremove("xlab") + rremove("ylab"),
  ggplot(MCMC_list[burnin_step,], aes(x = steps, y = kappa_v)) + 
    geom_line(color = "grey24") + rremove("ylab"),
  nrow = 8,
  ncol = 1,
  labels = c(sprintf("v[m = %d]",r_cols[1]), sprintf("v[m = %d]",r_cols[2]), sprintf("v[m = %d]",r_cols[3]),
             sprintf("u[m = %d]",r_cols[1]), sprintf("u[m = %d]",r_cols[2]), sprintf("u[m = %d]",r_cols[3])
             ,"kappa_u","kappa_v"),
  font.label = list(size = 10, color = "firebrick1"),
  label.x = 0,
  label.y = 1,
  hjust= -1.7,
  vjust = 1.9
)
fig_3a

## ---- break
ggsave("../figures/trace.pdf", plot = fig_3a, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- 3b
# autocorrelation plots with burn in removed
MCMC_acf <- data.frame(
  lag = acf(MCMC_list$kappa_v[burnin_step],plot = F,lag.max = 200)$lag, 
  v1 = acf(MCMC_list$v1[burnin_step],plot = F,lag.max = 200)$acf,
  v2 = acf(MCMC_list$v2[burnin_step],plot = F,lag.max = 200)$acf, 
  v3 = acf(MCMC_list$v3[burnin_step],plot = F,lag.max = 200)$acf,
  u1 = acf(MCMC_list$u1[burnin_step],plot = F,lag.max = 200)$acf,
  u2 = acf(MCMC_list$u2[burnin_step],plot = F,lag.max = 200)$acf, 
  u3 = acf(MCMC_list$u3[burnin_step],plot = F,lag.max = 200)$acf,
  kappa_u = acf(MCMC_list$kappa_u[burnin_step],plot = F,lag.max = 200)$acf,
  kappa_v = acf(MCMC_list$kappa_v[burnin_step],plot = F,lag.max = 200)$acf
)


fig_3b <- ggarrange(
  ggplot(MCMC_acf,aes(x = lag, y = v1)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") +
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = v2)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = v3)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = u1)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = u2)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = u3)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = kappa_u)) + rremove("x.text") + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(MCMC_acf,aes(x = lag, y = kappa_v)) + 
    geom_hline(aes(yintercept = 0)) + rremove("xlab") + rremove("ylab") + 
    geom_segment(aes(xend = lag, yend = 0)),
  nrow = 8,
  ncol = 1,
  labels = c(sprintf("v[m = %d]",r_cols[1]), sprintf("v[m = %d]",r_cols[2]), sprintf("v[m = %d]",r_cols[3]),
             sprintf("u[m = %d]",r_cols[1]), sprintf("u[m = %d]",r_cols[2]), sprintf("u[m = %d]",r_cols[3])
             ,sprintf("  kappa_u"),sprintf("  kappa_v")),
  font.label = list(size = 10, color = "firebrick1"),
  label.x = 0,
  label.y = 1,
  hjust= -1.7,
  vjust = 1.9
)
fig_3b

## ---- break
ggsave("../figures/acf.pdf", plot = fig_3b, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)


## ---- 3c

library(coda)
library(kableExtra)
z_scores <- geweke.diag(MCMC_list[,2:9], frac1=0.1, frac2=0.5)$z
geweke_diag <- data.frame(
  z_scores = z_scores,
  p_values = 2*pnorm(abs(z_scores),lower.tail = FALSE)
)
kable(geweke_diag,caption = "Result of Gweke Statistic on the parameters",booktabs = T) 
## ---- break
save(file = "../code/data/geweke_diag.Rdata",geweke_diag)
## ---- tableGweke
load("../code/data/geweke_diag.Rdata")
kable(geweke_diag,caption = "Result of Gweke Statistic on the parameters",booktabs = T) 

## ---- break
## ---- plotgeweke
testBurningGeweke <- function(MCMC_list, M){
  z_scores_burnin <- data.frame()
  burn = numeric()
  z_name = numeric()
  z = numeric()
  for (i in seq(1,11)){
    burnin[i] = 500*(i-1)
    z_score<-geweke.diag(MCMC_list[seq(burnin[i],M),2:9], frac1=0.1, frac2=0.5)$z
    burn = c(burn,rep(burnin[i],8))
    z = c(z,as.vector(z_score))
    z_name = c(z_name,names(z_score))
  }
  return(data.frame(burnin = burn, z_statistic = z, Parameter = z_name))
}
z_scores_burnin<-testBurningGeweke(MCMC_list, M)
burnin_test_plot<-ggplot()+
  geom_point(data= z_scores_burnin,aes(x=z_statistic,y=burnin,color=Parameter)) +
  geom_rect(aes(xmin=-1.6,xmax = 1.6, ymin=-100,ymax=5100 ), fill = "blue",alpha = 0.2)+ 
  ylab("Burn-in")+
  xlab("Z-statistic")+
  labs(colour="Parameter")+
  xlim(c(-2,2))
burnin_test_plot
## ---- break
ggsave("../figures/test_burnin.pdf", plot = burnin_test_plot, device = NULL, path = NULL,
       scale = 1, width = 8, height = 8, units = "in",
       dpi = 300, limitsize = TRUE)
  