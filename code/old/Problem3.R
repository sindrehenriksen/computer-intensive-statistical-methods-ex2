#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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

## ---- break
v <- data.frame(
  steps = steps,
  v1 = samples$eta[steps,r_cols[1]] - samples$u[steps,r_cols[1]],
  v2 = samples$eta[steps,r_cols[2]] - samples$u[steps,r_cols[2]],
  v3 = samples$eta[steps,r_cols[3]] - samples$u[steps,r_cols[3]],
  is_burnin = c(rep(TRUE,burnin),rep(FALSE,M-burnin))
)

u <- data.frame(
  steps = steps,
  u1 = samples$u[steps,r_cols[1]],
  u2 = samples$u[steps,r_cols[2]],
  u3 = samples$u[steps,r_cols[3]],
  is_burnin = c(rep(TRUE,burnin),rep(FALSE,M-burnin))
)

kappa <- data.frame(
  steps = steps,
  kappa_u = samples$kappa_u[steps],
  kappa_v = samples$kappa_v[steps],
  is_burnin = c(rep(TRUE,burnin),rep(FALSE,M-burnin))
)

## ---- 3a
# trace plots
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
ggsave("../figures/trace.pdf", plot = fig_3a, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- tracev
fig_v <- ggarrange(
  ggplot(v, aes(x = steps, y = v1, color = is_burnin)) + 
    geom_line()+
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  ggplot(v, aes(x = steps, y = v2, color = is_burnin)) +
    geom_line()+
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  ggplot(v, aes(x = steps, y = v3, color = is_burnin)) + 
    geom_line()+
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  nrow = 3,
  ncol = 1,
  common.legend = TRUE,
  legend = "bottom",
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3])),
  label.x = 0,
  label.y = 1,
  hjust= -1.7,
  vjust = 1.9
)
fig_v <- annotate_figure(fig_v,
                top = text_grob("Traceplot of v", color = "black", face = "bold", size = 14)
)
## ---- break
ggsave("../figures/trace_v.pdf", plot = fig_v, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- traceu
fig_u <- ggarrange(
  ggplot(u, aes(x = steps, y = u1, color = is_burnin)) + 
    geom_line() +
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  ggplot(u, aes(x = steps, y = u2, color = is_burnin)) + 
    geom_line() + 
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  ggplot(u, aes(x = steps, y = u3, color = is_burnin)) + 
    geom_line() + 
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  nrow = 3,
  ncol = 1,
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3]))
)
fig_u <- annotate_figure(fig_u,
                top = text_grob("Traceplot of u", color = "black", face = "bold", size = 14)
)
fig_u
## ---- break
ggsave("../figures/trace_u.pdf", plot = fig_u, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- tracekappa
fig_kappa <- ggarrange(
  ggplot(kappa, aes(x = steps,y = kappa_u, color = is_burnin)) + 
    geom_line() +
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  ggplot(kappa, aes(x = steps, y = kappa_v, color = is_burnin)) + 
    geom_line() + 
    geom_vline(xintercept = burnin, color="firebrick1")+
    scale_color_manual(
      name = "",
      labels = c("samples","burnin"),
      values = c("grey24","grey47")),
  nrow = 2,
  ncol = 1,
  labels = c("kappa_u", "kappa_v")
)
fig_kappa <- annotate_figure(fig_kappa,
                         top = text_grob("Traceplot of kappa_u and kappa_v", color = "black", face = "bold", size = 14),
                         fig.lab = "trace_kappa", fig.lab.face = "bold"
)
fig_kappa
## ---- break
ggsave("../figures/trace_kappa.pdf", plot = fig_kappa, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- 3b
# autocorrelation plots

MCMC_acf <- data.frame(
  lag = acf(v$v1[burnin_step],plot = F)$lag, 
  v1 = acf(v$v1[burnin_step],plot = F)$acf,
  v2 = acf(v$v2[burnin_step],plot = F)$acf, 
  v3 = acf(v$v3[burnin_step],plot = F)$acf,
  u1 = acf(u$u1[burnin_step],plot = F)$acf,
  u2 = acf(u$u2[burnin_step],plot = F)$acf, 
  u3 = acf(u$u3[burnin_step],plot = F)$acf,
  kappa_u = acf(kappa$kappa_u[burnin_step],plot = F)$acf,
  kappa_v = acf(kappa$kappa_v[burnin_step],plot = F)$acf
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

## ---- acfv
v_acf <- data.frame(
  lag = acf(v$v1,plot = F)$lag, 
  v1 = acf(v$v1,plot = F)$acf,
  v2 = acf(v$v2,plot = F)$acf, 
  v3 = acf(v$v3,plot = F)$acf
)


fig_v_acf <- ggarrange(
  ggplot(v_acf,aes(x = lag, y = v1)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(v_acf,aes(x = lag, y = v2)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(v_acf,aes(x = lag, y = v3)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  nrow = 3,
  ncol = 1,
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3])),
  label.x = 0,
  label.y = 1,
  hjust= -0.5,
  vjust = 3
)
fig_v_acf <- annotate_figure(fig_v_acf,
                             fig.lab = "acf_v", fig.lab.face = "bold"
)
fig_v_acf
## ---- break
ggsave("../figures/acf_v.pdf", plot = fig_v_acf, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- acfu
u_acf <- data.frame(
  lag = acf(u$u1,plot = F)$lag, 
  u1 = acf(u$u1,plot = F)$acf,
  u2 = acf(u$u2,plot = F)$acf, 
  u3 = acf(u$u3,plot = F)$acf
)

fig_u_acf <- ggarrange(
  ggplot(u_acf,aes(x = lag, y = u1)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(u_acf,aes(x = lag, y = u2)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(u_acf,aes(x = lag, y = u3)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  nrow = 3,
  ncol = 1,
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3])),
  label.x = 0,
  label.y = 1,
  hjust= -0.5,
  vjust = 3
  
)
fig_u_acf <- annotate_figure(fig_u_acf,
                             fig.lab = "acf_u", fig.lab.face = "bold"
)
fig_u_acf
## ---- break
ggsave("../figures/acf_u.pdf", plot = fig_u_acf, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)


## ---- acfkappa
kappa_acf <- data.frame(
  lag = acf(kappa$kappa_u,plot = F)$lag, 
  kappa_u = acf(kappa$kappa_u,plot = F)$acf,
  kappa_v = acf(kappa$kappa_v,plot = F)$acf
)


fig_kappa_acf <- ggarrange(
  ggplot(kappa_acf,aes(x = lag, y = kappa_u)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  ggplot(kappa_acf,aes(x = lag, y = kappa_v)) + 
    geom_hline(aes(yintercept = 0)) + xlab("") + 
    geom_segment(aes(xend = lag, yend = 0)),
  nrow = 2,
  ncol = 1
)
fig_kappa_acf
## ---- break
ggsave("../figures/acf_kappa.pdf", plot = fig_u_acf, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- 3c
library(coda)
z_scores <- geweke.diag(MCMC_list[,-1], frac1=0.1, frac2=0.5)$z
geweke_diag <- data.frame(
  params = names(z_scores),
  z_scores = z_scores,
  p_values = 2*pnorm(abs(z_scores),lower.tail = FALSE)
)
summary(geweke_diag)

u_z_scores <- geweke.diag(samples$u,frac1 = 0.1,frac2 = 0.5)$z
u_p_value_df <- enframe(2*pnorm(abs(u_z_scores),lower.tail = FALSE))
ggplot(u_p_value_df) + 
  geom_histogram(aes(x = value, y = ..density..),
                 bins = 50,
                 colour = "white",
                 fill = "cornflowerblue"
  ) 

v_z_scores <- geweke.diag(samples$eta - samples$u,frac1 = 0.1,frac2 = 0.5)$z
v_p_value_df <- enframe(2*pnorm(abs(v_z_scores),lower.tail = FALSE))
ggplot(v_p_value_df) + 
  geom_histogram(aes(x = value, y = ..density..),
                 bins = 50,
                 colour = "white",
                 fill = "cornflowerblue"
  ) 
