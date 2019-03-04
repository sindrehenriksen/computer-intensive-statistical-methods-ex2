#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# this is the solution to excercise 3
# it contains mostly plotting of data from excercise 2
library(tidyverse)
library(ggpubr)
load("input.Rdata")
load("samples.Rdata")
# Using a burning 
r_cols <- sample(1:input$n,3,replace = F)
M <- length(samples$eta[,1])
steps <- seq(burnin,M)
v <- data.frame(
  steps = steps,
  v1 = samples$eta[steps,r_cols[1]] - samples$u[steps,r_cols[1]],
  v2 = samples$eta[steps,r_cols[2]] - samples$u[steps,r_cols[2]],
  v3 = samples$eta[steps,r_cols[3]] - samples$u[steps,r_cols[3]]
)

u <- data.frame(
  steps = steps,
  u1 = samples$u[steps,r_cols[1]],
  u2 = samples$u[steps,r_cols[2]],
  u3 = samples$u[steps,r_cols[3]]
)

kappa <- data.frame(
  steps = steps,
  kappa_u = samples$kappa_u[steps],
  kappa_v = samples$kappa_v[steps]
)


## ---- 3a
# trace plots
# v
ggarrange(
  ggplot(v, aes(x = steps, y = v1)) + 
    geom_line(colour = "slateblue"),
  ggplot(v, aes(x = steps,y = v2)) + 
    geom_line(colour = "deeppink"),
  ggplot(v, aes(x = steps, y = v3)) + 
    geom_line(colour = "grey41"),
  nrow = 3,
  ncol = 1,
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3]))
)
ggsave("../figures/trace_v.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

# u
ggarrange(
  ggplot(u, aes(x = steps, y = u1)) + 
    geom_line(colour = "slateblue"),
  ggplot(u, aes(x = steps,y = u2)) + 
    geom_line(colour = "deeppink"),
  ggplot(u, aes(x = steps, y = u3)) + 
    geom_line(colour = "grey41"),
  nrow = 3,
  ncol = 1,
  labels = c(sprintf("m = %d",r_cols[1]), sprintf("m = %d",r_cols[2]), sprintf("m = %d",r_cols[3]))
)
ggsave("../figures/trace_u.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

# kappa
ggarrange(
  ggplot(kappa, aes(x = steps,y = kappa_u)) + 
    geom_line(colour = "slateblue"),
  ggplot(kappa, aes(x = steps, y = kappa_v)) + 
    geom_line(colour = "deeppink"),
  nrow = 2,
  ncol = 1
)
ggsave("../figures/trace_kappa.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

## ---- 3b
# autocorrelation plots
# v
