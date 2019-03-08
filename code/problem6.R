# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- 6_prelim
# Libraries
library(INLA)
library(spam)
library(fields, warn.conflict=FALSE)
library(colorspace)
library(ggplot2)
library(tibble)
library(gridExtra)
library(latex2exp)

# Load and prepare data
# str(Oral)
attach(Oral)
col = diverge_hcl(8)  # blue - red
germany.plot(Oral$Y/Oral$E, col=col, legend=TRUE)
g = system.file("demodata/germany.graph", package="INLA")

# Variables
alpha = 1
beta = 0.01
region = 1:length(Y)

## ---- 6a
# Define formula and run INLA
hyper = list(prec=list(prior="loggamma", param=c(alpha, beta)))
formula = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE,
    hyper=hyper) +
  f(region_random, model="iid", hyper=hyper)
data = list(Y=Y, E=E, region_struct=region, region_random=region)
result_a = inla(formula, family="poisson", data=data, E=E,
                control.compute=list(dic=TRUE))
summary(result_a)

# Plot median of posterior divided by E in map
# (Y/E = E*exp(eta)/E = exp(eta))
y_median_a = exp(result_a$summary.random$region_struct$`0.5quant`)
pdf("../figures/y_median_a.pdf")
germany.plot(y_median_a, col=col, legend=TRUE)
dev.off()

# Improve estimates of the posterior marginals
result_a = inla.hyperpar(result_a)

# Get and smooth marginals
rns = c(157, 222, 429) #######################################
kappa_u_marginal = result_a$marginals.hyperpar$
  `Precision for region_struct`
kappa_u_marginal = inla.smarginal(kappa_u_marginal)
kappa_v_marginal = result_a$marginals.hyperpar$
  `Precision for region_random`
kappa_v_marginal = inla.smarginal(kappa_v_marginal)
i1 = paste0("index.", rns[1])
i2 = paste0("index.", rns[2])
i3 = paste0("index.", rns[3])
u1_marginal = inla.smarginal(
  result_a$marginals.random$region_struct[[i1]])
u2_marginal = inla.smarginal(
  result_a$marginals.random$region_struct[[i2]])
u3_marginal = inla.smarginal(
  result_a$marginals.random$region_struct[[i3]])
v1_marginal = inla.smarginal(
  result_a$marginals.random$region_random[[i1]])
v2_marginal = inla.smarginal(
  result_a$marginals.random$region_random[[i2]])
v3_marginal = inla.smarginal(
  result_a$marginals.random$region_random[[i1]])

# Plot marginals
p_kappa_u = ggplot(as_tibble(kappa_u_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=expression(kappa[u]), y="") +
  theme(plot.margin = unit(c(0.1, 0.1, 0, -0.4), "cm"))
p_kappa_v = ggplot(as_tibble(kappa_v_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=expression(kappa[v]), y="") +
  xlim(c(min(kappa_u_marginal$x), 1000)) +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_u1 = ggplot(as_tibble(u1_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("u_{", rns[1], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm")) +
  xlim(c(min(u1_marginal$x), 1.5))
p_u2 = ggplot(as_tibble(u2_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("u_{", rns[2], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_u3 = ggplot(as_tibble(u3_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("u_{", rns[3], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_v1 = ggplot(as_tibble(v1_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("v_{", rns[1], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_v2 = ggplot(as_tibble(v2_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("v_{", rns[2], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_v3 = ggplot(as_tibble(v3_marginal)) +
  geom_line(aes(x, y)) +
  labs(x=TeX(paste0("v_{", rns[3], "}")), y="") +
  theme(plot.margin = unit(c(0, 0.1, 0, -0.4), "cm"))
p_marginals = grid.arrange(
  p_kappa_u, p_kappa_v, p_u1, p_v1, p_u2, p_v2, p_u3, p_v3,
  ncol=2)
ggsave("../figures/posterior_marginals.pdf", plot=p_marginals,
       width=6, height=7, units="in", dpi=300)

## ---- 6b
smoking = read.table("./data/ex2_additionalFiles/smoking.dat")
data["smoking"] = smoking

# Define formula and run INLA with smoking as a linear effect
formula_b_lin = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE,
    hyper=hyper) +
  f(region_random, model="iid", hyper=hyper) +
  smoking
result_b_lin = inla(formula_b_lin, family="poisson", data=data,
                    E=E, control.compute=list(dic=TRUE))
# summary(result_b_lin)

# Define formula and run INLA with smoking as a rw2 effect
formula_b_rw2 = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE,
    hyper=hyper) +
  f(region_random, model="iid", hyper=hyper) +
  f(smoking, model="rw2")
result_b_rw2 = inla(formula_b_rw2, family="poisson", data=data,
                    E=E, control.compute=list(dic=TRUE))
summary(result_b_rw2)

# Get DIC for all the models
smoking_none_DIC = result_a$dic$dic
smoking_lin_DIC = result_b_lin$dic$dic
smoking_rw2_DIC = result_b_rw2$dic$dic

# Plot posterior median and 95% CI of non-linear effect
smoking_beta = result_b_lin$summary.fixed$mean
smoking_val = result_b_rw2$summary.random$smoking$ID
smoking_eff_025 = result_b_rw2$summary.random$smoking$`0.025quant`
smoking_eff_5 = result_b_rw2$summary.random$smoking$`0.5quant`
smoking_eff_975 = result_b_rw2$summary.random$smoking$`0.975quant`

data = tibble(
  x = smoking_val,
  lin = smoking_beta * smoking_val + smoking_eff_5[1],
  med = smoking_eff_5,
  ci_low = smoking_eff_025,
  ci_high = smoking_eff_975
)
plot_smoking = ggplot(data) +
  geom_ribbon(aes(x, ymin=ci_low, ymax=ci_high), alpha=0.2) +
  geom_line(aes(x, y=lin, col="linear")) +
  geom_line(aes(x, y=med, col="rw2")) +
  geom_line(aes(x, y=ci_low, col="95CI")) +
  geom_line(aes(x, y=ci_high, col="95CI")) +
  # scale_color_manual(values=c("black", "pink", "yellow")) +
  labs(x="smoking", y="effect", col="")
ggsave("../figures/smoking_effect.pdf", plot=plot_smoking,
       width=5.5, height=3, units="in", dpi=300)
