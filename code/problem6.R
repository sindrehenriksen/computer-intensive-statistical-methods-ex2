setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- 6b_prelim
# Libraries
library(INLA)
library(spam)
library(fields, warn.conflict=FALSE)
library(colorspace)
library(ggplot2)

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
  f(region_struct, model="besag", graph.file=g, constr=FALSE, hyper=hyper) +
  f(region_random, model="iid", hyper=hyper)
data = list(Y=Y, E=E, region_struct=region, region_random=region)
result_a = inla(formula, family="poisson", data=data, E=E,
                control.compute=list(dic=TRUE))
summary(result_a)

# Plot median of estimate divided by E (Y/E = E*exp(eta)/E = exp(eta)) in map
y_median_a = exp(result_a$summary.random$region_struct$`0.5quant`)
pdf("../figures/y_median_a.pdf")
germany.plot(y_median_a, col=col, legend=TRUE)
dev.off()

# Improve estimates of the posterior marginals
result_a = inla.hyperpar(result_a)

# Plot
kappa_u_marginal = result_a$marginals.hyperpar$`Precision for region_struct`
kappa_v_marginal = result_a$marginals.hyperpar$`Precision for region_random`
plot(inla.smarginal(kappa_u_marginal))
plot(inla.smarginal(kappa_v_marginal))

## ---- 6b
smoking = read.table("./data/ex2_additionalFiles/smoking.dat")
data["smoking"] = smoking

# Define formula and run INLA with smoking as a linear effect
formula_b_lin = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE, hyper=hyper) +
  f(region_random, model="iid", hyper=hyper) +
  smoking
result_b_lin = inla(formula_b_lin, family="poisson", data=data, E=E,
                    control.compute=list(dic=TRUE))
# summary(result_b_lin)

# Define formula and run INLA with smoking as non-linear effect (rw2)
formula_b_rw2 = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE, hyper=hyper) +
  f(region_random, model="iid", hyper=hyper) +
  f(smoking, model="rw2")
result_b_rw2 = inla(formula_b_rw2, family="poisson", data=data, E=E,
                    control.compute=list(dic=TRUE))
summary(result_b_rw2)

# Get DIC for all the models
smoking_none_DIC = result_a$dic$dic
smoking_lin_DIC = result_b_lin$dic$dic
smoking_rw2_DIC = result_b_rw2$dic$dic

# Plot posterior median and 95% CI of non-linear effect, compare with linear
smoking_beta = result_b_lin$summary.fixed$mean
smoking_val = result_b_rw2$summary.random$smoking$ID
smoking_eff_025 = result_b_rw2$summary.random$smoking$`0.025quant`
smoking_eff_5 = result_b_rw2$summary.random$smoking$`0.5quant`
smoking_eff_975 = result_b_rw2$summary.random$smoking$`0.975quant`

data = data.frame(
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
# plot_smoking

## ---- save
ggsave("../figures/smoking_effect.pdf",
       plot = plot_smoking, scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)
