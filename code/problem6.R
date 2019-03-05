#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# Libraries
library(INLA)
library(spam)
library(fields, warn.conflict=FALSE)
library(colorspace)

# Load and prepare data
#str(Oral)
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
summary(result_b_lin)

# Define formula and run INLA with smoking as non-linear effect (rw2)
formula_b_rw2 = Y ~
  -1 +
  f(region_struct, model="besag", graph.file=g, constr=FALSE, hyper=hyper) +
  f(region_random, model="iid", hyper=hyper) +
  finla. model="rw2")
result_b_rw2 = inla(formula_b_rw2, family="poisson", data=data, E=E,
                    control.compute=list(dic=TRUE))
summary(result_b_rw2)

# Get DIC for all the models
smoking_none_DIC = result_a$dic$dic
smoking_lin_DIC = result_b_lin$dic$dic
smoking_rw2_DIC = result_b_rw2$dic$dic
print(smoking_none_DIC)
print(smoking_lin_DIC)
print(smoking_rw2_DIC)
