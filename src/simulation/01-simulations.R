# run simulations to explore scaling model
# author: sebastian daza


# libraries
library(data.table)
library(brms)
library(texreg)
library(ggplot2)
library(latex2exp)
library(future)
plan(multiprocess, workers = 10)
options(future.rng.onMisuse = "ignore")
source("src/simulation/utils.R")

# create data
E = rnorm(10000, 0, 1)
nreplicates = 100

datasets = list(
    dts = replicate(nreplicates, simScaling(E), simplify = FALSE),
    dti = replicate(nreplicates, simInteraction(E), simplify = FALSE),
    dtsi = replicate(nreplicates, simScalingInteraction(E), simplify = FALSE),
    dalt = replicate(nreplicates, simAlternative(E), simplify = FALSE)
)

summary(datasets[[1]][[1]])

# run simple linear regression models
cnames = c("Scaling", "Interaction", "Scaling + Interaction", "Alternative")
models = list()
models[[1]] = lm(y ~ g + E + g * E, data = datasets[[1]][[1]])
models[[2]] = lm(y ~ g + E + g * E, data = datasets[[2]][[1]])
models[[3]] = lm(y ~ g + E + g * E, data = datasets[[3]][[1]])
models[[4]] = lm(y ~ g + E + g * E, data = datasets[[4]][[1]])
screenreg(models, custom.model.names = cnames)

# distributional model using bayesian stats
model_names = paste0("b", seq_along(datasets))
models = list()
f = bf(y ~ g + E + g * E, sigma ~ 1 + E)
for (i in seq_along(datasets)) {
    print(paste0("::::: Data ", i, " ::::::"))
    models[[model_names[i]]] = brm_multiple(f, data = datasets[[i]], 
        family = brmsfamily("gaussian", link_sigma = "log"), chains = 1)
    class(models[[model_names[i]]]) = "brmsfit"
}

# convergence testing 
rhats = NULL
for (i in seq_along(models)) { 
    rhats  = c(rhats, as.vector(unlist(models[[i]]$rhats)))
}
table(rhats > 1.05 | rhats < 0.95)

# print table
cmap = list(
    "Intercept"  = "$\\tau_0$", 
    "E" = "$\\tau_1$", 
    "g" = "$\\pi_0$", 
    "g:E" = "$\\pi_1$", 
    "sigma_Intercept" = "$\\lambda_0$",
    "sigma_E" = "$\\lambda_1$")
 
suppressWarnings(
    texreg(models, 
        caption = "Distributional models",, 
        custom.coef.map = cmap,
        label = "tab:dist_models",
        custom.model.names = cnames,
        caption.above = TRUE,
        booktabs = TRUE,
        dcolumn = TRUE,
        use.packages = FALSE, 
        use.HDI = FALSE,
        include.rsquared = TRUE,
        include.nobs = TRUE,
        include.loo.ic = FALSE,
        include.waic = FALSE,
        file = "output/tables/model_tables.tex")
)
file.copy("output/tables/model_tables.tex", "manuscript/tables", 
    recursive = TRUE)

# scaling test and pltos
for (i in seq_along(models)) {
    t = scalingTest(models[[i]])
    plotScalingTest(t, paste0("output/plots/", model_names[i], "_test"))
    file.copy(paste0("output/plots/", model_names[i], "_test.pdf"), "manuscript/plots", 
        recursive = TRUE)
    plotDecomp(models[[i]], paste0("output/plots/", model_names[i], "_decomp"))
    file.copy(paste0("output/plots/", model_names[i], "_decomp.pdf"), "manuscript/plots", 
        recursive = TRUE)
}
print(":::::::: plots created ::::::")
