# run simulations to explore scaling model
# author: sebastian daza


# libraries
library(data.table)
library(brms)
library(texreg)
library(ggplot2)
library(patchwork)
library(latex2exp)
library(foreach)
library(doParallel)
source("src/simulation/utils.R")

set.seed(124907)

# create data
E = rnorm(10000, 0, 1)
nreplicates = 100
datasets = list(
    dts = replicate(nreplicates, simScaling(E), simplify = FALSE),
    dti = replicate(nreplicates, simInteraction(E), simplify = FALSE)
)
head(datasets[["dti"]][[sample(1:100, 1)]])

# summary(datasets[[1]][[1]])
# run simple linear regression models
# cnames = c("Scaling", "Interaction", "Scaling + Interaction", "Alternative")
cnames = c("Scaling", "GxE no scaling")
# models = list()
# models[[1]] = lm(y ~ g + E + g * E, data = datasets[[1]][[1]])
# models[[2]] = lm(y ~ g + E + g * E, data = datasets[[2]][[1]])
# models[[3]] = lm(y ~ g + E + g * E, data = datasets[[3]][[1]])
# models[[4]] = lm(y ~ g + E + g * E, data = datasets[[4]][[1]])
# screenreg(models, custom.model.names = cnames)

# distributional model using bayesian stats
model_names = paste0("m", seq_along(datasets))
# models = list()

# scaling
cl = makeCluster(10)
registerDoParallel(cl)
models = foreach(i = 1:nreplicates, 
    .packages = c("brms", "data.table")) %dopar% {
        f = bf(y ~ g + E + g * E, sigma ~ 1 + E)
        dat = datasets[["dts"]][[i]]
        m = brm(f, 
            data = dat, 
            family = brmsfamily("gaussian", link_sigma = "log"), 
            chains = 1, backend = "cmdstanr", cores = 1)
        return(m)
}
stopCluster(cl)
m1 = combine_models(mlist = models, check_data = FALSE)
rm(models)
summary(m1)

# interaction
cl = makeCluster(10)
registerDoParallel(cl)
models = foreach(i = 1:nreplicates, 
    .packages = c("brms", "data.table")) %dopar% {
        f = bf(y ~ g + E + g * E, sigma ~ 1 + E)
        dat = datasets[["dti"]][[i]]
        m = brm(f, 
            data = dat, 
            family = brmsfamily("gaussian", link_sigma = "log"), 
            chains = 1, backend = "cmdstanr", cores = 1)
        return(m)
}
stopCluster(cl)
m2 = combine_models(mlist = models, check_data = FALSE)
rm(models)

# convergence testing 
models = list(m1, m2)
rm(m1, m2)

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
        include.rsquared =  FALSE,
        include.nobs = TRUE,
        include.loo.ic = FALSE,
        include.waic = FALSE,
        file = "output/tables/model-tables.tex")
)
file.copy("output/tables/model-tables.tex", "manuscript/tables", 
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

# correlation vs slope

# scaling
test = data.table(datasets[["dts"]][[sample(1:100, 1)]])
test[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
# ggplot(test, aes(g, y, color = as.factor(qE))) + geom_point()

plots = list()
for (i in 1:10) {
    slope = specify_decimal(coef(lm(y ~ g, data = test[qE == i]))[2], 2)
    corr = specify_decimal(cor(test[qE == i, .(y, g)])[1, 2], 2)
    plots[[i]] = ggplot(test[qE == i], aes(g, y)) + geom_point(size = 0.5, 
        alpha = 0.1) +
        geom_smooth(method = "lm", color = "red", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("E", i), 
            subtitle = paste0("Correlation: ", corr, "; Slope: ", slope)) +
        theme_minimal()
}

savepdf("output/plots/slope-cor-gxe-scaling", 25, 20)
print(wrap_plots(plots, ncol = 3))
dev.off()
file.copy("output/plots/slope-cor-gxe-scaling.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)  

# no scaling
test = data.table(datasets[["dti"]][[sample(1:100, 1)]])
test[, qE := cut(E, quantile(E, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]

plots = list()
for (i in 1:10) {
    slope = specify_decimal(coef(lm(y ~ g, data = test[qE == i]))[2], 2)
    corr = specify_decimal(cor(test[qE == i, .(y, g)])[1, 2], 2)
    plots[[i]] = ggplot(test[qE == i], aes(g, y)) + geom_point(size = 0.5, 
        alpha = 0.1) +
        geom_smooth(method = "lm", color = "red", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("E", i), 
            subtitle = paste0("Correlation: ", corr, "; Slope: ", slope)) +
        theme_minimal()
}

savepdf("output/plots/slope-cor-gxe-no-scaling", 25, 20)
print(wrap_plots(plots, ncol = 3))
dev.off()
file.copy("output/plots/slope-cor-gxe-no-scaling.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)  
