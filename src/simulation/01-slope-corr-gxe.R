# slope vs correlation
# author: sebastian daza


library(data.table)
library(brms)
library(texreg)
library(ggplot2)
library(patchwork)
source("src/utils.R")

# create data
E = rnorm(10000, 0, 1)
nreplicates = 100

simHC = function(E) {
    N = length(E)
    G = rnorm(N,0,1)
    sigma = exp(0.2 + 0.5 * E) 
    y = rnorm(N, 0.2 + G * 0.5 + E * 0.4 + E * G * 0.2, sigma)
    df = data.frame(E = E, y = y, g = G, ys = scale(y))
}

test = data.table(simHC(E))
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

savepdf("output/plots/slope-cor-gxe", 25, 20)
print(wrap_plots(plots, ncol = 3))
dev.off()
file.copy("output/plots/slope-cor-gxe.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    

test[, ys := scale(y)]
test[, gs := scale(g)]
test[, Es := scale(E)]

f = bf(y ~ g + E + g * E, sigma ~ E)
m1 = brm(f, data = test, backend = "cmdstanr", cores = 4)

f = bf(ys ~ g + E + g * E, sigma ~ E)
m2 = brm(f, data = test, backend = "cmdstanr", cores = 4)

f = bf(ys ~ gs + E + gs * E, sigma ~ E)
m3 = brm(f, data = test, backend = "cmdstanr", cores = 4)

f = bf(ys ~ gs + Es + gs * Es, sigma ~ Es)
m4 = brm(f, data = test, backend = "cmdstanr", cores = 4)

models = list(m1, m2, m3, m4)
texreg(models)

screenreg(models,
    use.HDI = FALSE,
    include.random = TRUE,
    include.rsquared = FALSE,
    include.nobs = TRUE,
    include.loo.ic = FALSE,
    reloo = FALSE,
    include.waic = FALSE)

cnames = c("Unstandardized", "Standardized BMI")
custom_coeff_map = list(
    "Intercept" = "Constant", 
    "g" = "G", 
    "E" = "E", 
    "g:E" = "G x E", 
    "gs" = "G-zscore",
    "gs:E" = "G-zscore x E",
    "sigma_Intercept" = "Sigma Constant", 
    "sigma_E" = "Sigma E"
)

screenreg(models[1:2],
    use.HDI = FALSE,
    include.random = TRUE,
    include.rsquared = FALSE,
    include.nobs = TRUE,
    include.loo.ic = FALSE,
    reloo = FALSE,
    include.waic = FALSE,
    custom.coef.map = custom_coeff_map
    )

caption = paste0("GxE models")
groups = list("Coefficients" = 1:4, "Residuals" = 5:6)

texreg::texreg(models[1:2], 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:slope-cor-gxe",
    groups = groups,
    # custom.note = "\\item  $^*$ Null hypothesis value outside the confidence interval. 
    #     \\item All covariates are standardized. Life expectancy estimates were transformed using $ln\\left(-ln( 1-\\frac{e_0}{ (78.6 + 1.05)}\\right)$.",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    use.HDI = FALSE,
    include.random = TRUE,
    include.rsquared = FALSE,
    include.nobs = TRUE,
    include.loo.ic = FALSE,
    reloo = FALSE,
    include.waic = FALSE, 
    file = "output/tables/slope-cor-gxe.tex"
)    
file.copy("output/tables/slope-cor-gxe.tex", "manuscript/tables/", 
    recursive = TRUE)


