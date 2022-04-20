# scenario simulation
# author: sebastian daza


# setup
source("src/utils.R")

set.seed(124907)
colors =c("#e34a33", "#2b8cbe")

# model names
cnames = c("S1", "S2", "S3")
mnames = c("scaling", "gxe", "hc")

# create data
nreplicates = 100
n = 10000
datasets = list(
    dts = replicate(nreplicates, simScaling(rnorm(n, 0, 1)), simplify = FALSE),
    dti = replicate(nreplicates, simInteraction(rnorm(n, 0, 1)), simplify = FALSE), 
    dth = replicate(nreplicates, simHC(rnorm(n, 0, 1)), simplify = FALSE)
)

# explore data
head(datasets[["dts"]][[sample(1:100, 1)]])
head(datasets[["dti"]][[sample(1:100, 1)]])
head(datasets[["dth"]][[sample(1:100, 1)]])
length(datasets[["dts"]])

# scaling
models = runModels(datasets[["dts"]], clusters = 10)
m1_dist = models[["distributional"]]

# plot correlation vs slope
coeff_m1_slope = brmsVaryingCoefPlot(models[["slope"]], "b_zg", "r_qE[cohort,]", return_data = TRUE)
coeff_m1_correlation = brmsVaryingCoefPlot(models[["correlation"]], "b_gzg", "r_qE[cohort,]", return_data = TRUE)
coeff_m1_slope$type = "Slope"
coeff_m1_correlation$type = "Correlation"
coeffs = rbind(coeff_m1_slope, coeff_m1_correlation)

savepdf("output/plots/scaling-cor-slope", 16, 12)
print(ggplot(coeffs, aes(cohort, median, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_point(size = 0.5) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha =  0.15, linetype = 0) +
    labs(
        # title = "Correlation vs Slope", subtitle = "Scenario 1", 
        x = "\nE", y = "Coefficient\n", 
        caption = "Random coefficient model") +
    theme_minimal() +
    scale_color_manual(values = colors) + 
    theme() +
    theme(axis.text.x = element_text(angle = 0.0, size = 7, hjust = 0.5, vjust = 0.1), 
        legend.position = "top", legend.title = element_blank())  + 
    scale_x_continuous(breaks = 1:10) 
)
dev.off()

rm(models, coeff_m1_slope, coeff_m1_correlation)

# GxE, heterocesdasticity
models = runModels(datasets[["dti"]], clusters = 10)
m2_dist = models[["distributional"]]

# plot slopes vs correlation
coeff_m2_slope = brmsVaryingCoefPlot(models[["slope"]], "b_zg", "r_qE[cohort,]", return_data = TRUE)
coeff_m2_correlation = brmsVaryingCoefPlot(models[["correlation"]], "b_gzg", "r_qE[cohort,]", return_data = TRUE)
coeff_m2_slope$type = "Slope"
coeff_m2_correlation$type = "Correlation"
coeffs = rbind(coeff_m2_slope, coeff_m2_correlation)

savepdf("output/plots/gxe-cor-slope", 16, 12)
print(ggplot(coeffs, aes(cohort, median, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_point(size = 0.5) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha =  0.15, linetype = 0) +
    labs(
        # title = "Correlation vs Slope", subtitle = "Scenario 2", 
        x = "\nE", y = "Coefficient\n", 
        caption = "Random coefficient model") +
    theme_minimal() +
    scale_color_manual(values = colors) + 
    theme() +
    theme(axis.text.x = element_text(angle = 0.0, size = 7, hjust = 0.5, vjust = 0.1), 
        legend.position = "top", legend.title = element_blank())  + 
    scale_x_continuous(breaks = 1:10) 
)
dev.off()

rm(models, coeff_m2_slope, coeff_m2_correlation)

# GxE, no scaling, heterocedasticity
models = runModels(datasets[["dth"]], clusters = 10)
m3_dist = models[["distributional"]]

# plot slopes vs correlation
coeff_m3_slope = brmsVaryingCoefPlot(models[["slope"]], "b_zg", "r_qE[cohort,]", return_data = TRUE)
coeff_m3_correlation = brmsVaryingCoefPlot(models[["correlation"]], "b_gzg", "r_qE[cohort,]", return_data = TRUE)
coeff_m3_slope$type = "Slope"
coeff_m3_correlation$type = "Correlation"
coeffs = rbind(coeff_m3_slope, coeff_m3_correlation)

savepdf("output/plots/hc-cor-slope", 16, 12)
print(ggplot(coeffs, aes(cohort, median, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_point(size = 0.5) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha =  0.15, linetype = 0) +
    labs(
        # title = "Correlation vs Slope", subtitle = "Scenario 3", 
        x = "\nE", y = "Coefficient\n", 
        caption = "Random coefficient model") +
    theme_minimal() +
    scale_color_manual(values = colors) + 
    theme() +
    theme(axis.text.x = element_text(angle = 0.0, size = 7, hjust = 0.5, vjust = 0.1), 
        legend.position = "top", legend.title = element_blank())  + 
    scale_x_continuous(breaks = 1:10) 
)
dev.off()

rm(models, coeff_m3_slope, coeff_m3_correlation)

# convergence testing 
models = list(m1_dist, m2_dist, m3_dist)
testConvergence(models)

# create table
cmap = list(
    "Intercept"  = "$\\tau_0$", 
    "E" = "$\\tau_1$", 
    "g" = "$\\pi_0$", 
    "g:E" = "$\\pi_1$", 
    "sigma_Intercept" = "$\\lambda_0$",
    "sigma_E" = "$\\lambda_1$")
 
suppressWarnings(
    texreg(models, 
        caption = "Distributional models",
        custom.coef.map = cmap,
        label = "tab:sim-dist-models",
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
        file = "output/tables/sim-dist-models.tex")
)
file.copy("output/tables/sim-dist-models.tex", "manuscript/tables", 
    recursive = TRUE)

# scaling test and plots
for (i in seq_along(models)) {
    t = scalingTest(models[[i]])
    plotScalingTest(t, paste0("output/plots/", mnames[i], "-test"))
    file.copy(paste0("output/plots/", mnames[i], "-test.pdf"), "manuscript/plots", 
        recursive = TRUE)
    plotDecomp(models[[i]], paste0("output/plots/", mnames[i], "-decomp"))
    file.copy(paste0("output/plots/", mnames[i], "-decomp.pdf"), "manuscript/plots", 
        recursive = TRUE)
}

# scatter correlation vs slope plots
s = 60

# scaling 
temp = datasets[["dts"]][[s]]
slopeCorPlot(temp, "qE", "output/plots/scaling-scatter-cor-slope")
file.copy("output/plots/scaling-scatter-cor-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)  

# gxe, heterocedasticity
temp = data.table(datasets[["dti"]][[s]])
slopeCorPlot(temp, "qE", "output/plots/gxe-scatter-cor-slope")
file.copy("output/plots/gxe-scatter-cor-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)  

# GxE, no scaling, heterocedasticity
temp = data.table(datasets[["dth"]][[s]])
slopeCorPlot(temp, "qE", "output/plots/hc-scatter-cor-slope")
file.copy("output/plots/hc-scatter-cor-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE) 

# save simulated data
saveRDS(datasets, "output/data/simulated-data.rds", )