# ELSA GxE analysis
# author: sebastian daza


# setup
source("src//utils.R")
colors =c("#e34a33", "#2b8cbe")

# read data
dat = readRDS("data/elsa-renamed-data.rds")
names(dat)
dat[is.na(wave), wave := 0]
table(dat$birth_year)
dat = dat[birth_year >= 1920 & birth_year <= 1960]
table(dat$wave)
nrow(dat)
table(dat$proxy)

# select only one record per individual
setorder(dat, target_id, waveid)
dat[, seq := 1:.N, target_id]
vars = c("age", "sex", "white", "bmi_pgs", paste0("pc", 1:10))
dat[, (vars) := lapply(.SD, getFirstValue), target_id, .SDcols = vars]
dat[, srbmi := mean(bmi, na.rm = TRUE), target_id]
dat = dat[seq == 1]
nrow(dat)

# cohorts
table(dat$birth_year)
b = seq(1920, 1960, 4)
dat[, qbyr := cut(birth_year , b,
    labels = FALSE, include.lowest = TRUE)]   
t = dat[, .(min(birth_year), max(birth_year)), qbyr]
setorder(t, qbyr)

cohort_labels = paste0(t$V1, "-", t$V2)
cohort_labels = gsub("\\-19", "\\-", cohort_labels)
dat[, qbyr := factor(qbyr, labels = cohort_labels)] 

table(dat$qbyr)
dat[, zage := scale(age)]
dat[, zyear := scale(birth_year)]
dat[, zbmi := scale(srbmi)]
dat[, zbmipgs := scale(bmi_pgs)]
dat[, czbmipgs := scale(bmi_pgs), qbyr]
dat[, czbmi := scale(srbmi), qbyr]

# complete cases (listwise deletion)
anyDuplicated(dat$target_id)
vars = c("zbmipgs", "zbmi", "zbmipgs", "male", "white", "zage", paste0("pc", 1:10))
cc = complete.cases(dat[, ..vars])
table(cc)
dat = dat[cc]
# remove last cohort (too small sample size)
dat = dat[qbyr != "1957-60"]

# slope vs correlation plots using Bayesian models
f = bf(zbmi ~ male + zage + zbmipgs + 
     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (zbmipgs|qbyr) + (1|household_id))
m1 = brm(f, data = dat, backend = "cmdstanr", 
    cores = 4, 
    control = list(adapt_delta = 0.99))
summary(m1)

coeff_m1 = brmsVaryingCoefPlot(m1, "b_zbmipgs", "r_qbyr[cohort,]", return_data = TRUE)

f = bf(czbmi ~ male + zage + czbmipgs + 
     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (czbmipgs|qbyr) + (1|household_id))
m2 = brm(f, data = dat, backend = "cmdstanr", cores = 4, 
    control = list(adapt_delta = 0.99))
summary(m2)

coeff_m2 = brmsVaryingCoefPlot(m2, "b_czbmipgs", "r_qbyr[cohort,]", return_data = TRUE)

# plot slopes vs correlation
colors =c("#e34a33", "#2b8cbe")

coeff_m1$type = "Slope"
coeff_m2$type = "Correlation"
coeffs = rbind(coeff_m1, coeff_m2)
coeffs

savepdf("output/plots/elsa-slope-cor", 16, 12)
print(ggplot(coeffs, aes(cohort, median, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_point(size = 0.5) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha =  0.15, linetype = 0) +
    labs(
        # title = "Correlation vs Slope", subtitle = "ELSA", 
        x = "\nBirth cohort", y = "Coefficient\n", 
        caption = "Random coefficient model adjusting for PCAs, gender, age first interview, and household random effects") +
    theme_minimal() +
    scale_color_manual(values = colors) + 
    theme() +
    theme(axis.text.x = element_text(angle = 0.0, size = 7, hjust = 0.5, vjust = 0.1), 
        legend.position="top", legend.title=element_blank())
)
dev.off()

# scatter plots
pp = list()
v = 1
for (i in cohort_labels[-length(cohort_labels)]) { 
    temp = copy(dat[qbyr == i, .(zbmi, zbmipgs)])
    slope = specify_decimal(coef(lm(zbmi ~ zbmipgs, data = temp))[2], 2)
    corr = specify_decimal(cor(temp[, .(zbmi, zbmipgs)])[1, 2], 2)
    pp[[v]] = ggplot(temp, aes(zbmipgs, zbmi)) +
        geom_point(color = "#2b8cbe", size = 0.5, alpha = 0.2) +
        geom_smooth(method = "lm", color = "#e34a33", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("Cohort ", i), 
            subtitle = TeX(paste0("$\\rho$=", corr, ", $\\beta$=", slope)), 
            x = "BMI polygenic risk score", y = "BMI") +
        theme_minimal()
    v = v + 1;
}

savepdf("output/plots/elsa-slope-cor-scatter", 25, 30)
print(wrap_plots(pp, ncol = 3))
dev.off()
file.copy("output/plots/elsa-slope-cor-scatter.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    


# domingue's test
# using standardized measures for convergence
dat[, E := zyear]
dat[, g := zbmipgs]
dat[, y := zbmi]

# get residuals
df = as.data.frame(dat[, .(E, y, g)])
bf(y ~ male + zage + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 +
    (1|household_id))
t = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.99))
res = residuals(t)

# simple
est = mlest(df,hess=TRUE)
xi.test(est)

# y and yr are highly correlated
y = df$y
cor(y, res[, 1])

# residuals
df$y = res[, 1]
est = mlest(df,hess=TRUE)
xi.test(est)

# distributional model tables (standardized)
f = bf(zbmi ~ male + zage + zbmipgs + zbmipgs * zyear +
     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1|household_id), 
     sigma ~ zyear)
md1 = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.99))
summary(md1)

tabs = readRDS("output/data/tabs.rds")
tabs[[2]] = extractBRMS(md1)
saveRDS(tabs, "output/data/tabs.rds")

# create tables elsa and hrs
cnames = c("HRS", "ELSA")
custom_coeff_map = list(
    "Intercept" = "Intercept",
    "zyear" = "Birth year (z-score)",
    "zbmipgs" = "BMI polygenic risk score",
    "zbmipgs:zyear" = "BMI polygenic risk score x Birth year (z-scores)", 
    "sigma_Intercept" = "sIntercept", 
    "sigma_zyear" = "sBirth year (z-score)", 
    "sd(Intercept)" = "$\\sigma$ Intercept Households"
)

screenreg(tabs, custom.coef.map = custom_coeff_map)
caption = paste0("BMI and birth year, HRS-ELSA")

groups = list("Fixed effects" = 1:4,  "Residuals" = 5:6, "Random effects" = 7)

texreg::texreg(tabs, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:dist-model",
    groups = groups,
    custom.note = "\\item  $^*$ Null hypothesis value outside the confidence interval. HH = household. 
        \\item All the models adjust for gender, race, age first interview, population stratification (10 PCs), 
        and household random effects",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0("output/tables/hrs-elsa-dist-models.tex")
)    
file.copy("output/tables/hrs-elsa-dist-models.tex", "manuscript/tables/", 
    recursive = TRUE)