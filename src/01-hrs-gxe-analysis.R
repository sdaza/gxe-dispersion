# HRS GxE analysis
# author: sebastian daza


# setup
source("src/utils.R")
colors =c("#e34a33", "#2b8cbe")

# read data
dat = data.table(read_stata("data/hrs-obesity-pgs.dta"))
setnames(dat, names(dat), tolower(names(dat)))
names(dat)

# select only one record per individual
setorder(dat, hhidpn, wave)
dat[, seq := 1:.N, hhidpn]
dat = dat[seq == 1]
dat = dat[birthyr >= 1920 & birthyr <= 1960]
nrow(dat)

# define variables
table(dat$wave)
dat[, bmi_pgs := bmi_giant15]
dat[, zbmipgs := scale(bmi_pgs)]
dat[, zyear := scale(birthyr)]
dat[, zbmi := scale(srbmi)]
dat[, zage := scale(age)]
dat[, male := ifelse(gender == 1, 1, 0)]
dat[, white := ifelse(race == 1, 1, 0)]
countmis(dat)

# cohorts
table(dat$birthyr)

b = seq(1920, 1960, 4)
dat[, qbyr := cut(birthyr , b,
    labels = FALSE, include.lowest = TRUE)]    
t = dat[, .(min(birthyr), max(birthyr)), qbyr]
setorder(t, qbyr)

cohort_labels = paste0(t$V1, "-", t$V2)
cohort_labels = gsub("\\-19", "\\-", cohort_labels)
dat[, qbyr := factor(qbyr, labels = cohort_labels)]

# cohort standardization for correlations
dat[, czbmi := scale(srbmi), qbyr]
dat[, czbmipgs := scale(bmi_pgs), qbyr]

# complete cases (listwise deletion)
vars = c("zage", "male", "white", "zbmipgs", "zbmi")
cc = complete.cases(dat[, ..vars])
table(cc)
dat = dat[cc]
table(dat$qbyr)

# slope vs correlation plots using Bayesian models
f = bf(zbmi ~ male + white + zage + zbmipgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid) + (zbmipgs|qbyr))
m1 = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.99))
summary(m1)

coeff_m1 = brmsVaryingCoefPlot(m1, "b_zbmipgs", "r_qbyr[cohort,]", return_data = TRUE)

f = bf(czbmi ~ male + white + zage + czbmipgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid) + (czbmipgs|qbyr))
m2 = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    control = list(adapt_delta = 0.99))
summary(m2)

coeff_m2 = brmsVaryingCoefPlot(m2, "b_czbmipgs", "r_qbyr[cohort,]", return_data = TRUE)

# plot slopes vs correlation
coeff_m1$type = "Slope"
coeff_m2$type = "Correlation"
coeffs = rbind(coeff_m1, coeff_m2)
coeffs
savepdf("output/plots/hrs-slope-cor", 16, 12)
print(ggplot(coeffs, aes(cohort, median, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_point(size = 0.5) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha =  0.15, linetype = 0) +
    labs(
        # title = "Correlation vs Slope", subtitle = "HRS", 
        x = "\nBirth cohort", y = "Coefficient\n", 
        caption = "Random coefficient model adjusting for PCAs, gender, race, age first interview, and household random effects") +
    theme_minimal() +
    scale_color_manual(values = colors) + 
    theme() +
    theme(axis.text.x = element_text(angle = 0.0, size = 7, hjust = 0.5, vjust = 0.1), 
        legend.position="top", legend.title=element_blank())
    #scale_x_continuous(breaks = 1:10) + 
    # scale_y_continuous(breaks = seq(-1, 2.0, 0.3)) + 
    # geom_hline(yintercept = 0, size=0.5, color='red', alpha=0.8, linetype = 'dotted')) + 

)
dev.off()

# scatter plots
# slope and correlation plots
pp = list()
v = 1
for (i in cohort_labels) { 
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

savepdf("output/plots/hrs-slope-cor-scatter", 25, 30)
print(wrap_plots(pp, ncol = 3))
dev.off()
file.copy("output/plots/hrs-slope-cor-scatter.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    

# domingue's test
# using standardized measures to obtain convergence
dat[, E := zyear]
dat[, g := zbmipgs]
dat[, y := zbmi]
df = as.data.frame(dat[, .(E, y, g)])

# residuals
f = bf(zbmi ~ male + white + zage + 
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid))
t = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.90))
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
f = bf(zbmi ~ male + white + zage + zbmipgs + zbmipgs * zyear + 
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid))
md1 = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.99))

f = bf(zbmi ~ male + white + zage + zbmipgs + zbmipgs * zyear + 
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e, 
     sigma ~ zyear  + zbmipgs + male + white)
md2 = brm(f, data = dat, cores = 4, backend = "cmdstanr", 
    # threads = threading(8), 
    control = list(adapt_delta = 0.99))

tabs = readRDS("output/data/tabs.rds")
tabs[[1]] = extractBRMS(md1)
saveRDS(tabs, "output/data/tabs.rds")
