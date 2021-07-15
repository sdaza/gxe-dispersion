# explore HRS gxe analysis
# author: sebastian daza


library(haven)
library(data.table)
library(texreg)
library(brms)
library(tidybayes)
library(ggplot2)
library(patchwork)
source("src/utils.R")

# read data
dat = data.table(read_stata("data/hrs-obesity-pgs.dta"))
setnames(dat, names(dat), tolower(names(dat)))
names(dat)

# select only one record per individual
dat[, seq := 1:.N, hhidpn]
dat[, n := NULL]
dat = dat[seq == 1]
dat[, bmi_pgs := bmi_giant15]
dat[, zyear := scale(birthyr)]
dat = dat[birthyr >= 1920 & birthyr <= 1960]
dat[, qbyr := cut(birthyr , quantile(birthyr, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]

names(dat)
dim(dat)

# slope and correlation plots
plots = list()
for (i in 1:10) {
    slope = specify_decimal(coef(lm(srbmi ~ bmi_pgs, data = dat[qbyr == i]))[2], 2)
    corr = specify_decimal(cor(dat[qbyr == i, .(srbmi, bmi_pgs)])[1, 2], 2)
    plots[[i]] = ggplot(dat[qbyr == i], aes(bmi_pgs, srbmi)) + geom_point(size = 0.5, 
        alpha = 0.1) +
        geom_smooth(method = "lm", color = "red", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("Cohort ", i), 
            subtitle = paste0("Correlation: ", corr, "; Slope: ", slope), 
            x = "BMI polygenic score", y "Self-reported BMI") +
        theme_minimal()
}

savepdf("output/plots/hrs-slope-cor-gxe", 25, 20)
print(wrap_plots(plots, ncol = 3))
dev.off()
file.copy("output/plots/hrs-slope-cor-gxe.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    


m1 = lm(srbmi ~ gender + bmi_pgs, data = dat)
summary(m1)
# m2 = lm(srbmi ~ gender + bmi_pgs * zyear, data = dat)

f = bf(srbmi ~ gender + bmi_pgs)

f = bf(srbmi ~ gender + age + bmi_pgs * zyear, sigma ~ zyear)
m1 = brm(f, data = dat, backend = "cmdstanr", cores = 4)

f = bf(srbmi ~ gender + age + bmi_pgs * zyear +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d +  pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d +  pc6_10e, sigma ~ zyear)
m1 = brm(f, data = dat, backend = "cmdstanr", cores = 4)
summary(m2)

f = bf(srbmi ~ gender + bmi_pgs + (1 + bmi_pgs|birthyr), sigma ~ birthyr)
m3 = brm(srbmi ~ gender + bmi_pgs + (1 + bmi_pgs|birthyr), data = dat, 
    backend = "cmdstanr", cores = 4)

summary(m3)
screenreg(list(m1, m2))
str(dat)
dim(dat)


s = data.table(spread_draws(m3, r_birthyr[bmi_pgs,term], b_bmi_pgs))
s = s[term == "bmi_pgs"]
setnames(s, "bmi_pgs", "byear")

s[, slope := r_birthyr + b_bmi_pgs]
ss = s[, .(m = median(slope), 
    l = quantile(slope, probs = 0.025), 
    h = quantile(slope, probs = 0.975)), byear]

savepdf("testing")
ggplot(ss, aes(byear, m)) + geom_line(color='#2b8cbe', size = 0.4) +
    geom_ribbon(aes(ymin = l, ymax = h), fill = '#a6bddb', alpha=0.2) + 
    theme_minimal() + 
    labs(y = "BMI PGS slope", x = "Birth year")
dev.off()

