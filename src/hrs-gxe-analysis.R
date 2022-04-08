# explore HRS gxe analysis
# author: sebastian daza


library(haven)
library(data.table)
library(texreg)
library(brms)
library(tidybayes)
library(ggplot2)
library(patchwork)
# library(mice)
source("src/utils.R")

# read data
dat = data.table(read_stata("data/hrs-obesity-pgs.dta"))
setnames(dat, names(dat), tolower(names(dat)))
names(dat)


# select only one record per individual
setorder(dat, hhidpn, wave)
dat[, seq := 1:.N, hhidpn]
dat = dat[seq == 1]

table(dat$wave)

dat[, bmi_pgs := bmi_giant15]
dat[, zbmipgs := scale(bmi_pgs)]
dat = dat[birthyr >= 1920 & birthyr <= 1960]
dat[, zyear := scale(birthyr)]

# self-reported BMI
dat[, zbmi := scale(srbmi)]
dat[, zage := scale(age)]
dat[, male := ifelse(gender == 1, 1, 0)]
dat[, white := ifelse(race == 1, 1, 0)]

countmis(dat)

# definition of cohorts
table(dat$birthyr)
b = seq(1920, 1960, 4)
dat[, qbyr := cut(birthyr , b,
    labels = FALSE, include.lowest = TRUE)]    
t = dat[, .(min(birthyr), max(birthyr)), qbyr]
setorder(t, qbyr)
t

# cohort standardization
dat[, czbmi := scale(srbmi), qbyr]
dat[, czbmipgs := scale(bmi_pgs), qbyr]


t = tidybayes::spread_draws(m1, b_zbmipgs, r_qbyr[cohort,])

# initial model 
f = bf(zbmi ~ male + white + zage + zbmipgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid) + (zbmipgs|qbyr))
m1 = brm(f, data = dat, cores = 4, backend = "cmdstanr")
summary(m1)

coeff_m1 = varyingCoefPlot(m1, "b_zbmipgs", "r_qbyr[zbmipgs,]", 
    file = "testing", return_data = TRUE)


f = bf(czbmi ~ male + white + zage + czbmipgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid) + (czbmipgs|qbyr))
m2 = brm(f, data = dat, cores = 4, backend = "cmdstanr")
summary(m2)

coeff_m2 = varyingCoefPlot(m2, "b_czbmipgs", "r_qbyr[czbmipgs,]", 
    file = "testing2", return_data = TRUE)

a = data.table(coeff_m1)
a[, type := "slope"][, cohort := zbmipgs]

dat[, E := zyear]
dat[, g := zpgsbmi]
dat[, y := zbmi]
df = as.data.frame(dat[, .(E, y, g)])
library(scalingGxE)
est = mlest(df,hess=TRUE)
xi.test(est)

# distributional model 
f = bf(y ~ g*E, 
     sigma ~ E)
m1 = brm(f, data = df, backend = "cmdstanr", cores = 4)
summary(m1)
tt = scalingTest(m1)
plotScalingTest(tt)
plotDecomp(m1)

df$group = rbinom(nrow(df), 5, 0.5)
table(df$group)

f = bf(y ~ g*E + (g|group))

m1 = brm(f, data = df, backend = "cmdstanr", cores = 4)
summary(m1)


get_variables(m1)


dat[, qbyr := cut(birthyr , quantile(birthyr, probs = 0:10/10),
        labels = FALSE, include.lowest = TRUE)]
dat[, cbirthyr := scale(birthyr, scale = FALSE)]

dim(dat)
countmis(dat)

# correlation spouses
dat[, id := 1:.N, hhid]
testing = dcast(dat, hhid ~ id, 
    value.var = c("rbmi", "bmi_pgs", "raeduc", "rshlt", "rmbmi"))
testing

cor(testing[, .(rbmi_1, rbmi_2)], use = "complete.obs")
cor(testing[, .(rmbmi_1, rmbmi_2)], use = "complete.obs")
cor(testing[, .(bmi_pgs_1, bmi_pgs_2)], use = "complete.obs")
cor(testing[, .(raeduc_1, raeduc_2)], use = "complete.obs")
cor(testing[, .(rshlt_1, rshlt_2)], use = "complete.obs")
cor(testing[, .(rbmi_1, raeduc_1)], use = "complete.obs")

table()
summary(testing)

dcast()
names(dat)
t = lm(srbmi ~ birthyr, data = dat)
dat[, resd := resid(t)]
savepdf("testing")
ggplot(dat, aes(age, resd)) + geom_point()
dev.off()

names(dat)
# dat = dat[!is.na(pmbmi)]
# dat[, zpmbmi := scale(pmbmi)]




test = dcast(dat, hhidpn ~ age, value.var = "birthyr")
test[, apply(.SD, 2, function(x) sum(!is.na(x)))]
names(dat)
dim(dat)

tdat = copy(dat)
tdat = tdat[age >=50 & age <= 70]
table(tdat$birthyr)

table(tdat[, .(birthyr, age)])

# slope and correlation plots
plots = list()
for (i in 1:10) {
    slope = specify_decimal(coef(lm(pmbmi ~ bmi_pgs, data = dat[qbyr == i]))[2], 2)
    corr = specify_decimal(cor(dat[qbyr == i, .(srbmi, bmi_pgs)])[1, 2], 2)
    plots[[i]] = ggplot(dat[qbyr == i], aes(bmi_pgs, srbmi)) + geom_point(size = 0.5, 
        alpha = 0.1) +
        geom_smooth(method = "lm", color = "red", alpha = 0.2, size = 0.3) + 
        labs(title = paste0("Cohort ", i), 
            subtitle = paste0("Correlation: ", corr, "; Slope: ", slope), 
            x = "BMI polygenic score", y = "BMI") +
        theme_minimal()
}

savepdf("output/plots/hrs-slope-cor-gxe", 25, 20)
print(wrap_plots(plots, ncol = 3))
dev.off()
file.copy("output/plots/hrs-slope-cor-gxe.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    



# unstandardized
f = bf(pmbmi ~ gender + race + zage + bmi_pgs * zyear +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid), 
     sigma ~ zyear)
m2 = brm(f, data = dat, backend = "cmdstanr", cores = 4)
summary(m2)

tabs = list()
tabs[[1]] = extractBRMS(m2)


# standardized
f = bf(zpmbmi ~ gender + race + zage + bmi_pgs * zyear +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + (1|hhid), 
     sigma ~ zyear)
m3 = brm(f, data = dat, backend = "cmdstanr", cores = 4)
summary(m3)
tabs[[2]] = extractBRMS(m3)

screenreg(tabs)

# tables
cnames = c("Unstandardized BMI", "Standardized BMI")
custom_coeff_map = list(
    "Intercept" = "Intercept",
    "bmi_pgs:zyear" = "BMI x Birth year (z-score)", 
    "sigma_Intercept" = "sIntercept", 
    "sigma_zyear" = "Birth year (z-score)", 
    "sd(Intercept)" = "$\\sigma$ Intercept HH"
)

screenreg(tabs, custom.coef.map = custom_coeff_map)
caption = paste0("BMI and birth year, HRS")

groups = list("Fixed effects" = 1:2,  "Residuals" = 3:4, "Random effects" = 5)

texreg::texreg(tabs, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:hrs-model",
    groups = groups,
    custom.note = "\\item  $^*$ Null hypothesis value outside the confidence interval. HH = household. 
        \\item All the models adjust for gender, race, age first interview, population stratification (10 PCs)",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0("output/tables/hrs-model.tex")
)    
file.copy("output/tables/hrs-model.tex", "manuscript/tables/", 
    recursive = TRUE)

# compute R2
s = data.table(spread_draws(m3, r_birthyr[bmi_pgs,term], b_bmi_pgs))
s = s[term == "bmi_pgs"]
setnames(s, "bmi_pgs", "byear")
s[, slope := r_birthyr + b_bmi_pgs]
ss = s[, .(m = median(slope), 
    l = quantile(slope, probs = 0.025), 
    h = quantile(slope, probs = 0.975)), byear]

savepdf("output/plots/hrs-slope")
ggplot(ss, aes(byear, m)) + geom_line(color='#2b8cbe', size = 0.4) +
    geom_ribbon(aes(ymin = l, ymax = h), fill = '#a6bddb', alpha=0.2) + 
    theme_minimal() +
     
    labs(y = "BMI PGS slope", x = "Birth year")
dev.off()
file.copy("output/plots/hrs-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    

brms::bayes_R2(m3)
ypred = posterior_epred(m3)
bayes_r2(dat$srbmi, ypred)
r2_m3 = bayes_r2_group(dat$srbmi, ypred, dat$birthyr)
savepdf("output/plots/hrs-r2")
ggplot(r2_m4, aes(group, m)) + geom_line(color='#2b8cbe', size = 0.4) +
    geom_ribbon(aes(ymin = l, ymax = h), fill = '#a6bddb', alpha=0.2) + 
    theme_minimal() + 
    labs(y = "R2", x = "Birth year")
dev.off()
file.copy("output/plots/hrs-r2.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)   


lambda0 = 0.2
lambda1 = 0.5
pi0 = 0.5
pi1 = 0.2

pi0 / lambda0
pi1 / lambda1

% pi0 = G
% pi1 = EG
% lambda0 = intercept sigma
% labdaa1 = effect sigma
% pi0/lambda0 = pi1 / lambda1
%  0.5/0.4 = 0.3 / 0.24

# additional testing models
# ndat = copy(dat)
# ndat[, sq := 1:.N, hhid]
# ndat = ndat[sq == 1]

f = bf(srbmi ~ gender + race + zage + bmi_pgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + 
     (1|hhid),
     sigma ~ (1|birthyr))
m0 = brm(f, data = dat, backend = "cmdstanr", cores = 4)
summary(m0)

f = bf(srbmi ~ gender + race + zage + bmi_pgs +
     pc1_5a + pc1_5b + pc1_5c + pc1_5d + pc1_5e + 
     pc6_10a + pc6_10b + pc6_10c + pc6_10d + pc6_10e + 
     (1|hhid) + (bmi_pgs|birthyr), 
     sigma ~ (1|birthyr))
m1 = brm(f, data = dat, backend = "cmdstanr", cores = 4)
summary(m1)

newdata = dat[birthyr == 1959 & gender == 1 & race == 1]
table(newdata$age)
newdata = newdata[age == 51]

p1 = array(predict(m1, newdata, summary = FALSE))
quantile(p1, probs = c(0.5, 0.025, 0.975))

newdata[, birthyr := 1921]
p2 = array(predict(m1, newdata, summary = FALSE))
quantile(p2, probs = c(0.5, 0.025, 0.975))


r2_0 = brms::bayes_R2(m0)
r2_1 = brms::bayes_R2(m1)
loo0 = loo(m0)
loo1 = loo(m1)

# test = dat[, quantile(age, probs = c(0.5)), birthyr]
# setorder(test, birthyr)
# test
r2_0
r2_1

comp = loo_compare(loo0, loo1)
print(comp, simplify = FALSE, digits = 3)

s = data.table(spread_draws(m1, r_birthyr[bmi_pgs,term], b_bmi_pgs))
s = s[term == "bmi_pgs"]
setnames(s, "bmi_pgs", "byear")
s[, slope := r_birthyr + b_bmi_pgs]
ss = s[, .(m = median(slope), 
    l = quantile(slope, probs = 0.025), 
    h = quantile(slope, probs = 0.975)), 
    byear]

savepdf("output/plots/hrs-slope")
ggplot(ss, aes(byear, m)) + geom_line(color='#2b8cbe', size = 0.4) +
    geom_ribbon(aes(ymin = l, ymax = h), fill = '#a6bddb', alpha=0.2) + 
    theme_minimal() +     
    labs(y = "BMI PGS slope", x = "Birth year")
dev.off()
file.copy("output/plots/hrs-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    


s = data.table(spread_draws(m0, r_birthyr[zage,term], b_zage))
s = s[term == "zage"]
setnames(s, "zage", "byear")
s[, slope := r_birthyr + b_zage]
ss = s[, .(m = median(slope), 
    l = quantile(slope, probs = 0.025), 
    h = quantile(slope, probs = 0.975)), 
    byear]

savepdf("output/plots/hrs-slope")
ggplot(ss, aes(byear, m)) + geom_line(color='#2b8cbe', size = 0.4) +
    geom_ribbon(aes(ymin = l, ymax = h), fill = '#a6bddb', alpha=0.2) + 
    theme_minimal() +     
    labs(y = "BMI PGS slope", x = "Birth year")
dev.off()
file.copy("output/plots/hrs-slope.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    
