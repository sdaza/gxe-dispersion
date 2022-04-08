# exploratory ELSA gxe analysis
# author: sebastian daza


library(data.table)
library(texreg)
library(brms)
# library(tidybayes)
library(ggplot2)
library(patchwork)
source("src//utils.R")

# read data
dat = readRDS("data/elsa-renamed-data.rds")
names(dat)
dat[is.na(wave), wave := 0]
table(dat$wave)
nrow(dat)
table(dat$proxy)

# bmi
summary(dat$bmi)

table(dat$birth_year)
dat = dat[birth_year >= 1920 & birth_year <= 1960]
summary(dat$bmi)
# select only one record per individual
setorder(dat, target_id, waveid)
dat[, seq := 1:.N, target_id]
vars = c("age", "sex", "bmi_pgs", paste0("pc", 1:10))
dat[, (vars) := lapply(.SD, getFirstValue), target_id, .SDcols = vars]
dat[, srbmi := mean(bmi, na.rm = TRUE), target_id]

dat = dat[seq == 1]
dim(dat)

dat[, .(target_id, bmi, srbmi)]
dat[, qbyr := cut(birth_year , quantile(birth_year, probs = 0:10/10),
    labels = FALSE, include.lowest = TRUE)]

table(dat$birth_year)
dat[, .(qbyr, birth_year)]
tt = dat[, .(.N, sd = sd(srbmi, na.rm = TRUE)), birth_year]
setorder(tt, birth_year)
ggplot(tt, aes(birth_year, sd)) + geom_line() + theme_minimal()

dat[, zage := scale(age)]
dat[, zyear := scale(birth_year)]
dat[, male := ifelse(sex == 1, 1, 0)]
dat[, bmipgs := scale(bmi_pgs)]
dat[, zbmipgs := scale(bmipgs), qbyr]
dat[, zbmi := scale(srbmi), qbyr]
dat[, zsrbmi := scale(srbmi)]

anyDuplicated(dat$target_id)
table(dat$male)
summary(dat$zage)
vars = c("bmipgs", "zbmi", "zbmipgs", "male", "zage", "srbmi", paste0("pc", 1:10))
table(complete.cases(dat[, ..vars]))

# create plots using brms
countmis(dat,vars)

f = bf(srbmi ~ male + zage + bmipgs + 
     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (bmipgs|qbyr))
m1 = brm(f, data = dat, backend = "cmdstanr", cores = 4, 
    control = list(adapt_delta = 0.9))
summary(m1)

f = bf(zbmi ~ male + zage + zbmipgs + 
     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (zbmipgs|qbyr))
m2 = brm(f, data = dat, backend = "cmdstanr", cores = 4, 
    control = list(adapt_delta = 0.9))
summary(m2)

getRandomCoeff = function(model, random, coef) {
    tab = data.table(ranef(model)[[random]][, , coef], keep.rownames = TRUE)
    ff = fixef(model)[coef, "Estimate"]
    setnames(tab, names(tab), c("rn", "est", "se", "lo", "up"))
    tab[, `:=` ( est = ff + est, lo = lo + ff, up = up + ff)]
    return(tab)
}

ss = getRandomCoeff(m1, "qbyr", "bmipgs")
ss[, type := factor("slope")]
cc = getRandomCoeff(m2, "qbyr", "zbmipgs")
cc[, type := factor("correlation")]

tab = rbind(ss, cc)
tab[, rn := as.numeric(rn)]
tab

savepdf("output/plots/elsa-slope-cor")
print(ggplot(tab, aes(rn, est, group = type, color = type, fill = type)) + 
    geom_line(size = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = up), alpha =  0.2, linetype = 0) +
    labs(title = "Slope vs Correlation", subtitle = "ELSA", x = "\nBirth cohort", y = "Coefficient/correlation\n", 
        caption = "Random coefficient model adjusting for PCAs, gender, age first interview") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle = 0, hjust = 0.1, vjust = 0.1)) + 
    scale_fill_brewer(palette="Dark2") +
    scale_color_brewer(palette="Dark2") +
    scale_x_continuous(breaks = 1:10) + 
    scale_y_continuous(breaks = seq(-1, 2.0, 0.2)) + 
    # geom_hline(yintercept = 0, size=0.5, color='red', alpha=0.8, linetype = 'dotted')) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
)
dev.off()

file.copy("output/plots/elsa-slope-cor.pdf", 
    "manuscript/plots/", 
    recursive = TRUE)    

# models

# unstandardized
f = bf(srbmi ~ male + zage + bmipgs * zyear + 
    pc1 + pc2 + pc3 + pc4 + pc5 + 
    pc6 + pc7 + pc8 + pc9 + pc10 + (1|household_id), 
    sigma ~ zyear)
m1 = brm(f, data = dat, backend = "cmdstanr", cores = 4, 
    control = list(adapt_delta = 0.9))
summary(m1)

tabs = list()
tabs[[1]] = extractBRMS(m1)

# standardized
f = bf(zsrbmi ~ male + zage + bmipgs * zyear +
    pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + 
    pc9 + pc10+ (1|household_id), 
    sigma ~ zyear)
m2 = brm(f, data = dat, backend = "cmdstanr", cores = 4, 
    control = list(adapt_delta = 0.9))
summary(m2)
tabs[[2]] = extractBRMS(m2)

screenreg(tabs)

# tables
cnames = c("Unstandardized BMI", "Standardized BMI")
custom_coeff_map = list(
    "Intercept" = "Intercept",
    "bmipgs:zyear" = "BMI x Birth year (z-score)", 
    "sigma_Intercept" = "sIntercept", 
    "sigma_zyear" = "Birth year (z-score)", 
    "sd(Intercept)" = "$\\sigma$ Intercept HH"
)

screenreg(tabs, custom.coef.map = custom_coeff_map)
caption = paste0("BMI and birth year, ELSA")

groups = list("Fixed effects" = 1:2,  "Residuals" = 3:4, "Random effects" = 5)

texreg::texreg(tabs, 
    caption = caption, 
    custom.model.names = cnames, 
    custom.coef.map = custom_coeff_map,
    label = "tab:elsa-model",
    groups = groups,
    custom.note = "\\item  $^*$ Null hypothesis value outside the confidence interval. HH = household. 
        \\item All the models adjust for gender, age first interview, population stratification (10 PCs)",
    scalebox = 0.7,
    center = TRUE,
    dcolumn = TRUE, 
    use.packages = FALSE, 
    threeparttable = TRUE, 
    caption.above = TRUE, 
    file = paste0("output/tables/elsa-model.tex")
)    
file.copy("output/tables/elsa-model.tex", "manuscript/tables/", 
    recursive = TRUE)
