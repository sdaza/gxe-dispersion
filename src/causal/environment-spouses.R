#######################################################
# exploring casual model environament spouses and GxE
#######################################################


library(data.table)
library(brms)
library(performance)

dat = fread("models/BMI-spouses-GxE/output/bmi_0_0.csv")


summary(dat$bmi)
summary(dat$bmi_spouse)

# hist(dat$bmi)
# hist(dat$bmi_spouse)
# hist(dat$pgs)
# hist(dat$pgs_spouse)

cor(dat$bmi, dat$bmi_spouse)
cor(dat$pgs, dat$pgs_spouse)
cor(dat$environment, dat$environment_spouse)

dat[, couple := paste0(sort(.SD), collapse = ""), .SDcols = c("id", "spouse_id"), 
    by = 1:nrow(dat)]
setorder(dat, couple)


m0 = brm(bmi ~ pgs + environment + pgs * environment, 
    data = dat, iter = 2000, chains = 1)
summary(m0)

m1 = brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
    data = dat, iter = 2000, chains = 1)
summary(m1)


icc(m1, by_group = FALSE)
