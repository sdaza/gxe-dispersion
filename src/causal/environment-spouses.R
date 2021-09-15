#######################################################
# exploring casual model environament spouses and GxE
#######################################################


library(data.table)
library(brms)
library(texreg)
library(foreach)
library(doParallel)

path = "models/BMI-spouses-GxE/output/"
clusters = 3

# remove files
# files = paste0(path, list.files(path, pattern = "*.csv"))
# sapply(files, unlink, recursive = TRUE)

files = paste0(path, list.files(path, pattern = "bmi"))
dat = rbindlist(lapply(files, fread))

files = paste0(path, list.files(path, pattern = "parameter"))
params = rbindlist(lapply(files, fread))

# dat[, couple := paste0(sort(.SD), collapse = ""), 
    # .SDcols = c("id", "spouse_id"), 
    # by = .(1:nrow(dat))]
# setorder(dat, couple)

table(dat$iteration) 
table(dat$replicate)

cor(dat$bmi, dat$bmi_spouse)
cor(dat$pgs, dat$pgs_spouse)
table(dat$susceptability)

cl = makeCluster(clusters)
registerDoParallel(cl)

output = foreach(i = 1:max(params$replicate),
        .packages = c("brms", "cmdstanr", "doParallel", 
            "data.table", "rstan", "rethinking")) %dopar% {
    
    options(mc.cores = 1)
    temp = copy(dat[iteration == 1 & replicate == i])
    brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
        data = temp, iter = 2000, chains = 1, backend = "cmdstan")
}
stopCluster(cl)

baseline = combine_models(mlist = output, check_data = FALSE)

summary(baseline)
params

m2 = brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
    data = dat, iter = 2000, chains = 1, backend = "cmdstan")
summary(m2)
 

m0 = brm(bmi ~ pgs + environment + pgs * environment, 
    data = dat, iter = 2000, chains = 1, backend = "cmdstan")
summary(m0)

m1 = brm(bmi ~ pgs + environment + pgs * environment + unobserved +
    pgs * unobserved, 
    data = dat, iter = 2000, chains = 1, backend = "cmdstan")
summary(m1)

m2 = brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
    data = dat, iter = 2000, chains = 1, backend = "cmdstan")
summary(m2)
 

screenreg(list(m0, m1))
params

# hist(dat$bmi)
# hist(dat$bmi_spouse)
# hist(dat$pgs)
# hist(dat$pgs_spouse)

models = list()
for (i in 1:100) {
    temp = dat[iteration == 1 & replicate == i]
    models[[i]] =  brm(bmi ~ pgs + environment + pgs * environment, 
    data = temp, iter = 2000, chains = 1, backend = "cmdstan")

}
m1 = combine_models(mlist = models, check_data = FALSE)
summary(m1)


cor(dat$bmi, dat$bmi_spouse)
cor(dat$pgs, dat$pgs_spouse)
cor(dat$environment, dat$environment_spouse)


dat[, couple := paste0(sort(.SD), collapse = ""), .SDcols = c("id", "spouse_id"), 
    by = 1:nrow(dat)]
setorder(dat, couple)
dat

m0 = brm(bmi ~ pgs + environment + pgs * environment, 
    data = dat, iter = 2000, chains = 1)
summary(m0)

m1 = brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
    data = dat, iter = 2000, chains = 1)
summary(m1)



