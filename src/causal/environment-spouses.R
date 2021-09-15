#######################################################
# exploring casual model environament spouses and GxE
#######################################################


library(data.table)
library(brms)
library(texreg)
library(foreach)
library(doParallel)
source("src/utils.R")

path = "models/BMI-spouses-GxE/output/"
clusters = 5

# remove files
# files = paste0(path, list.files(path, pattern = "*.csv"))
# sapply(files, unlink, recursive = TRUE)

files = paste0(path, list.files(path, pattern = "bmi"))
dat = rbindlist(lapply(files, fread))

files = paste0(path, list.files(path, pattern = "parameter"))
params = rbindlist(lapply(files, fread))
params = unique(params, by = c("iteration"))
params

# dat[, couple := paste0(sort(.SD), collapse = ""), 
    # .SDcols = c("id", "spouse_id"), 
    # by = .(1:nrow(dat))]
# setorder(dat, couple)

table(dat$iteration) 
table(dat$replicate)

cor(dat$bmi, dat$bmi_spouse)
cor(dat$pgs, dat$pgs_spouse)
summary(dat$susceptability)

estimateModel = function(dat, iter = NULL, clusters = 4) {

    cl = makeCluster(clusters)
    registerDoParallel(cl)

    output = foreach(i = 1:max(params$replicate),
        .packages = c("brms", "cmdstanr", "doParallel", 
            "data.table")) %dopar% {
    
    temp = copy(dat[iteration == iter & replicate == i])
    brm(bmi ~ pgs + bmi_spouse + pgs * bmi_spouse, 
        data = temp, iter = 2000, cores = 1, chains = 1, backend = "cmdstan")
    }
    stopCluster(cl)
    baseline = combine_models(mlist = output, check_data = FALSE)
    return(baseline)
}

m1 = estimateModel(dat, iter = 1)
m2 = estimateModel(dat, iter = 2)
m3 = estimateModel(dat, iter = 3)
m4 = estimateModel(dat, iter = 4)

params

screenreg(list(m1, m2, m3, m4),
    custom.model.names = c("none", "total", "random", "alternate"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)


dat[iteration == 1]
dat[iteration == 2]
dat[iteration == 3]
