#######################################################
# exploring casual model environament spouses and GxE
#######################################################


library(data.table)
library(brms)
library(texreg)
library(foreach)
library(doParallel)
library(rethinking)
source("src/utils.R")

path = "models/BMI-spouses-GxE/output/"
clusters = 2

# # remove files
# files = paste0(path, list.files(path, pattern = "*.csv"))
# sapply(files, unlink, recursive = TRUE)

files = paste0(path, list.files(path, pattern = "bmi"))
dat = rbindlist(lapply(files, fread))

files = paste0(path, list.files(path, pattern = "parameter"))
params = rbindlist(lapply(files, fread))
params = unique(params, by = c("iteration"))
params

dat[, seq := 1:.N, .(couple_id, iteration, replicate)]
table(dat$seq)

# correlations
test = dat[seq == 2]
cor(test$ego_bmi_married, test$alter_bmi_married)
cor(test$ego_bmi_single, test$alter_bmi_single)
cor(test$ego_pgs, test$alter_pgs)
cor(test$ego_environment, test$alter_environment)

# create basic model
temp = copy(dat[replicate == 1 & iteration == 1])
couple_data =
    list(
        N = nrow(temp),
        couple_id = temp$couple_id, 
        bmiA = temp$ego_bmi_married, 
        bmiB = temp$alter_bmi_married, 
        pgsA = temp$ego_pgs, 
        pgsB = temp$alter_pgs,
        pgsA_bmiB = temp$ego_pgs * temp$alter_bmi_married, 
        pgsB_bmiA = temp$alter_pgs * temp$ego_bmi_married
)

srm = alist(
    bmiA ~ normal(muA, sigmaA),
    bmiB ~ normal(muB, sigmaB),

    muA <- aA + ba_pgsA * pgsA +  ba_bmiB * bmiB + ba_pgsAbmiB * pgsA_bmiB,
    muB <- aB + bb_pgsB * pgsB + bb_bmiA * bmiA + bb_pgsBbmiA * pgsB_bmiA,

    c(ba_bmiB, ba_pgsA, ba_pgsAbmiB)  ~ normal(0, 1),
    c(bb_bmiA, bb_pgsB, bb_pgsBbmiA)  ~ normal(0, 1),
    c(aA, aB) ~ normal(0, 1),
    c(sigmaA, sigmaB) ~ exponential(1)

)


srmC = alist(
    bmiA ~ normal(muA, sigmaA),
    bmiB ~ normal(muB, sigmaB),

    muA <- aA + ba_pgsA * pgsA +  ba_bmiB * bmiB + ba_pgsAbmiB * pgsA_bmiB + d[couple_id, 1],
    muB <- aB + bb_pgsB * pgsB + bb_bmiA * bmiA + bb_pgsBbmiA * pgsB_bmiA + d[couple_id, 2],

    c(ba_bmiB, ba_pgsA, ba_pgsAbmiB)  ~ normal(0, 1),
    c(bb_bmiA, bb_pgsB, bb_pgsBbmiA)  ~ normal(0, 1),
    c(aA, aB) ~ normal(0, 1),
    c(sigmaA, sigmaB) ~ exponential(1),

    ## dyad effects
    transpars> matrix[N,2]:d <-
    compose_noncentered(rep_vector(sigma_d, 2), L_Rho_d, z), 
    matrix[2,N]:z ~ normal(0, 1),
    cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky(8), 
    sigma_d ~ exponential(1),
    
    ## compute correlation matrix for dyads
    gq> matrix[2, 2]:Rho_d <<- Chol_to_Corr(L_Rho_d)

)

m0 = ulam(srm, data = couple_data, chains = 1, cores = 1, 
            iter = 2000, cmdstan = TRUE
)
precis(m0)

m1= ulam(srmC, data = couple_data, chains = 1, cores = 1, 
            iter = 2000, cmdstan = TRUE
)
precis(m1)
precis(m1, depth = 3, pars = c("Rho_d", "sigma_d"))


# no random mating
no_random_mating = runModel(srm, replicates = 20, data = dat[iteration == 1 & seq == 1])
precis(no_random_mating$model, prob = 0.95)
no_random_mating$rhat
no_random_mating$neff

# random mating
random_mating = runModel(srm, replicates = 20, data = dat[iteration == 2 & seq == 1])
precis(random_mating$model, prob = 0.95)
random_mating$rhat
random_mating$neff





# other models

m14.7 <- 
  ulam( 
    alist(
      giftsAB ~ poisson(lambdaAB),
      giftsBA ~ poisson(lambdaBA),
      log(lambdaAB) <- a + gr[hidA, 1] + gr[hidB, 2] + d[did, 1] , 
      log(lambdaBA) <- a + gr[hidB, 1] + gr[hidA, 2] + d[did, 2] , 
      a ~ normal(0, 1),
      
      ## gr matrix of varying effects
      vector[2]:gr[N_households] ~ multi_normal(0, Rho_gr, sigma_gr), 
      Rho_gr ~ lkj_corr(4),
      sigma_gr ~ exponential(1),
      
      ## dyad effects
      transpars> matrix[N,2]:d <-
        compose_noncentered(rep_vector(sigma_d, 2), L_Rho_d, z), 
      matrix[2,N]:z ~ normal(0, 1),
      cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky(8), 
      sigma_d ~ exponential(1),
      
      ## compute correlation matrix for dyads
      gq> matrix[2, 2]:Rho_d <<- Chol_to_Corr(L_Rho_d)
    ), 
    data = kl_data, 
    chains = 4, cores = 4, iter = 2000
  )

f1 = formula(ego_bmi_married ~ ego_pgs + alter_bmi_married + ego_pgs * alter_bmi_married +
    (1|couple_id))




table(dat$iteration)
m1.1 = estimateModel(f1, dat, iter = 1)
summary(m1.1)


s = dat[seq == 1]
ego_model = bf(ego_bmi_married ~ 1+  ego_pgs * alter_bmi_married)
alter_model = bf(alter_bmi_married ~ 1 + alter_pgs * ego_bmi_married)

fit = brm(ego_model + alter_model, data = s[replicate == 1], 
           seed = 111, # set a seed to ensure reproducibility
           iter = 2000, # number of iterations
           chain = 1, 
           backend = "cmdstan"
        )

summary(fit)

dat[replicate == 1 & couple_id == 3]

# testing environment
f1 = formula(bmi_single ~ pgs + bmi_single_spouse + pgs * bmi_single_spouse +
    (1|couple_id))
f2 = formula(bmi_married ~ pgs + bmi_married_spouse + pgs * bmi_married_spouse +
    (1|couple_id))
f3 = formula(bmi_married ~ pgs + environment +  pgs * environment +
     household + pgs * household + (1|couple_id))

m1.1 = estimateModel(f1, dat, iter = 1)
m1.2 = estimateModel(f2, dat, iter = 1)
m1.3 = estimateModel(f3, dat, iter = 1)
screenreg(list(m1.1, m1.2, m1.3),
    custom.model.names = c("single BMI", "married BMI", "complete model"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)

# testing susceptibility 
f1 = formula(bmi_single ~ pgs + bmi_single_spouse + pgs * bmi_single_spouse)
f2 = formula(bmi_married ~ pgs + bmi_married_spouse + pgs * bmi_married_spouse)

m4.1 = estimateModel(f1, dat[susceptibility < 0.5], iter = 6)
m4.2 = estimateModel(f2, dat[susceptibility < 0.5], iter = 6)
m4.3 = estimateModel(f1, dat[susceptibility > 0.5], iter = 6)
m4.4 = estimateModel(f2, dat[susceptibility > 0.5], iter = 6)
screenreg(list(m4.1, m4.2, m4.3, m4.4),
    custom.model.names = c("S=0: single BMI", "S=0: married BMI", 
        "S=1: single BMI", "S=1: married BMI"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)

# testing betas
dat[, beta5 := (bmi_married - (1 + susceptibility * 0.01) * bmi_single_spouse) /
     (bmi_single_spouse * susceptibility * pgs_spouse)
]
mean(dat[iteration == 3 & susceptibility < 0.5]$beta5)


# correlations
cor(dat$bmi_married, dat$bmi_married_spouse)
cor(dat$bmi_single, dat$bmi_single_spouse)
cor(dat$pgs, dat$pgs_spouse)
cor(dat$environment, dat$environment_spouse)

plot(dat[, .(environment_spouse, environment)])
plot(dat[, .(bmi_single, bmi_single_spouse)])


test = dat$bmi_single
abs(sample(test, 1)  - sample(test, 1))

summary(dat$bmi_single)
summary(dat$environment)



table(dat$iteration) 
table(dat$replicate)

summary(dat$bmi_single)
summary(dat$bmi_single_spouse)

cor(dat$bmi_married, dat$bmi_married_spouse)
cor(dat$bmi_single, dat$bmi_single_spouse)

dat
cor(dat$environment, dat$environment_spouse)
cor(dat$pgs, dat$pgs_spouse)
summary(dat$susceptibility)
summary(dat$household)

cor(dat$environment, dat$household)
table(dat[, .N, household]$N)

params

# testing suscetibility 
testing = dat[iteration == 3]
s1 = testing[susceptibility == 0]
s2 = testing[susceptibility == 1]

cor(s1[, .(bmi_single, bmi_single_spouse)])
cor(s1[, .(bmi_married, bmi_married_spouse)])
cor(s2[, .(bmi_single, bmi_single_spouse)])
cor(s2[, .(bmi_married, bmi_married_spouse)])

# summary(lm(bmi_single ~ pgs * bmi_single_spouse, data = s1))
# summary(lm(bmi_single ~ pgs * bmi_single_spouse, data = s2))

# a = lm(bmi_married ~ pgs * bmi_married_spouse, data = s1)
# b = lm(bmi_married ~ pgs * bmi_married_spouse, data = s2)
# screenreg(list(a,b))

table(dat[iteration == 3]$susceptibility)






f3 = formula(bmi_married ~ pgs + household + pgs * household)
f4 = formula(bmi_married ~ pgs + environment + pgs * environment)

test = estimateModel(f3, dat, iter = 1)
summary(test)

test = estimateModel(f4, dat, iter = 1)
summary(test)
params

m1.1 = estimateModel(f1, dat, iter = 1)
m1.2 = estimateModel(f2, dat, iter = 1)
screenreg(list(m1.1, m1.2),
    custom.model.names = c("single", "married"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)

m2.1 = estimateModel(f1, dat, iter = 2)
m2.2 = estimateModel(f2, dat, iter = 2)
screenreg(list(m2.1, m2.2),
    custom.model.names = c("single", "married"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)

m3.1 = estimateModel(f1, dat, iter = 3)
m3.2 = estimateModel(f2, dat, iter = 3)
screenreg(list(m3.1, m3.2),
    custom.model.names = c("single", "married"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)


names(dat)


sd(s1$bmi_married)
sd(s2$bmi_married)
sd(s1$bmi_married_spouse)
sd(s2$bmi_married_spouse)

cor(s2[, .(bmi_single, bmi_single_spouse)])
cor(s2[, .(bmi_married, bmi_married_spouse)])
summary(lm(bmi_single ~ bmi_single_spouse, data = s2))
summary(lm(bmi_married ~ bmi_married_spouse, data = s2))
sd(s2$bmi_married_spouse)



m2 = estimateModel(dat, iter = 2)
m3 = estimateModel(dat, iter = 3)
m4 = estimateModel(dat, iter = 4)

params
table(dat[iteration == 1, susceptibility])

screenreg(list(m1, m2, m3, m4),
    custom.model.names = c("none", "total", "random", "alternate"), 
    include.rsquared = FALSE,
    include.loo.ic = FALSE)


dat[iteration == 1]
dat[iteration == 2]
dat[iteration == 3]
