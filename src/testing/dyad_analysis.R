
library(data.table)
library(brms)


# read data
dat = fread("data/dyads.csv")

# format dataset
ch8_f = copy(dat[female == 1])
colnames(ch8_f) <- c("coupleid", "f_personid", "time", 
                     "time7c", "gender", "female", "male", "f_reldis",
                     "f_wrkstrs", "f_wrkstrsc", "f_wrkstrscb", "f_wrkstrscw")
head(ch8_f)


ch8_m = copy(dat[male == 1, .SD, .SDcols = !names(dat) %in% c("time7c", "gender", "male", "female")]) 
colnames(ch8_m) <- c("coupleid", "m_personid", "time", 
                     "m_reldis", "m_wrkstrs", "m_wrkstrsc", "m_wrkstrscb", "m_wrkstrscw")
head(ch8_m)

ch8 = merge(ch8_f, ch8_m, by = c("coupleid", "time"))
head(ch8)

dim(ch8)
# model
f_model = bf(f_reldis ~ 0 + intercept + time7c + f_wrkstrscw + f_wrkstrscb + (time7c + f_wrkstrscw |p| coupleid))
m_model = bf(m_reldis ~ 0 + intercept + time7c + m_wrkstrscw + + m_wrkstrscb + (time7c + m_wrkstrscw |p| coupleid))

fit <- brm(f_model + m_model, data = ch8, 
           seed = 111, # set a seed to ensure reproducibility
           iter = 8000, # number of iterations
           autocor = cor_ar(~ time7c | coupleid, p = 1), 
           backend = "cmdstan")


library(rethinking)
data(KosterLeckie)

kl_data = 
  list(
    N            = nrow(kl_dyads),
    N_households = max(kl_dyads$hidB), 
    did          = kl_dyads$did,
    hidA         = kl_dyads$hidA,
    hidB         = kl_dyads$hidB,
    giftsAB      = kl_dyads$giftsAB, 
    giftsBA      = kl_dyads$giftsBA
  )

kl_data

m14.7 =
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

precis(m14.7, 2)
