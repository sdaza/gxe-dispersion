// generated with brms 2.15.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_vest;  // number of population-level effects
  matrix[N, K_vest] X_vest;  // population-level design matrix
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[K_vest] b_vest;  // population-level effects
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_vest = X_vest * b_vest;
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize non-linear predictor term
    vector[N] sigma;
    for (n in 1:N) {
      // compute non-linear predictor values
      sigma[n] = (sqrt(nlp_vest[n])) ^ 2;
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += student_t_lpdf(b_vest[2] | 7, 0, 8);
  target += student_t_lpdf(Intercept | 3, -0.1, 2.5);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

