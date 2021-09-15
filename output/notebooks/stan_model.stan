 
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real alpha_sigma;
  real beta_sigma;
}
model {
    vector[N] sigma = x * beta_sigma  + alpha_sigma;
    for (n in 1:N) {
         sigma[n] = (sqrt(sigma[n]))^2;
    }
    y ~ normal(x * beta + alpha, sigma);
}


