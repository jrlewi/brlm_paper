data {
 int<lower = 0> N;
 vector[N] y;
 vector[N] x;
 vector[2] mu0;
 matrix[2,2] Sigma0;
 real a0;
 real b0;
}

parameters {
  vector[2] beta;
  real sigma2;
}

transformed parameters{
  real sigma;
  // real<lower=0> sigma_2;
    sigma =  pow(sigma2, 0.5);
}

model {
 sigma2 ~  inv_gamma(a0, 0.6*b0);
 beta ~  multi_normal(mu0,   Sigma0);
 for (n in 1:N)
  y[n] ~ student_t(5, beta[1] + x[n]*beta[2], sigma);
}
