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
  ordered[2] beta0;
  ordered[2] beta1;
  vector[2] sigma2;
  // real<lower=0> sigma2_2;
  real<lower=0, upper=1> p;
}

transformed parameters{
  vector[2] sigma;
  // real<lower=0> sigma_2;
  for(i in 1:2)
    sigma[i] =  pow(sigma2[i], 0.5);
}

model {
 sigma2 ~  inv_gamma(a0, b0);
 beta0 ~ normal(mu0[1], pow(Sigma0[1,1], .5));
 beta1 ~ normal(mu0[2], pow(Sigma0[2,2], .5));
 p ~ beta(5, 1);
 for (n in 1:N)
   target += log_mix(p,
                     normal_lpdf(y[n] |beta0[1] + beta1[1]*x[n], sigma[1]),
                     normal_lpdf(y[n] | beta0[2] + beta1[2]*x[n], sigma[2]));
}
