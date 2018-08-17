data {
 int<lower = 0> n;
 vector[n] y;
 real eta;
 real tau;
 real a;
 real b;
}

parameters {
  real beta;
  real<lower = 0> sigma2;
}

transformed parameters{
  real<lower=0> sigma;
    sigma =  pow(sigma2, 0.5);
}

model {
 sigma2 ~  inv_gamma(a, b);
 beta ~ normal(eta, tau);
 y ~ normal(beta, sigma);  
}
