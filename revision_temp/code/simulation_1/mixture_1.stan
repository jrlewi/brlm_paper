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
  real sigma2;
  real<lower = 1> c;
  real<lower=0, upper=1> p;
}

transformed parameters{
  real<lower=0> sigma;
  sigma =  pow(sigma2, 0.5);
}

model {
 sigma2 ~  inv_gamma(a, b);
 beta ~ normal(eta, tau);
 p ~ beta(9, 1);
 c ~ normal(1, 4);
 for (i in 1:n)
   target += log_mix(p,
                     normal_lpdf(y[i] | beta, sigma),
                     normal_lpdf(y[i] | beta, c*sigma));
}
