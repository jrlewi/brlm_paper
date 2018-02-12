data {
 int<lower = 0> N;
 vector[N] y;
 real eta;
 real tau;
 real a;
 real b;
}

parameters {
  real beta;
  real sigma2;
  //real<lower = 1> c2;
  real<lower=0, upper=1> p;
}

transformed parameters{
  real<lower=0> sigma;
  //real c;
    sigma =  pow(sigma2, 0.5);
   // c = pow(c2, 0.5);
}

model {
 sigma2 ~  inv_gamma(a, b);
 beta ~ normal(eta, tau);
 p ~ beta(20, 1);
 //c2 ~ normal(5, 5);
 for (n in 1:N)
   target += log_mix(p,
                     normal_lpdf(y[n] | beta, sigma),
                     normal_lpdf(y[n] | beta, pow(10, 0.5)*sigma));
}
