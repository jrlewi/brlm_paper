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
  real<lower = 3> nu;
}

transformed parameters{
  real<lower=0> sigma;
  sigma =  pow(sigma2, 0.5);
}

model {
  nu ~ normal(3, 3);
  sigma2 ~  inv_gamma(a, (nu-2)*b/nu);
  beta ~ normal(eta, tau);
  y ~ student_t(nu, beta, sigma);
}
