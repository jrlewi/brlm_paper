# Fitting several models to speed of light data

library(brlm)
library(MASS)
library(mvtnorm)
library(coda)
library(MCMCpack)

data(newcomb)
y <- newcomb

# prior parameters ----
eta <- 23.6; tau <- 2.04; alpha <- 5; beta <- 10;


# fitting parameters for rl_direct ----
mu_lims <- c(20, 35); length_mu <- 300;
sigma2_lims <- c(1, 60); length_sigma2 <- 500;
N <- 5000


# Tukey Fit ---- 
rlm_tukey <- function(y){
  fit <- rlm(y~1,psi = psi.bisquare,  scale.est = 'proposal 2', maxit=50)
  c(fit$coefficients, fit$s)
}

set.seed(1)
fit_tukey <- brlm:::rl_direct(y, statistic = rlm_tukey, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)

saveRDS(fit_tukey, 'out/fit_tukey.rds')

# Huber Fit ----
rlm_huber <- function(y){
  fit <- rlm(y~1,psi = psi.huber,  scale.est = 'proposal 2', maxit=50)
  c(fit$coefficients, fit$s)
}


set.seed(1)
fit_huber <- brlm:::rl_direct(y, statistic = rlm_huber, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)
saveRDS(fit_huber, 'out/fit_huber.rds')

# LMS fit -----
rlm_lms <- function(y){
  fit <- lqs(y~1, method='lms')
  c(fit$coefficients, fit$scale[1])
}


set.seed(1)
fit_lms <- brlm:::rl_direct(y, statistic = rlm_lms, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)
saveRDS(fit_lms, 'out/fit_lms.rds')


# LTS fit -----
rlm_lts <- function(y){
  n <- length(y)
  quantile<-floor(.95*n)
  fit <- lqs(y~1, method='lts', quantile = quantile)
  c(fit$coefficients, fit$scale[1])
}

set.seed(1)
fit_lts <- brlm:::rl_direct(y, statistic = rlm_lts, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)
saveRDS(fit_lts, 'out/fit_lts.rds')


# fit normal model ----
set.seed(1)
X <- matrix(rep(1, length(y)))
fit_normal <- brlm::bayesLm(y, X, mu0 = eta, Sigma0 = tau^2, a0 = alpha, b0 = beta, sigma2Int = mad(y)^2, nkeep = 10000, nburn = 5000)
traceplot(fit_normal$mcmc)

saveRDS(fit_normal$mcmc, 'out/fit_normal.rds')


# fit t-model -----

set.seed(1)
nu <- 5 # degrees of freedom for the t
beta_t <- beta*(nu-2)/nu #want the prior on the var(Y)=(nu/(n-2))sigma2 ~IG(alpha, beta) so sigma2=((n-2)/nu)Var(Y)~IG(alpha,((n-2)/nu) beta)

fit_t <- brlm::bayesTdistLm(y, X, mu0 = eta, Sigma0 = tau^2, a0 = alpha, b0 = beta_t, nu = nu, nkeep = 20000, nburn = 10000, rwTune = 6)

traceplot(fit_t$mcmc)
mean(fit_t$acceptSigma2)

saveRDS(fit_t$mcmc, 'out/fit_t.rds')



rlm_tukey(y)
rlm_huber(y)
rlm_lts(y)
rlm_lms(y)
mean(y)
