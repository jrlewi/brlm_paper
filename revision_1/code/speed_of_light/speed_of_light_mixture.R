library(MASS)
library(rstan)
library(bayesplot)

data(newcomb)
y <- newcomb
N <- length(y)
# prior parameters ----
eta <- 23.6; tau <- 2.04; a <- 5; b <- 10;


stan_rdump(c("N", "y", 'eta', 'tau', 'a', 'b'), file="mix.data.R")

input_data <- read_rdump("mix.data.R")

fit <- stan(file='speed_of_light_mixture.stan', data=input_data,
                 chains=4, seed=123, iter = 1e4, refresh=2000)
fit
plot(fit)
mcmc_intervals(fit)
stan_trace(fit)
stan_scat(fit, pars = c("beta", 'p'))
stan_scat(fit, pars = c("sigma", "beta"))
stan_scat(fit, pars = c("p", "sigma"))
stan_dens(fit)

mcmc_samps <- rstan::extract(fit)
plot(density(mcmc_samps$p))
lines(density(rbeta(1000, 20, 1)), col = 4)
abline(v = 64/N)

pr_p <- rbeta(10000, 20, 1)
boxplot(pr_p, mcmc_samps$p)
quantile(mcmc_samps$p, c(.025, .975))
quantile(pr_p, c(.025, .975))
mean(mcmc_samps$p)
mean(pr_p)
median(mcmc_samps$p)
64/N
sd(mcmc_samps$p)/sd(pr_p)



y_tilde <- rnorm(length(mcmc_samps$beta), mcmc_samps$beta, mcmc_samps$sigma)
hist(y_tilde, xlim = range(y, y_tilde))
points(y, rep(1, N), pch = 19, col = 2)

#predictive checks----
y_tilde <- t(sapply(1:10000, function(a) {
  ind <- sample(1: length(mcmc_samps$beta),1)
  x <- rnorm(N, mcmc_samps$beta[ind], mcmc_samps$sigma[ind])
  c(mean(x), sd(x))
}))
plot(y_tilde, pch = 19, cex = .5)
points(median(y), mad(y), col = 2, pch = 19)
points(mean(sort(y)[-c(1:2)]), sd(sort(y)[-c(1:2)]), col = 3, pch = 19)

quantile(mcmc_samps$beta, c(.05, .95))
quantile(y_tilde[,1], c(.05, .95))
plot(density(mcmc_samps$beta))
lines(density(y_tilde[,1]))
quantile(mcmc_samps$sigma, c(.05, .95))
quantile(y_tilde[,2], c(.05, .95))
plot(density(mcmc_samps$sigma))
lines(density(y_tilde[,2]))


plot(density(mcmc_samps$p))
lines(density(rbeta(1000, 20, 1)))


saveRDS(mcmc_samps, 'out/fit_mixture.rds')
