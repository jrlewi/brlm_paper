# Analyze the Belgium Calls Data

library(tidyverse)
library(MASS)
library(rstan)

log_scale  <- 1

phones$year<- phones$year - mean(phones$year)
phones <- bind_cols(phones)
if(log_scale) phones$calls <- log(phones$calls)
plot(phones$year, phones$calls)

# mn_sd <- summarise_all(phones, funs(mean, sd))
# phones <- as_tibble(scale(phones))


# Restricted fit with Tukey locataion and MAD for scale-----

n_prior <- 3
(fit <- lm(calls ~ year, data = phones[1:n_prior,]))


(lm_fit <- lm(calls ~ year, data = phones))
# prior parameters -----
g <- nrow(phones) - n_prior
mu0 <- fit$coefficients
X <- model.matrix(fit)
Sigma0 <- g*sigma(fit)^2*solve(t(X)%*%X)
a0 <- 2
b0 <- 1 

mu0 + 3*sqrt(diag(Sigma0))
mu0 - 3*sqrt(diag(Sigma0))
mu0 + 3*(diag(Sigma0))
mu0 - 3*(diag(Sigma0))


y <- phones$calls[-c(1:n_prior)]
X_fit <- cbind(1, phones$year[-c(1:n_prior)])
x <- X_fit[,2]
N <- length(y)
stan_rdump(c("N", "y", 'x', 'mu0', 'Sigma0', 'a0', 'b0'), file="mix.data.R")
library(MASS)
summary(rlm(y~x,  psi = 'psi.bisquare'))
summary(lm(y~x))
input_data <- read_rdump("mix.data.R")

fit_mix <- stan(file='calls_mixture.stan', data=input_data,
            chains=4, seed=123, iter = 1e4, refresh=2000)

saveRDS(fit_mix, 'out/mixture_stan.rds')
fit_mix <- readRDS('out/mixture_stan.rds')
rstan::get_posterior_mean(fit_mix)

plot(fit_mix)
stan_trace(fit_mix, window = c(9000,10000))
# stan_scat(fit_mix, pars = c("beta0[1]", "beta0[2]"))
# stan_scat(fit_mix, pars = c("beta1[1]", "beta1[2]"))
# stan_scat(fit_mix, pars = c("sigma[1]", "sigma[2]"))

# stan_scat(fit_mix, pars = c("sigma", "beta[2]"))
# stan_scat(fit_mix, pars = c("p", "sigma"))
stan_dens(fit_mix, fill = 'white')
rlmf <- rlm(y~x, psi = psi.bisquare)
mcmc_samps <- rstan::extract(fit_mix)
plot(x, y)
abline(mean(mcmc_samps$beta0[,1]),mean(mcmc_samps$beta1[,1]))
abline(coef(rlmf)[1],coef(rlmf)[2], col = 2)


saveRDS(mcmc_samps, 'out/mixture_fit.rds')
