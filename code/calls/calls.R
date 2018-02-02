# Analyze the Belgium Calls Data

library(tidyverse)
library(brlm)
library(MASS)
library(mvtnorm)
library(coda)
library(MCMCpack)

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


y <- phones$calls[-c(1:n_prior)]
X_fit <- cbind(1, phones$year[-c(1:n_prior)])

# Importance Sampling 
# define statistic
rlm_tukey <- function(y, X){
  fit <- rlm(X, y, psi = psi.bisquare,  scale.est = 'MAD', maxit=1e4)
  c(fit$coefficients, fit$s)
}
# importance distribution parameters 
rlm_tukey(y, X_fit)

rlm_fit <- rlm(calls ~ year, data = phones[-c(1:n_prior),], psi = psi.bisquare, scale.est = 'MAD')
c(rlm_fit$coefficients, rlm_fit$s)
plot(phones$year, phones$calls)
lines(phones$year[-c(1:n_prior)], fitted(rlm_fit))


i <- 2
par(mfrow = c(2,1))
smp1 <- rnorm(1000, mu0[i], sqrt(Sigma0[i,i]))
smp2 <- rnorm(1000, rlm_fit$coefficients[i], sqrt(vcov(rlm_fit)[i,i]))
hist(smp1, xlim = range(smp1, smp2))
hist(smp2, xlim = range(smp1, smp2))
par(mfrow = c(1,1))



cov_b <- vcov(rlm_fit)
scale <- .25*rlm_fit$s

N <- 1e5
Nins <- 1e5
set.seed(125)
rl_fit <- rl_importance(y, X_fit, statistic = rlm_tukey, mu0 = mu0, Sigma0 = Sigma0, alpha = a0, beta = b0, cov_b = cov_b, scale = scale, smooth = 1, N = N,Nins = Nins)

# check weights -----
# plot(log(rl_fit$w + 1e-20), cex = .2)
plot(rl_fit$w, cex = .2)


write_rds(rl_fit, path = file.path(here::here(), 'rl_fit.rds'))



#fit t - distribution ----
nkeep <- 1e5
nburn <- 3e4
nu <- 5
set.seed(126)
t_fit <-bayesTdistLm(y=y, X_fit, mu0, Sigma0, a0=a0, b0=((nu-2)/nu)*b0, parInit=c(rlm_fit$coefficients, (nu-2)*sigma(fit)^2/nu), nu=nu, nkeep=nkeep, nburn=nburn, rwTune=sigma(lm_fit)*5)

coda::traceplot(t_fit$mcmc)

write_rds(t_fit, path = file.path(here::here(), 't_fit.rds'))

# fit normal model on subset of the data.
nkeep <- 1e4
nburn <- 2e3
normal_ind <- c(4:13,21:24)
y_normal <- phones$calls[normal_ind]
X_normal <- cbind(1, phones$year[normal_ind])
set.seed(127)
normal_fit <- bayesLm(y_normal, X_normal, mu0, Sigma0, a0, b0, sigma2Int = rlm_fit$s^2, nkeep = nkeep, nburn = nburn)

coda::traceplot(normal_fit$mcmc,smooth = TRUE)
write_rds(normal_fit, path = file.path(here::here(), 'normal_fit.rds'))
