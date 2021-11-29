library(brlm)
library(MASS)
library(expectreg)
# install.packages("expectreg")
library(coda)




z <- seq(-2, 2, length.out = 1000)
plot(z, brlm:::fn.psi.quantile(z, tau = .9, deriv = 0), type = "l", ylim = c(-4,4))
lines(z, z*psi.bisquare(z, deriv = 0), type = "l", col = 2)
lines(z, z*psi.huber(z), type = "l", col = 4)




plot(z, brlm:::psi.quantile(z, tau = .9, deriv = 0), type = "l", ylim = c(-4,4))
lines(z, psi.bisquare(z, deriv = 0), type = "l", col = 2)
lines(z, psi.huber(z), type = "l", col = 4)

n <- 10000
x <- rnorm(n, mean = 0, sd = 2)
hist(x, breaks = 100)
rlm(rep(1, n), x, psi =  brlm:::psi.quantile, tau = .9,  maxit = 10000)
enorm(.9, 0, 2)

brlm:::psi.quantile


#####
n <- 1000
#plot(x, mean)
shape = 5
scale = 1
population <- rgamma(10000, shape, scale = scale)
median_est = median(population)
mu0 = mean(population) - median_est
sigma2Int = var(population)


y <- rgamma(n, shape = shape, scale = scale) - median_est
hist(y, breaks = 50)
tau = 0.99
# there is a bug in egamma, must specify rate parameter
(true_expectile = egamma(asy = tau, shape = shape, rate = 1/scale) - median_est)

quant_fit1 <- rlm(y ~ 1 , psi = brlm:::psi.quantile, tau = tau, maxit = 2000)
quant_fit1$coefficients


X <- matrix(rep(1, n), ncol = 1)
Sigma0 = t(X)%*%X*n
a0 = 1
b0 = 1


restricted_fit <- restrictedBayesLm(y, X,
                                    regEst = 'quantile',
                                    scaleEst = 'Huber',
                                    mu0 = mu0,
                                    Sigma0 = Sigma0,
                                    a0 = a0,
                                    b0 = b0,
                                    sigma2Int = sigma2Int,
                                    nkeep = 1e3,
                                    nburn = 1e3,
                                    maxit = 400)


apply(restricted_fit$mcmc,2, mean)
plot(restricted_fit$mcmc)

normal_fit <- bayesLm(y, X, 
                      mu0, 
                      Sigma0, 
                      a0 = a0, 
                      b0 = b0,
                      sigma2Int = sigma2Int, 
                      nkeep=1e4, 
                      nburn=1e3)

apply(normal_fit$mcmc,2, mean)
plot(normal_fit$mcmc)



compute_normal_expectile <- function(mcmc){
  apply(mcmc, 1,  function(x){
    enorm(asy = tau, m = x[1], sd = sqrt(x[2]))
  })
}

normal_fit_expectiles = compute_normal_expectile(normal_fit$mcmc)
resticted_fit_expectiles = compute_normal_expectile(restricted_fit$mcmc)

hist(normal_fit_expectiles)
abline(v = true_expectile, lwd = 2)
hist(resticted_fit_expectiles)
abline(v = true_expectile, lwd = 2)
