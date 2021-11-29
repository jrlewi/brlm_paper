library(brlm)
library(MASS)





z <- seq(-2, 2, length.out = 1000)
plot(z, brlm:::psi.quantile(z, tau = .1), type = "l")
plot(z, psi.bisquare(z), type = "l")
plot(z, psi.huber(z), type = "l")


n <- 10000
x <- rnorm(n, mean = 10, sd = 2)

rlm(rep(1, n), x, psi = brlm:::psi.quantile, tau = .1,  maxit = 100)
qnorm(.1, 10, 2)





#####
n <- 10000
x <- runif(n)
beta = 20
mean <- beta*x #cos(8*pi*x) + 3*sin(8*pi*x)

#plot(x, mean)
shape = 5
scale = 2
median_est = median(rgamma(10000, shape, scale = scale))
error <- rgamma(n, shape = shape, scale = scale) - median_est

quantile_est = quantile(rgamma(10000, shape, scale = scale), probs = .9) - median_est

y <- mean + error
plot(x, y, cex = .1)
points(x, mean, cex  = .5, col = 2)
points(x, mean + quantile_est, col = 4, cex = .5)

#cos(8*pi*x) + sin(8*pi*x)
quant_fit1 <- rlm(y ~ x , psi = psi.quantile, tau = .9, maxit = 2000)

lines(sort(x), predict(quant_fit1, data.frame(x = sort(x))), col = "green")


normal_fit1 <- lm(y ~ x)
lines(sort(x), predict(normal_fit1, 
                       data.frame(x = sort(x))) + qnorm(.9)*sigma(normal_fit1),
      col = "yellow")
      
lines(sort(x), predict(normal_fit1, 
                       data.frame(x = sort(x))),
      col = "yellow")



X <- matrix(x, ncol = 1)

restricted_fit <-
  restrictedBayesLm(y, X,
                    regEst = "quantile",
                    mu0 = beta,
                    Sigma0 = n*(t(X)%*%X)^-1,
                    a0 = 1,
                    b0 = 1,
                    sigma2Int = 1,
                    nkeep = 10,
                    nburn = 1)
      