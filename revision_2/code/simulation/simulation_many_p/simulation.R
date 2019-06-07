# Consider regression problem with many p, larger n,
# and outliers

library(tidyverse)
library(brlm)
library(MCMCpack)

data_generation <- function(n, p,
                               a_0,
                               b_0,
                               prob_outlier,
                               outlier_contam){
  X <- matrix(runif(n*p), n, p)
  beta <- rnorm(p, 0, 1)
  sigma2 <- rinvgamma(1, a_0, b_0)
  expected <- X%*%beta
  y <- sapply(expected, function(mn){
    outlier <- rbinom(1,1, prob_outlier)
    vr_out <- (outlier_contam)*sigma2
    err <- if(outlier){
      abs(rnorm(1, 0, sd = sqrt(vr_out)))
    } else {
      rnorm(1, 0, sd = sqrt(sigma2))
    }
    mn + err
  })
  data <- cbind(y, X)
  params <- c(beta, sigma2)
  list(data = data, params = params)
}

data <- data_generation(1000, 3, 
                5,
                5, 
                .1,
                10)

pairs(data$data)
data$params
