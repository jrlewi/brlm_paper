# We simulate data from a skewed (shifted) gamma, interest is in an expectile
# two models both with base normal
# 1. typical normal theory Bayes model
# 2. same, but restricted - condition on the expectile of interest

library(brlm)
library(MASS)
library(expectreg)
library(coda)
library(tidyverse)

# similation parameters
n_sims <- 200
set.seed(321)
n <- 1000
shape <- 9
scale <- 5
rate <- 1/scale
tau <- 0.99 # expectile 

# we will shift the gamma so it has negative values.
population <- rgamma(100000, shape, scale = scale)
median_est <- median(population)
# note there is a bug in egamma, must specify rate parameter
(true_expectile = egamma(asy = tau, shape = shape, rate = 1/scale) - median_est)

# prior parameters
mu0 = mean(population) - median_est
sigma2Int = var(population)
X <- matrix(rep(1, n), ncol = 1)
Sigma0 <- t(X)%*%X*n
a0 <- 1
b0 <- 1


## Functions
compute_normal_expectile <- function(mcmc){
  # compute posterior of expectile assuming
  # normal distribution.
  # mcmc is a 2 column matrix: with samples of mean and sigma2 from posterior.

  apply(mcmc, 1,  function(x){
    enorm(asy = tau, m = x[1], sd = sqrt(x[2]))
  })
}

quantile_based_ci <- function(sample, level = 0.95){
  alpha <- (1-level)/2
  quantile(sample, probs = c(alpha, 1 - alpha))
}

expectile_ci <- function(mcmc, level = 0.95){
  sample <- compute_normal_expectile(mcmc)
  quantile_based_ci(sample, level = level)
}

gen_data <- function(n, shape, scale, median_est){
  rgamma(n, shape = shape, scale = scale) - median_est
}



one_fit <- function(){

  y <- gen_data(n, shape, scale, median_est)

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
                                      maxit = 400,
                                      tau = tau)

  normal_fit <- bayesLm(y, X,
                        mu0,
                        Sigma0,
                        a0 = a0,
                        b0 = b0,
                        sigma2Int = sigma2Int,
                        nkeep=1e3,
                        nburn=1e3)
  list(restricted_fit = restricted_fit,
       normal_fit = normal_fit, y = y)
}

one_sim <- function(dummy_argument){
  # run one simulation - extract the ci's, estimates for each model
  one_example <- one_fit()
  restricted_posterior <- compute_normal_expectile(one_example$restricted_fit$mcmc)
  restricted_estimate <- mean(restricted_posterior)
  normal_posterior <- compute_normal_expectile(one_example$normal_fit$mcmc)
  normal_estimate <- mean(normal_posterior)

  restricted_ci <- as.numeric(expectile_ci(one_example$restricted_fit$mcmc))
  normal_ci <- as.numeric(expectile_ci(one_example$normal_fit$mcmc))
  rbind(
    data.frame(posterior = "restricted",
               estimate = restricted_estimate,
               lower = restricted_ci[1],
               upper = restricted_ci[2]),
    data.frame(posterior = "normal",
               estimate = normal_estimate,
               lower = normal_ci[1],
               upper = normal_ci[2])
  )
}

simulation <- function(n_sims){
  results <- lapply(1:n_sims, FUN = one_sim)
  bind_rows(results, .id = "simulation") %>%
    arrange(posterior, simulation)
}



#####
# One example
####

one_example <- one_fit()
y <- one_example$y
y <- gen_data(n, shape, scale, median_est)

SAVE_PATH = "/Users/johnlewis/personal/brlm_paper/revision_2/code/expectile_example"
ggplot(data.frame(y = y), aes(x = y)) + geom_histogram(fill = "white", col = 1, bins = 50)
ggsave(file.path(SAVE_PATH, paste0("data_hist", "_", n, "_", tau, ".png")))

true_expectile
expectile_ci(one_example$restricted_fit$mcmc)
expectile_ci(one_example$normal_fit$mcmc)

#expectile posteriors
restricted_posterior <- compute_normal_expectile(one_example$restricted_fit$mcmc)
normal_posterior <- compute_normal_expectile(one_example$normal_fit$mcmc)

post_df <- rbind(data.frame(posterior = "restricted_posterior", expectile = restricted_posterior),
                 data.frame(posterior = "normal_posterior", expectile = normal_posterior))


ggplot(post_df, aes(x=expectile, group = posterior, col = posterior)) +
  geom_density(lwd = 1.1) + scale_color_brewer(palette ="Set2") +
  geom_vline(xintercept = true_expectile, lty = 2) +
  ggtitle("Posterior of expectile under the two models:
          vertical dashed line is the true expectile")
ggsave(file.path(SAVE_PATH, paste0("expectile_posterior_example", "_", n, "_", tau, ".png")))


# simulation
(strt <- Sys.time())
results <- simulation(n_sims = n_sims)
(end <- Sys.time() - strt)

saveRDS(results, file = file.path(SAVE_PATH,
                                  paste0("simulation_results", "_", n, "_", tau, ".rds")))
saveRDS(true_expectile, file = file.path(SAVE_PATH,
                                         paste0("true_expectile", "_", n, "_", tau, ".rds")))

# compute MSE
summary_stats <- function(results){
  results %>%
    group_by(posterior) %>%
    summarise(mse = mean((estimate - true_expectile)^2),
              bias = mean(estimate - true_expectile),
              coverage = mean((lower <= true_expectile) &
                                (true_expectile <= upper)))
}

(sum_stats <- summary_stats(results))
saveRDS(sum_stats, file = file.path(SAVE_PATH,
                                    paste0("summary_stats", "_", n, "_", tau, ".rds")))


#detecting deviation from normality

n = 100
shape = 9
scale = 5


perform_norm_test <- function(aa){
  x <- rgamma(n, shape = shape, scale = scale)
  shap_test <- shapiro.test(x)
  shap_test$p.value
}


p_vals <- sapply(1:10000, perform_norm_test)
mean(p_vals > 0.05)

