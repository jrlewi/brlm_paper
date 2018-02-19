# Simulation
  # simulate data from a 'true model'
  # fit various models. 
  # estimate mean
  # repeat



#Note - when re-run, previous results are deleted

library(brlm)
library(rstan)
library(tidyverse)
library(MASS)
library(MCMCpack)

set.seed(123)

# choose true model ----
#either normal, mixture 1, mixture 2, or t
true_model <- 't'

#simulation parameters ----
n_sim <- 50
samp_sizes <- c(25, 50, 100)


# create data_gen function----
#generates true data
if(true_model == 'normal'){
data_gen <- function(n){
  rnorm(n)
}
} else if(true_model == 'mixture 1'){
  data_gen <- function(n){
    x <- rbinom(n, 1, .9)
    sapply(x, function(a) ifelse(a, rnorm(1), rnorm(1, sd = 3)))
  }
} else if(true_model == 'mixture 2'){
  data_gen <- function(n){
  x <- rbinom(n, 1, .9)
  sapply(x, function(a) ifelse(a, rnorm(1), rnorm(1, mean = 3,  sd = 3)))
  }
} else if(true_model == 't'){
  data_gen <- function(n){
  rt(n,df = 3)/sqrt(3)
  }
}




# prior parameters ----
eta <- 0
tau <- 3
a <- 7
b <- 10


# set up empty tibbles to store results----
empty_mat <- matrix(NA, n_sim, length(samp_sizes))
colnames(empty_mat) <- samp_sizes
empty_mat <- as_tibble(empty_mat)

normal_post_means <- mixture_1_post_means <- mixture_2_post_means <- t_post_means <- tukey_post_means <- huber_post_means <- empty_mat 

# initiate log file ----
out <- here::here(paste0('results_','true_model_', true_model))

if(dir.exists(out)){
 unlink(file.path(out, '*')) #removes previous results
} else {
  dir.create(out)
}

writeLines(paste('n  sim' , 'at', Sys.time(), '\n'),file.path(out, 'log.txt'))


# begin for loop ----
for(i in seq_along(samp_sizes)){
n <- samp_sizes[i]
stan_data <- list(n = n, eta = eta, tau = tau, a = a, b = b)
  for(j in seq(n_sim)){

#generare data ----
y <- data_gen(n)
stan_data$y <- y

# fit normal model ----
fit_normal <- stan(file='standard_normal.stan', data=stan_data, chains=1, iter = 4e3, refresh = -1, verbose = FALSE)
normal_post_means[j,i] <- get_posterior_mean(fit_normal, pars = 'beta')
write_rds(normal_post_means,file.path(out, 'normal_post_means.rds'))

# fit mixture 1 ----
fit_mixture_1 <- stan(file='mixture_1.stan', data=stan_data, chains=1, iter = 4e3, refresh = -1, verbose = FALSE)

mixture_1_post_means[j,i] <- get_posterior_mean(fit_mixture_1, pars = 'beta')
write_rds(mixture_1_post_means, file.path(out, 'mixture_1_post_means.rds'))

# fit mixture 2 ---
fit_mixture_2 <- stan(file='mixture_2.stan', data=stan_data, chains=1, iter = 4e3, refresh = -1, verbose = FALSE)

mixture_2_post_means[j,i] <- get_posterior_mean(fit_mixture_2, pars = 'beta[1]')
write_rds(mixture_2_post_means,  file.path(out, 'mixture_2_post_means.rds'))


# fit t model -----
fit_t <- stan(file='t.stan', data=stan_data, chains=1, iter = 4e3, refresh = -1, verbose = FALSE)
stan_dens(fit_t)
t_post_means[j,i] <- get_posterior_mean(fit_t, pars = 'beta')
write_rds(t_post_means, file.path(out, 't_post_means.rds'))



# fit restricted likelihood Tukey ----
X <- matrix(rep(1, length(y)), ncol = 1)
fit_tukey <- restrictedBayesLm(y, X, regEst = 'Tukey', mu0 = eta, Sigma0 = tau^2, a0 = a, b0 = b, sigma2Int = var(y), nkeep = 4000, nburn = 2000)
tukey_post_means[j,i] <- mean(fit_tukey$mcmc[,1])
write_rds(tukey_post_means, file.path(out, 'tukey_post_means.rds'))


# fit restricted likelihood Huber ----
fit_huber <- restrictedBayesLm(y, X, regEst = 'Huber', mu0 = eta, Sigma0 = tau^2, a0 = a, b0 = b, sigma2Int = var(y),  nkeep = 4000, nburn = 2000)
huber_post_means[j,i] <- mean(fit_huber$mcmc[,1])
write_rds(huber_post_means, file.path(out, 'huber_post_means.rds'))

cat(paste(samp_sizes[i], j, '  at', Sys.time(), '\n'), file = file.path(out, 'log.txt'), append = TRUE)
  }
}


means_all <- bind_rows('Normal' = normal_post_means, 'Mixture 1' = mixture_1_post_means, 'Mixture 2' = mixture_2_post_means, 'Student-t' = t_post_means, 'Tukey' = tukey_post_means, 'Huber' = huber_post_means, .id = 'Fitted Model')
write_rds(means_all, file.path(out, 'means_all_wide.rds'))
means_all <- gather(means_all, n, mean, -'Fitted Model') %>% 
  group_by(`Fitted Model`, n) %>% 
  mutate(simulation = 1:n())
means_all
write_rds(means_all, file.path(out, 'means_all_long.rds'))
