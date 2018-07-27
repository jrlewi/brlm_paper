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
# true_model <- 't'


source('source.R')
#simulation parameters ----
n_sim <- 25
samp_sizes <- c(50, 100, 200)


# prior parameters ----
eta <- 0
tau <- 4
a <- 3
b <- 2
ngridpts <- 500 #for KL computations

iter <- 4000 #total length of MCMC iterations


for(true_model in c("Recording Error")){#"Normal", "Mixture 1", "Mixture 2",, "Student-t"))

p <- .1
  
# create data_gen function----
#generates true data
if(true_model == 'Normal'){
data_gen <- function(n){
  rnorm(n)
}
} else if(true_model == 'Mixture 1'){
  data_gen <- function(n){
    x <- rbinom(n, 1, 1-p)
    sapply(x, function(a) ifelse(a, rnorm(1), rnorm(1, sd = 3)))
  }
} else if(true_model == 'Mixture 2'){
  data_gen <- function(n){
  x <- rbinom(n, 1, 1-p)
  sapply(x, function(a) ifelse(a, rnorm(1), rnorm(1, mean = 3,  sd = 3)))
  }
} else if(true_model == 'Recording Error'){
  data_gen <- function(n){
  x <- rnorm(n)
  ind <- sample(1:length(x), size = floor(n*p))
  x[ind] <- 10*x
   # end <- floor(length(ind)/2)
   # x[ind[1:end]] <- -3
   # x[ind[(end+1):length(ind)]] <- 4
  # plot(x)
  x
  }
}





# set up empty tibbles to store results----
empty_mat <- matrix(NA, n_sim, length(samp_sizes))
colnames(empty_mat) <- samp_sizes
empty_mat <- as_tibble(empty_mat)

normal_post_means <- mixture_1_post_means <- mixture_2_post_means <- t_post_means <- tukey_post_means <- huber_post_means <- normal_KL <- mixture_1_KL <- mixture_2_KL <- t_KL <- tukey_KL <- huber_KL <- empty_mat 

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
fit_normal <- stan(file='standard_normal.stan', data=stan_data, chains=1, iter = iter, refresh = -1, verbose = FALSE)

samps <- rstan::extract(fit_normal, pars = c('beta','sigma2'))
normal_KL[j,i] <- compute_KL(thetaSamples = samps$beta, sigma2Samples = samps$sigma2, nuSamples = NULL, ngridpts = ngridpts,prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1)
write_rds(normal_KL,file.path(out, 'normal_KL.rds'))


normal_post_means[j,i] <- get_posterior_mean(fit_normal, pars = 'beta')
write_rds(normal_post_means,file.path(out, 'normal_post_means.rds'))

# fit mixture 1 ----
fit_mixture_1 <- stan(file='mixture_1.stan', data=stan_data, chains=1, iter = iter, refresh = -1, verbose = FALSE)

samps <- rstan::extract(fit_mixture_1, pars = c('beta','sigma2'))
mixture_1_KL[j,i] <- compute_KL(thetaSamples = samps$beta, sigma2Samples = samps$sigma2, nuSamples = NULL, ngridpts = ngridpts, prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1)
write_rds(mixture_1_KL,file.path(out, 'mixture_1_KL.rds'))


mixture_1_post_means[j,i] <- get_posterior_mean(fit_mixture_1, pars = 'beta')
write_rds(mixture_1_post_means, file.path(out, 'mixture_1_post_means.rds'))

# fit mixture 2 ---
fit_mixture_2 <- stan(file='mixture_2.stan', data=stan_data, chains=1, iter = iter, refresh = -1, verbose = FALSE)

samps <- rstan::extract(fit_mixture_2, pars = c('beta[1]','sigma2'))
mixture_2_KL[j,i] <- compute_KL(thetaSamples = samps$`beta[1]`, sigma2Samples = samps$sigma2, nuSamples = NULL, ngridpts = ngridpts, prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1)
write_rds(mixture_2_KL,file.path(out, 'mixture_2_KL.rds'))


mixture_2_post_means[j,i] <- get_posterior_mean(fit_mixture_2, pars = 'beta[1]')
write_rds(mixture_2_post_means,  file.path(out, 'mixture_2_post_means.rds'))


# fit t model -----
fit_t <- stan(file='t.stan', data=stan_data, chains=1, iter = iter, refresh = -1, verbose = FALSE)
stan_dens(fit_t)

samps <- rstan::extract(fit_t, pars = c('beta','sigma2', 'nu'))

t_KL[j,i] <- compute_KL(thetaSamples = samps$beta, sigma2Samples = samps$sigma2,nuSamples = samps$nu, ngridpts = ngridpts, prediction_dist = 't', thetaTrue = 0, sig2True = 1)
write_rds(t_KL,file.path(out, 't_KL.rds'))



t_post_means[j,i] <- get_posterior_mean(fit_t, pars = 'beta')
write_rds(t_post_means, file.path(out, 't_post_means.rds'))



# fit restricted likelihood Tukey ----
X <- matrix(rep(1, length(y)), ncol = 1)
fit_tukey <- restrictedBayesLm(y, X, regEst = 'Tukey', mu0 = eta, Sigma0 = tau^2, a0 = a, b0 = b, sigma2Int = var(y), nkeep = iter, nburn = floor(iter*0.5))

samps <- list(beta = as.numeric(fit_tukey$mcmc[,1]), sigma2 = as.numeric(fit_tukey$mcmc[,2]))
tukey_KL[j,i] <- compute_KL(thetaSamples = samps$beta, sigma2Samples = samps$sigma2, nuSamples = NULL, ngridpts = ngridpts, prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1)
write_rds(tukey_KL,file.path(out, 'tukey_KL.rds'))

tukey_post_means[j,i] <- mean(fit_tukey$mcmc[,1])
write_rds(tukey_post_means, file.path(out, 'tukey_post_means.rds'))


# fit restricted likelihood Huber ----
fit_huber <- restrictedBayesLm(y, X, regEst = 'Huber', mu0 = eta, Sigma0 = tau^2, a0 = a, b0 = b, sigma2Int = var(y),  nkeep = iter, nburn = floor(iter*0.5))

samps <- list(beta = as.numeric(fit_huber$mcmc[,1]), sigma2 = as.numeric(fit_huber$mcmc[,2]))
huber_KL[j,i] <- compute_KL(thetaSamples = samps$beta, sigma2Samples = samps$sigma2, nuSamples = NULL, ngridpts = ngridpts, prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1)
write_rds(huber_KL,file.path(out, 'huber_KL.rds'))


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


KL_all <- bind_rows('Normal' = normal_KL, 'Mixture 1' = mixture_1_KL, 'Mixture 2' = mixture_2_KL, 'Student-t' = t_KL, 'Tukey' = tukey_KL, 'Huber' = huber_KL, .id = 'Fitted Model')
write_rds(KL_all, file.path(out, 'KL_all_wide.rds'))
KL_all <- gather(KL_all, n, KL, -'Fitted Model') %>% 
  group_by(`Fitted Model`, n) %>% 
  mutate(simulation = 1:n())
KL_all
write_rds(KL_all, file.path(out, 'KL_all_long.rds'))
}
