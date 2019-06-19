# Consider regression problem with many p, larger n,
# and outliers

library(tidyverse)
library(brlm)
library(MCMCpack)
set.seed(123)
#Simulation Parameters ----- 
n_sims <- 30
n <- 500 #sample size
p <- 3 #active covariates
beta <- rep(1, p)
sigma2 <- 2
p_extra <- 30 - p #inactive covariates
num_corr <- 6 # number of extra covariates to correlate with the active ones
a_0 <- 5
b_0 <- 8 
sd_cor_1 <- 2 # sd for cor between the x's
sd_cor_2 <- 1 # sd for cor between true x's and extra x's

prob_outlier <- .2
outlier_contam <- 25
mu_0 <-  rep(0, p + p_extra)
prior_sd <- seq(.4, 1.4, by = .2) # c(.4, 1) # 
# Sigma_0 = diag(p+p_extra),
nkeep <- 1000
nburn <- 1500
maxit <- 1000
nu <- 5 #df for t_fit


# helper functions -----
data_generation <- function(n, p, p_extra,
                               beta = NULL,
                               sigma2 = NULL,
                               a_0,
                               b_0,
                               prob_outlier,
                               outlier_contam, 
                               sd_cor_1,
                               sd_cor_2,
                               num_corr){
  #num_corr is the number of X_extras to correlate with one of the active covariates
  
  if(is.null(beta)){
  beta <- rnorm(p, 0, 1)
  } else {
    p <- length(beta)
  }
  if(is.null(sigma2)){
  sigma2 <- rinvgamma(1, a_0, b_0)
  }
  
  
  
  X <- matrix(NA, n, p)
  x1 <- rnorm(n)
  X[,1] <- x1
  for(ii in 2:p){
    X[,ii] <- x1 + rnorm(n, sd = sd_cor_1) 
  }
  # consider some correlation structure for the extra covariates
  #n_groups <- floor(p_extra/p)
  #extra <- lapply(1:n_groups, function(a){
  
  if (num_corr > p_extra) { 
    num_corr <- p_extra 
    }
  X_extra <- matrix(NA, n, num_corr)
which_active <- rep(1:p, length.out = num_corr)
  for (jj in 1:num_corr) {
  ind <- which_active[jj]
  slope  <- 1 #rnorm(1, 0, sd = sd_slopes)
  X_extra[,jj] <- X[,ind] * slope + rnorm(n, sd =  sd_cor_2)
  }
  p_remain <-  p_extra - num_corr
  if (p_remain > 0) {
  X_extra <- cbind(X_extra, matrix(rnorm(n*p_remain), n,  p_remain))
  }
  # cor(X_extra[,5], X[, 2])
  # X_extra <- matrix(runif(n*p_extra), n, p_extra)
  # X <- matrix(rnorm(n*p), n, p)
  # X_extra <- matrix(rnorm(n*p_extra), n, p_extra)
  # X_all <- cbind(X, X_extra)
  # X_all_corr <- 
  
  X <- scale(X)
  X_extra <- scale(X_extra) #mean zero, variance 1
  
  #cor(X_extra[,5], X[, 2])
  
  expected <- X %*% beta
  vr_out <- (outlier_contam)*sigma2
  
  y <- sapply(expected, function(mn){
    outlier <- rbinom(1,1, prob_outlier)
    err <- if (outlier) {
      abs(rnorm(1, 0, sd = sqrt(vr_out)))
    } else {
      rnorm(1, 0, sd = sqrt(sigma2))
    }
    c(mn + err, outlier)
  })
  # X_extra <- matrix(runif(n*p_extra), n, p_extra)
  data <- cbind(y[1,], X)
  params <- c(beta, sigma2)
  outlier <- y[2,]
  #pairs(data, col = outlier + 1)
  list(data = data, params = params, outlier = outlier, 
       X_extra = X_extra)
}


# tmp <- data_generation(n, p, p_extra,
#                             beta = beta,
#                             sigma2 = sigma2,
#                             a_0,
#                             b_0,
#                             prob_outlier,
#                             outlier_contam,
#                             sd_cor_1,
#                             sd_cor_2,
#                             num_corr)
# # # 
# # pairs(tmp$data, col = tmp$outlier + 1)



#marginal computations ----
get_marginal_rest <- function(y_gen, X_gen, mcmc, point_est = FALSE){
  samps_coef <- t(mcmc[, 1:ncol(X_gen)])
  samps_sigma <- as.numeric(mcmc[, (ncol(X_gen) + 1)])^0.5
  if(!point_est){
  mns <- X_gen %*% samps_coef
  sds <- matrix(samps_sigma, ncol = nkeep, nrow = nrow(X_gen),
                byrow = TRUE) 
  
  pdfs <- dnorm(y_gen, mns, sds)
  out <- log(rowMeans(pdfs))
  } 
  if(point_est){
    post_mean_betas <- apply(samps_coef, 1, mean)
    post_mn <- X_gen %*% post_mean_betas
    post_mean_sigma <- mean(samps_sigma)
    out <- log(dnorm(y_gen, post_mn, post_mean_sigma))
  }
  out
}

get_marginal_t <- function(y_gen, X_gen, mcmc, point_est = FALSE, nu){
  samps_coef <- t(mcmc[, 1:ncol(X_gen)])
  samps_sigma <- as.numeric(mcmc[, (ncol(X_gen) + 1)])^0.5
  if(!point_est){
  mns <- X_gen%*%samps_coef
  sds <- matrix(samps_sigma, ncol = nkeep, nrow = nrow(X_gen),
                byrow = TRUE) 
  
  pdfs <- brlm::tdensity(y_gen, mns, sds, nu = nu)
  out <- log(rowMeans(pdfs))
  }
  if(point_est){
    post_mean_betas <- apply(samps_coef, 1, mean)
    post_mn <- X_gen %*%  post_mean_betas
    post_mean_sigma <- mean(samps_sigma)
    out <- log(brlm::tdensity(y_gen, post_mn, post_mean_sigma, nu = nu))
  }
  as.numeric(out)
}


one_sim <- function(
                    a_0,
                    b_0,
                    mu_0 = rep(0, p + p_extra),
                    Sigma_0 = diag(p + p_extra),
                    nkeep,
                    nburn,
                    maxit,
                    nu,
                    data){
  #p number of active parameters
  #p_extra number of non-active
  
  y <- data$data[,1]
  X <- data$data[,-1]
  p <- ncol(X)
  n <- length(y)
  outlier <- data$outlier
  params <- data$params
  
  #random generation of extra covariates
  X_extra <- data$X_extra
  p_extra <- ncol(X_extra)
  X_all <- cbind(X, X_extra)
  
  #model fits -----
  
  #rlm_tukey ---
  rlm_tukey <- rlm(X_all, y, 
                  psi = psi.bisquare, scale.est = 'Huber', maxit = maxit)
  rlm_estimates <- coef(rlm_tukey) 
  rlm_sig2_hat <- rlm_tukey$s^2
  sigma2Int <- rlm_sig2_hat
  #restricted tukey
  rest_tukey <- brlm::restrictedBayesLm(y, 
                          X_all, 
                          regEst = 'Tukey',
                          scaleEst = 'Huber',
                          mu0 = mu_0,
                          Sigma0 = Sigma_0,
                          a0 = a_0,
                          b0 = b_0,
                          sigma2Int = sigma2Int,
                          nkeep = nkeep,
                          nburn = nburn,
                          maxit = maxit)
  
  rest_estimates <- colMeans(rest_tukey$mcmc)   
  rest_sig2_est <- rest_estimates[p+p_extra+1]
  rest_estimates <- rest_estimates[1:(p + p_extra)]
  y_accept <- mean(rest_tukey$yAccept)
  
  t_fit <- brlm::bayesTdistLm(y, 
                               X_all,
                               mu0 = mu_0,
                               Sigma0 = Sigma_0,
                               a0 = a_0,
                               b0 =  ((nu-2)/nu)*b_0,
                               nu = nu,
                               nkeep = nkeep,
                               nburn = nburn)
  t_estimates <- colMeans(t_fit$mcmc)   
  t_sig2_est <- t_estimates[p+p_extra+1]
  t_estimates <- t_estimates[1:(p + p_extra)]
  
  truth <- c(params[1:p], rep(0, p_extra))
  truth_sigma2 <- params[p+1]
  estimates_all <- rbind(rlm_estimates, 
                         rest_estimates, 
                         t_estimates,
                         truth)
  estimates_sigma2 <- c(rlm_sig2_hat,
                        rest_sig2_est,
                        t_sig2_est,
                        truth_sigma2)
  
  
  #compute marginals for non-outliers
  y_gen <- y[!outlier]
  X_gen <- X_all[!outlier, ] 
  
  #rlm ----
  mn <- X_gen %*% rlm_estimates
  rlm_margninals <- 
    dnorm(y_gen, as.numeric(mn), sd = rlm_sig2_hat^.5, log = TRUE)
  
  #restricted ---
  rest_marginals <-
    get_marginal_rest(y_gen, X_gen, mcmc = rest_tukey$mcmc)
  
  rest_marginals_point <-
    get_marginal_rest(y_gen, X_gen, mcmc = rest_tukey$mcmc, point_est = TRUE)
  
  #t-dist
  t_marginals <- get_marginal_t(y_gen, X_gen, mcmc = t_fit$mcmc, nu = nu)
  t_marginals_point <- get_marginal_t(y_gen, X_gen, mcmc = t_fit$mcmc, point_est = TRUE, nu = nu)
  
  
  marginals_all <- rbind(rlm_margninals, 
                         rest_marginals, 
                         rest_marginals_point,
                         t_marginals, 
                         t_marginals_point)
  
  list(estimates = estimates_all, 
       estimates_sigma2 = estimates_sigma2,
       marginals = marginals_all, 
       y_accept = y_accept, 
       rest_mcmc = rest_tukey$mcmc,
       t_mcmc = t_fit$mcmc)
}


# Start Simulation ------
I_mat <- diag(p + p_extra) 

rlm_estimates <- rest_estimates <- t_estimates <- 
  truth <-
  array(NA, c(n_sims, length(prior_sd), p + p_extra))

rlm_marginals <- rest_marginals <- rest_marginals_point <- 
  t_marginals <- t_marginals_point <- y_accept <- 
  array(NA, c(n_sims, length(prior_sd)))

mcmc_samples <- mcmc_samples_t <- array(NA, 
                      c(n_sims, length(prior_sd), 
                        nkeep, p + p_extra + 1))

estimates_sigma2 <- array(NA, c(n_sims, length(prior_sd), 4))

strt <- Sys.time()
data_save <- outlier_save <- X_extra_save <- vector('list', n_sims)
for(i in 1:n_sims) {
  data <- data_generation(n, p,p_extra,
                          beta, 
                          sigma2,
                          a_0,
                          b_0,
                          prob_outlier,
                          outlier_contam, 
                          sd_cor_1,
                          sd_cor_2,
                          num_corr)
  
  data_save[[i]] <- data$data
  X_extra_save[[i]] <- data$X_extra
  outlier_save[[i]] <- data$outlier
  

  for(j in seq_along(prior_sd)){

    p_sd <- prior_sd[j]
    prior_cov <- p_sd^2 * I_mat  
#strt <- Sys.time()
    result <- one_sim(
          a_0,
          b_0,
          mu_0 = mu_0,
          Sigma_0 = prior_cov,
          nkeep = nkeep,
          nburn = nburn,
          maxit = maxit,
          nu, 
          data)
#Sys.time() - strt
rlm_estimates[i,j,] <- result$estimates[1,]
rest_estimates[i,j,] <- result$estimates[2,]
t_estimates[i,j,] <- result$estimates[3,]
truth[i,j,] <- result$estimates[4,]

estimates_sigma2[i,j,] <- result$estimates_sigma2


rlm_marginals[i,j] <- mean(result$marginals[1,])
rest_marginals[i,j] <- mean(result$marginals[2,])
rest_marginals_point[i,j] <- mean(result$marginals[3,])
t_marginals[i,j] <- mean(result$marginals[4,])
t_marginals_point[i,j] <- mean(result$marginals[5,])

y_accept[i,j] <- result$y_accept

mcmc_samples[i,j,,] <- result$rest_mcmc
mcmc_samples_t[i,j,,] <- result$t_mcmc
  print(j)
  }
  
  print(i)
}

Sys.time() - strt

# Tibblize and Save output -----

n_variables <- (p + p_extra)
n_methods <- 3

estimates_list <- marginals_list <- 
  vector('list', n_sims*length(prior_sd))  
k <- 1
for(i in 1:n_sims){
  for(j in seq_along(prior_sd)){
    
    estimates_list[[k]] <- tibble(simulation = i, 
           prior_sd = prior_sd[j], 
           variable = rep(c(1:n_variables, 'sigma2'), n_methods),
           estimates = c(rlm_estimates[i,j,], estimates_sigma2[i,j,1],
                         rest_estimates[i,j,], estimates_sigma2[i,j,2],
                         t_estimates[i,j,], estimates_sigma2[i,j,3]),
           method = c(rep('rlm', n_variables + 1),
                      rep('restricted', n_variables + 1),
                      rep('t', n_variables + 1)),
           true_value = rep(c(truth[i,j,], estimates_sigma2[i,j,4]), n_methods))
    
    marginals_list[[k]] <- tibble(simulation = i, 
            prior_sd = prior_sd[j],
            mean_log_predictive = c(
              rlm_marginals[i,j],
              rest_marginals[i,j],
              t_marginals[i,j]),
            mean_log_predictive_point = c(
              rlm_marginals[i,j],
              rest_marginals_point[i,j],
              t_marginals_point[i,j]),
              method = c('rlm', 'restricted', 't'))
      k <- k + 1    
  }
}

estimates <- bind_rows(estimates_list)
marginals <- bind_rows(marginals_list)

colnames(y_accept) <- prior_sd

y_accept_tibble <- bind_cols(
  tibble(simulation = 1:n_sims),
  as_tibble(y_accept)) %>% 
  gather(prior_sd, y_accept,  -simulation)

out <- list(estimates = estimates, marginals = marginals,
            y_accept = y_accept_tibble, 
            mcmc_samples = mcmc_samples, 
            mcmc_samples_t = mcmc_samples_t,
            mcmc_samples =  mcmc_samples,
            mcmc_samples_t = mcmc_samples_t,
            data_save = data_save,
            X_extra_save = X_extra_save,
            outlier_save = outlier_save
            )

write_rds(out,
          file.path(here::here(),
                    '1_simulation_out_5.rds'))



# make plots ----

mse  <- out$estimates %>% filter(variable != 'sigma2') %>% 
      mutate(sq_error = (true_value - estimates)^2) %>%
      group_by(prior_sd, method) %>%
      summarize(mse = mean(sq_error)) %>%
      ungroup() %>%
      mutate(prior_sd = as_factor(as.character(prior_sd)))

mse  <- out$estimates %>% filter(variable != 'sigma2') %>% 
  mutate(sq_error = (true_value - estimates)^2) %>% 
  group_by(simulation, prior_sd, method) %>% 
  summarize(mse_sim = mean(sq_error)) %>% 
  ungroup() %>% group_by(prior_sd, method) %>% 
  summarize(mse = mean(mse_sim), se = sd(mse_sim)/sqrt(n())) %>% ungroup() %>% 
  mutate(prior_sd = as_factor(as.character(prior_sd)))

ggplot(mse, aes(prior_sd, mse, group = method, col = method)) +
  geom_line()



mse_sub <- out$estimates %>% filter(variable %in% c('1', '2', '3')) %>% 
  mutate(sq_error = (true_value - estimates)^2) %>%
  group_by(prior_sd, method) %>%
  summarize(mse = mean(sq_error)) %>%
  ungroup() %>%
  mutate(prior_sd = as_factor(as.character(prior_sd)))
ggplot(mse_sub, aes(prior_sd, mse, group = method, col = method)) +
  geom_line()



ave_margs <- out$marginals %>% 
  group_by(prior_sd, method) %>%
  summarise(mean = mean(mean_log_predictive_point),
            se = sd(mean_log_predictive)/sqrt(n()))


ggplot(ave_margs, aes(prior_sd, mean, group = method,
                      col = method)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0)
