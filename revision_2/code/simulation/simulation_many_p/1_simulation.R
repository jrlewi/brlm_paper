# Consider regression problem with many p, larger n,
# and outliers

library(tidyverse)
library(brlm)
library(MCMCpack)

#Simulation Parameters ----- 
n_sims <- 50
n <- 500 #sample size
p <- 3 #active covariates
p_extra <- 20 - p #inactive covariates
a_0 <- 5
b_0 <- 5 
prob_outlier <- .2
outlier_contam <- 10
mu_0 <-  rep(0, p+p_extra)
prior_sd <- seq(.2, 1, by = .2) 
# Sigma_0 = diag(p+p_extra),
nkeep <- 1000
nburn <- 2000
maxit <- 1000
nu <- 5 #df for t_fit


# helper functions -----
data_generation <- function(n, p, p_extra,
                               a_0,
                               b_0,
                               prob_outlier,
                               outlier_contam){
  X <- matrix(runif(n*p), n, p)
  beta <- rnorm(p, 0, 1)
  sigma2 <- rinvgamma(1, a_0, b_0)
  expected <- X %*% beta
  vr_out <- (outlier_contam)*sigma2
  y <- sapply(expected, function(mn){
    outlier <- rbinom(1,1, prob_outlier)
    err <- if (outlier) {
      rnorm(1, 0, sd = sqrt(vr_out))
    } else {
      rnorm(1, 0, sd = sqrt(sigma2))
    }
    list(mn + err, outlier)
  })
  X_extra <- matrix(runif(n*p_extra), n, p_extra)
  data <- cbind(unlist(y[1,]), X)
  params <- c(beta, sigma2)
  outlier <- unlist(y[2,])
  list(data = data, params = params, outlier = outlier, 
       X_extra = X_extra)
}

#marginal computations ----
get_marginal_rest <- function(y_gen, X_gen, mcmc){
  samps_coef <- t(mcmc[, 1:ncol(X_gen)])
  mns <- X_gen%*%samps_coef
  samps_sigma <- as.numeric(mcmc[, (ncol(X_gen) + 1)])^0.5
  sds <- matrix(samps_sigma, ncol = nkeep, nrow = nrow(X_gen),
                byrow = TRUE) 
  
  pdfs <- dnorm(y_gen, mns, sds)
  log(rowMeans(pdfs))
}

get_marginal_t <- function(y_gen, X_gen, mcmc){
  samps_coef <- t(mcmc[, 1:ncol(X_gen)])
  mns <- X_gen%*%samps_coef
  samps_sigma <- as.numeric(mcmc[, (ncol(X_gen) + 1)])^0.5
  sds <- matrix(samps_sigma, ncol = nkeep, nrow = nrow(X_gen),
                byrow = TRUE) 
  
  pdfs <- brlm::tdensity(y_gen, mns, sds, nu = nu)
  log(rowMeans(pdfs))
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
                  psi = psi.bisquare, maxit = maxit)
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
  
  estimates_all <- rbind(rlm_estimates, 
                         rest_estimates, 
                         t_estimates,
                         truth)
  
  #compute marginals for non-outliers
  y_gen <- y[!outlier]
  X_gen <- X_all[!outlier, ] 
  
  #rlm ----
  mn <- X_gen%*%rlm_estimates
  rlm_margninals <- 
    dnorm(y_gen, mn, sd = rlm_sig2_hat^.5, log = TRUE)
  
  #restricted ---
  rest_marginals <-
    get_marginal_rest(y_gen, X_gen, mcmc = rest_tukey$mcmc)
  
  #t-dist
  t_marginals <- get_marginal_t(y_gen, X_gen, mcmc = t_fit$mcmc)

  marginals_all <- rbind(rlm_margninals,
                         rest_marginals, 
                         t_marginals)
  
  list(estimates = estimates_all, 
       marginals = marginals_all, 
       sigma2 = params[p + 1], 
       y_accept = y_accept, 
       rest_mcmc = rest_tukey$mcmc)
}


# Start Simulation ------
I_mat <- diag(p + p_extra) 

rlm_estimates <- rest_estimates <- t_estimates <- 
  truth <-
  array(NA, c(n_sims, length(prior_sd), p + p_extra))

rlm_marginals <- rest_marginals <- 
  t_marginals <- y_accept <- 
  array(NA, c(n_sims, length(prior_sd)))

mcmc_samples <- array(NA, 
                      c(n_sims, length(prior_sd), 
                        nkeep, p + p_extra + 1))

strt <- Sys.time()

for(i in 1:n_sims){
  data <- data_generation(n, p,p_extra,
                          a_0,
                          b_0,
                          prob_outlier,
                          outlier_contam)
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

rlm_marginals[i,j] <- mean(result$marginals[1,])
rest_marginals[i,j] <- mean(result$marginals[2,])
t_marginals[i,j] <- mean(result$marginals[3,])

y_accept[i,j] <- result$y_accept

mcmc_samples[i,j,,] <- result$rest_mcmc
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
           variable = rep(1:n_variables, n_methods),
           estimates = c(rlm_estimates[i,j,],
                         rest_estimates[i,j,],
                         t_estimates[i,j,]),
           method = c(rep('rlm', n_variables),
                      rep('restricted', n_variables),
                      rep('t', n_variables)),
           true_value = rep(truth[i,j,], n_methods))
    
    marginals_list[[k]] <- tibble(simulation = i, 
            prior_sd = prior_sd[j],
            mean_log_predictive = c(
              rlm_marginals[i,j],
              rest_marginals[i,j],
              t_marginals[i,j]),
              method = c('rlm', 'restricted', 't')
            )
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
            mcmc_samples = mcmc_samples)

write_rds(out, 
          file.path(here::here(), 
                    '1_simulation_out.rds'))



# make plots ----

mse  <- out$estimates %>%
      mutate(sq_error = (true_value - estimates)^2) %>% 
      group_by(prior_sd, method) %>% 
      summarize(mse = mean(sq_error)) %>% 
      ungroup() %>% 
      mutate(prior_sd = as_factor(as.character(prior_sd)))
  
ggplot(mse, aes(prior_sd, mse, group = method, col = method)) +
  geom_line()



ave_margs <- out$marginals %>% 
  group_by(prior_sd, method) %>% 
  summarise(mean = mean(mean_log_predictive),
            se = sd(mean_log_predictive)/sqrt(n()))


ggplot(ave_margs, aes(prior_sd, mean, group = method, 
                      col = method)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0)
