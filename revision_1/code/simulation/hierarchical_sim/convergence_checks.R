#convergence checks -----

library(MCMCpack)
library(tidyverse)



sims <- 1:10

ass.vec <- c(1.25, 5, 10) #c(1.25, 2.5,5,10,20) 
scale_vec <- c(0.5,1, 2) #round(c(0.5, 1/sqrt(2),1, sqrt(2), 2),2)

ss <- 4 #sigma2
convergence_scores <- array(NA, c(length(sims),3,  length(ass.vec), length(scale_vec), 182)) #sims, models, ass.vec, scale_vec, parameters

for(data_sim in sims){
  # load the data, along with the factor levels, and true values of theta
  data <- read_rds(file.path(getwd(),'data_sig2_4', paste0('data_', data_sim, '.rds')))
  n_groups <- length(data$theta)
  results_df <- file.path(getwd(), 'results')

    for(as in ass.vec){
      for(sc in scale_vec){
        rds_name <- paste0("_",data_sim, "__as_", as,"__scale_as_", sc, "__sig2_", ss, ".rds")
        
        fit <- read_rds(file.path(results_df, paste0("huber", rds_name)))
        th <-  geweke.diag(mcmc(fit$theta))$z
        sig2 <-  geweke.diag(mcmc(fit$sigma2))$z
        tauMu <- geweke.diag(mcmc(cbind(fit$mu, fit$tau2)))$z
        convergence_scores[data_sim, 1, which(as == ass.vec), which(sc == scale_vec),] <- c(th, sig2, tauMu)
        
        fit <- read_rds(file.path(results_df, paste0("tukey", rds_name)))
        th <-  geweke.diag(mcmc(fit$theta))$z
        sig2 <-  geweke.diag(mcmc(fit$sigma2))$z
        tauMu <- geweke.diag(mcmc(cbind(fit$mu, fit$tau2)))$z
        convergence_scores[data_sim, 2, which(as == ass.vec), which(sc == scale_vec),] <- c(th, sig2, tauMu)

        
        fit <- read_rds(file.path(results_df, paste0("normal", rds_name)))
        th <-  geweke.diag(mcmc(fit$theta))$z
        sig2 <-  geweke.diag(mcmc(fit$sigma2))$z
        tauMu <- geweke.diag(mcmc(cbind(fit$mu, fit$tau2)))$z
        convergence_scores[data_sim, 3, which(as == ass.vec), which(sc == scale_vec),] <- c(th, sig2, tauMu)
      }
    }
  print(data_sim)
}


plot(convergence_scores, pch = 19, cex = .3)
mean(abs(convergence_scores)>2)
2*(1-pnorm(2))
range(convergence_scores[,1,,,])
range(convergence_scores[,2,,,])
range(convergence_scores[,3,,,])
