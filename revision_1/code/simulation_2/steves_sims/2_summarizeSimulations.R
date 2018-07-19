#
# Summarize the simulations ---
#
library('tidyverse')


# load the data, along with the factor levels, and true values of theta
data <- read_rds(file.path(getwd(),'data_sig2_4', 'data.rds'))
n_groups <- length(data$theta)
results_df <- file.path(getwd(), 'results')

ass.vec <- c(1.25, 2.5,5,10,20) 
scale_vec <- round(c(0.5, 1/sqrt(2),1, sqrt(2), 2),2)

sig2 <- 4

dfs <- vector('list', length(ass.vec)*length(scale_vec)*length(sig2))
i <- 1
for(ss in sig2){
  for(as in ass.vec){
    for(sc in scale_vec){
      rds_name <- paste0("__as_", as,"__scale_as_", sc, "__sig2_", ss, ".rds")
    
    huber <- read_rds(file.path(results_df, paste0("huber", rds_name)))
    theta_huber_rest <- apply(huber$theta, 2, mean)
    theta_huber_rlm  <- sapply(huber$robustFits, function(x) x$coefficients)
    huber_df <- tibble(theta = rep(data$theta,2), 
                       theta_hat = c(theta_huber_rest, theta_huber_rlm), 
                       statistic = rep('Huber', 2*n_groups), 
                       method = c(rep('restricted', n_groups), rep('rlm', n_groups)),
                       a_s = rep(as, 2*n_groups),
                       scale_as = rep(sc, 2*n_groups),
                       sigma2 = rep(ss, 2*n_groups)
                       )
    
    
    tukey <- read_rds(file.path(results_df, paste0("tukey", rds_name)))
    theta_tukey_rest <- apply(tukey$theta, 2, mean)
    theta_tukey_rlm  <- sapply(tukey$robustFits, function(x) x$coefficients)
    
    tukey_df <- tibble(theta = rep(data$theta,2), 
                       theta_hat = c(theta_tukey_rest, theta_tukey_rlm), 
                       statistic = rep('Tukey', 2*n_groups), 
                       method = c(rep('restricted', n_groups), rep('rlm', n_groups)),
                       a_s = rep(as, 2*n_groups),
                       scale_as = rep(sc, 2*n_groups),
                       sigma2 = rep(ss, 2*n_groups)
    )
    
    
    
    
    normal <- read_rds(file.path(results_df, paste0("normal", rds_name)))
    theta_normal <- apply(normal$theta, 2, mean)

    
    
    normal_df <- tibble(theta = rep(data$theta), 
                       theta_hat = c(theta_normal), 
                       statistic = rep('Normal', n_groups), 
                       method =  rep('Normal', n_groups),
                       a_s = rep(as, n_groups),
                       scale_as = rep(sc, n_groups),
                       sigma2 = rep(ss, n_groups)
    )
    
    
    df <- bind_rows(huber_df, tukey_df, normal_df)
    dfs[[i]] <- df
    i <- i+1
    }
  }
}

df_estimates <- bind_rows(dfs)



df_mse <- df_estimates %>% 
  group_by(statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(MSE = mean((theta - theta_hat)^2))



ggplot(df_mse, aes(x = as.factor(a_s), y = MSE, col = interaction(method,statistic), group = interaction(method,statistic))) +
  geom_line() + geom_point() +
  facet_wrap(~scale_as)
  


ggplot(df_mse, aes(x = as.factor(scale_as), y = MSE, col = interaction(method,statistic), group = interaction(method,statistic))) +
  geom_line() + geom_point() +
  facet_wrap(~a_s)


