library(tidyverse)


out <- read_rds(file.path(here::here(), 
                          '1_simulation_out_5.rds'))

nu <- 5
mse  <- out$estimates %>% filter(variable != 'sigma2') %>% 
  mutate(sq_error = (true_value - estimates)^2) %>% 
  group_by(simulation, prior_sd, method) %>% 
  summarize(mse_sim = mean(sq_error)) %>% 
  ungroup() %>% group_by(prior_sd, method) %>% 
  summarize(mse = mean(mse_sim), se = sd(mse_sim)/sqrt(n())) %>% ungroup() %>% 
  mutate(prior_sd = as_factor(as.character(prior_sd)))

ggplot(mse %>%  filter(prior_sd != 0.2), aes(prior_sd, mse, group = method, col = method)) +
  geom_point(position = position_dodge(width = .75)) + 
  geom_errorbar(mapping = aes(ymin = (mse - se), 
                              ymax = (mse + se)),
                width = 0.05, position  = position_dodge(width = .75), 
                linetype = 1) +
  theme_bw() +  theme(text = element_text(family = 'Times')) + 
  ylab('MSE') + xlab('Prior Standard Deviation') +
  guides(col=guide_legend(title="Method"))
ggsave(file.path(getwd(),
                 "..", "..","..", "figs", 
                 'mse_sim_many_p.png'), width = 6, height = 4)

#mse on just the active parameters

mse_sub  <- out$estimates %>% filter(variable %in% c('1', '2', '3'), prior_sd != 0.2) %>% 
  mutate(sq_error = (true_value - estimates)^2) %>% 
  group_by(simulation, prior_sd, method) %>% 
  summarize(mse_sim = mean(sq_error)) %>% 
  ungroup() %>% group_by(prior_sd, method) %>% 
  summarize(mse = mean(mse_sim), se = sd(mse_sim)/sqrt(n())) %>% ungroup() %>% mutate(prior_sd = as_factor(as.character(prior_sd)))

ggplot(mse_sub, aes(prior_sd, mse, group = method, col = method)) +
  geom_point(position = position_dodge(width = .75)) + 
  geom_errorbar(mapping = aes(ymin = (mse - se), 
                              ymax = (mse + se)),
                width = 0.05, position  = position_dodge(width = .75), 
                linetype = 1) +
  theme_bw() +  theme(text = element_text(family = 'Times')) + 
  ylab('MSE') + xlab('Prior Standard Deviation') +
  guides(col=guide_legend(title="Method"))
ggsave(file.path(getwd(),
                 "..", "..","..", "figs", 
                 'mse_sub_sim_many_p.png'), width = 6, height = 4)


ave_margs <- out$marginals %>% 
  group_by(prior_sd, method) %>% 
  summarise(mn = mean(mean_log_predictive_point),
            se = sd(mean_log_predictive)/sqrt(n())) %>% 
  ungroup() %>%  mutate(prior_sd = as_factor(as.character(prior_sd)))


#, method != "rlm"
ggplot(ave_margs %>% filter(prior_sd != 0.2), aes(prior_sd, -mn, group = method, 
                      col = method)) +
  geom_point(position = position_dodge(width = .75)) + 
  geom_errorbar(mapping = aes(ymin = (-mn - se), 
                              ymax = (-mn + se)),
                width = 0.05, position  = position_dodge(width = .75), 
                linetype = 1) + 
  theme_bw() +  theme(text = element_text(family = 'Times')) + 
  ylab('Ave. Neg. Log Predictive') + xlab('Prior Standard Deviation') +
  guides(col=guide_legend(title="Method"))
ggsave(file.path(getwd(),
                 "..", "..","..", "figs", 
                 'negll_sim_many_p.png'), width = 6, height = 4)


# Predictive MSE --------


# gen_good_data <- function(true_value, n, p, p_extra, sd_slopes, num_corr){
#   #true_value - vector c(betas, sigma2)
#   # n - number of samples to generate. 
#   X <- matrix(runif(n*p), n, p)
#   
#   if (num_corr > p_extra) { 
#     num_corr <- p_extra 
#   }
#   X_extra <- matrix(NA, n, num_corr)
#   which_active <- rep(1:p, length.out = num_corr)
#   for (jj in 1:num_corr) {
#     ind <- which_active[jj]
#     slope  <- rnorm(1, 0, sd = sd_slopes)
#     X_extra[,jj] <- X[,ind] * slope + rnorm(n, sd =  sd_slopes/3)
#   }
#   p_remain <-  p_extra - num_corr
#   if (p_remain > 0) {
#     X_extra <- cbind(X_extra, matrix(rnorm(n*p_remain), n,  p_remain))
#   }
#   
#   
#   beta <- true_value[1:p]
#   sigma <- sqrt(true_value[p + 1])
#   expected <- X %*% beta
#   y <- rnorm(n, expected, sigma)
#   cbind(y, X)
# }

# predict_full_model <- function(estimates, gen_out, 
#                       uncertainty = TRUE, model = 'Normal', nu = NULL){
#   #vector  of estimates c(betas, sigma2)
#   #gen_out: output of gen_good_data
#   p <- length(estimates) - 1
#   X <- gen_out[, -1]
#   y <- gen_out[, 1]
#   n <- length(y)
#   beta_hat <- estimates[1:p]
#   mn <- X %*% beta_hat
#   if(!uncertainty) {
#     preds <- mn 
#     } else {
#       sigma_hat <- sqrt(estimates[p+1])
#       if(model == 'Normal'){
#         error <- rnorm(n, 0, sigma_hat)
#       } else {
#         if(model == 't') {
#           error <- sigma_hat*rt(n, df = nu)
#         } else {
#           error('model must be Normal or t')
#         }
#       }
#       preds <- mn + error
#     }
#   as.numeric(preds)
# }


# predict_full_model <- function(estimates, data, X_extra, outlier, 
#                                uncertainty = TRUE, model = 'Normal', nu = NULL){
#   #vector  of estimates c(betas, sigma2)
#   #gen_out: output of gen_good_data
#   p <- length(estimates) - 1
#   X <- data[!outlier, -1]
#   y <- data[!outlier, 1]
#   n <- length(y)
#   beta_hat <- estimates[1:p]
#   X_full <- cbind(X, X_extra[!outlier, ])
#   mn <- X_full %*% beta_hat
#   if(!uncertainty) {
#     preds <- mn 
#   } else {
#     sigma_hat <- sqrt(estimates[p+1])
#     if(model == 'Normal'){
#       error <- rnorm(n, 0, sigma_hat)
#     } else {
#       if(model == 't') {
#         error <- sigma_hat*rt(n, df = nu)
#       } else {
#         error('model must be Normal or t')
#       }
#     }
#     preds <- mn + error
#   }
#   as.numeric(preds)
# }
# 
# 
# 
# 
# tmp_split <- out$estimates %>% 
#   split(list(.$simulation, .$prior_sd, .$method))
# 
# sim_one_pred <- function(split_df, n, nu){
#   #split_df - one element of list created above
#   true_value <- split_df$true_value
#   estimates <- split_df$estimates
#   gen_out <- gen_good_data(true_value, n)
#   y <- gen_out[, 1]
#   method <- split_df$method[1]
#   if(method == 'restricted' | method == 'rlm') model <- 'Normal'
#   if(method == 't') model <- 't'
#   
# preds <- predict_full_model(estimates, gen_out, uncertainty = FALSE, model = model, nu = nu) 
#   cbind(y, preds)
# }
# 
# n_preds <- 500
# 
# preds_list <- tmp_split %>% map(.f = function(df){
#   preds <- sim_one_pred(df, n = n_preds, nu = nu)
#   y <- preds[,1]
#   y_hat <- preds[,2]
#   simulation <- df$simulation[1]
#   prior_sd <- df$prior_sd[1]
#   method <- df$method[1]
#   sq_error <- (y - y_hat)^2
#   tibble(simulation = df$simulation[1],prior_sd = prior_sd, 
#          method = method, y = y, y_hat = y_hat, sq_error = sq_error)
# })
# preds <- bind_rows(preds_list)  
# 
# pmse <- preds %>% group_by(simulation, prior_sd, method) %>% 
#   summarize(ave_per_sim = mean(sq_error)) %>% 
#   ungroup() %>% group_by(prior_sd, method) %>% 
#   summarise(pmse = mean(ave_per_sim), se = sd(ave_per_sim)/sqrt(n())) %>% 
#   ungroup() %>% 
#   mutate(prior_sd = as_factor(as.character(prior_sd)))
# 
# 
# ggplot(pmse %>% filter(prior_sd != 0.2, method != 'rlm'), aes(prior_sd, pmse, group = method, col = method)) +
#   geom_point(position = position_dodge(width = .75)) + 
#   geom_errorbar(mapping = aes(ymin = (pmse - se), 
#                               ymax = (pmse + se)),
#                 width = 0.05, position  = position_dodge(width = .75), 
#                 linetype = 1) +
#   theme_bw() +  theme(text = element_text(family = 'Times')) + 
#   ylab('PMSE') + xlab('Prior Standard Deviation') +
#   guides(col=guide_legend(title="Method"))
# 
# ggsave(file.path(getwd(),
#                  "..", "..","..", "figs", 
#                  'pmse_sim_many_p.png'), width = 6, height = 4)
# 
# 
# 
# # split_df <- tmp_split$`1.0.8.t`
# # pred_tmp <- sim_one_pred(split_df, n = 10000, nu = 5)
# # summary((pred_tmp[,1] - pred_tmp[,2])^2)
# # 
# # n_pred <- 1000
# # true_value <- out$estimates$true_value[1:31]
# # estimates <- out$estimates$estimates[1:31]
# # tmp <- gen_good_data(true_value, n_pred)
# # preds_normal <- predict_full_model(estimates, tmp)
# # preds_t <- predict_full_model(estimates, tmp, model = 't', nu = 5)
# # summary((preds_normal[,1] - preds_normal[,2])^2 )
# # summary((preds_t[,1] - preds_t[,2])^2 )
# 
# # tmp <- data_generation(500, 3, 27, 5,1, .1, 25)
# # tmp$params
# # 
# # pairs(tmp$data[,1:4])
# # library(MASS)
# # fit_tmp <- rlm(cbind(tmp$data[,-1], tmp$X_extra), tmp$data[,1], psi = psi.bisquare)
# # summary(fit_tmp)

# acceptance rates-----

out$y_accept
ggplot(out$y_accept, aes(prior_sd, y_accept)) +
  geom_jitter() + theme_bw()



# mcmc chains and convergence metrics ----

dim(out$mcmc_samples)
abs(geweke.diag(out$mcmc_samples[1,1,,])$z)
all.equal(out$mcmc_samples[1,1,,], out$mcmc_samples[2,3,,])


geweke_test <- apply(out$mcmc_samples, c(1,2), function(x){
  abs(geweke.diag(x)$z)
})
dim(geweke_test)
apply(geweke_test, c(3), max)
plot(geweke_test[,,-1])
abline(h = 2)
mean(geweke_test[,,-1] <= 2)
1-2*(1 - pnorm(2))

mean(geweke_test[,,-1] <= 3)
1-2*(1 - pnorm(3))
