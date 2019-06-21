library(tidyverse)
library(MCMCpack)
out <- read_rds(file.path(here::here(), 
                          '1_simulation_out.rds'))
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

mse_sub  <- out$estimates %>% filter(variable %in% c('1', '2', '3')) %>% 
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
ggplot(ave_margs, aes(prior_sd, -mn, group = method, 
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



# acceptance rates-----

range(out$y_accept$y_accept)
ggplot(out$y_accept, aes(prior_sd, y_accept)) +
  geom_jitter() + theme_bw()
ggsave(file.path(getwd(),
                 "..", "..","..", "figs", 
                 'yaccept_sim_many_p.png'), width = 6, height = 4)



# mcmc chains and convergence metrics ----

dim(out$mcmc_samples)
abs(geweke.diag(out$mcmc_samples[1,1,,])$z)


geweke_test <- apply(out$mcmc_samples, c(1,2), function(x){
  abs(geweke.diag(x)$z)
})
dim(geweke_test)
apply(geweke_test, c(3), max)
plot(geweke_test)
abline(h = 2)
mean(geweke_test <= 2)
1-2*(1 - pnorm(2))

mean(geweke_test[,,-1] <= 3)
1-2*(1 - pnorm(3))


plot(geweke_test[31,,])
abline(h = 2)
mean(geweke_test[31,,] <= 2)


matplot(out$mcmc_samples[5,1,,31], type = 'l')


#correlation ----

X1 <- out$data_save[[1]][,-1]
X2 <- out$X_extra_save[[1]]

cor_mat <- cor(cbind(X1,X2)) %>% round(2)
dim(cor_mat)
cor_mat <- cor_mat %>% as.data.frame()
names(cor_mat) <- rownames(cor_mat) <-  1:ncol(cor_mat)
cor_mat <- cor_mat %>% mutate(`variable 1` = as.character(1:ncol(cor_mat)))
cor_df <- cor_mat %>% gather(`variable 2`, correlation, -`variable 1`) %>% as_tibble()
cor_df <- cor_df %>% mutate_at(.vars = c(1,2), factor, levels = 1:ncol(cor_mat))

ggplot(cor_df, aes(x= `variable 1`, y= `variable 2`, fill = correlation)) + 
  geom_tile() + theme_bw() +  theme(text = element_text(family = 'Times')) + xlab('variable 1') + ylab('variable 2')
ggsave(file.path(getwd(),
                 "..", "..","..", "figs", 
                 'cor_plot_sim_many_p.png'), width = 6, height = 4)

