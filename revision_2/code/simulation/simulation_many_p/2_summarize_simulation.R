library(tidyverse)


out <- read_rds(file.path(here::here(), 
                    '1_simulation_out.rds'))


mse  <- out$estimates %>% filter(variable != 'sigma2') %>% 
  mutate(sq_error = (true_value - estimates)^2) %>% 
  group_by(simulation, prior_sd, method) %>% 
  summarize(mse_sim = mean(sq_error)) %>% 
  ungroup() %>% group_by(prior_sd, method) %>% 
  summarize(mse = mean(mse_sim), se = sd(mse_sim)/sqrt(n())) %>% ungroup() %>% 
  mutate(prior_sd = as_factor(as.character(prior_sd)))

ggplot(mse, aes(prior_sd, mse, group = method, col = method)) +
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

  


ave_margs <- out$marginals %>% 
  group_by(prior_sd, method) %>% 
  summarise(mn = mean(mean_log_predictive),
            se = sd(mean_log_predictive)/sqrt(n())) %>% 
  ungroup() %>%  mutate(prior_sd = as_factor(as.character(prior_sd)))



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
out$y_accept
ggplot(out$y_accept, aes(prior_sd, y_accept)) +
  geom_jitter()



#mcmc chains and convergence metrics ----

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
