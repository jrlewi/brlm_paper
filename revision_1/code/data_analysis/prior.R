# Script to set up prior distributions ----

library(MASS)
library(tidyverse)


prior_data <- read_rds(file.path(here::here(), 'data', 'prior_data.rds'))
prior_data <- prior_data %>% 
  mutate(sqrt_count_2008 = sqrt(Count_2008), sqrt_count_2010 = sqrt(Count_2010))


# pooled regression analysis ----
prior_fit <- MASS::rlm(sqrt_count_2010 ~ sqrt_count_2008 - 1, scale.est = 'Huber', data =  prior_data)

theme_set(theme_bw(base_family = 'Times'))
ggplot(prior_data, aes(x = sqrt_count_2008, y = sqrt_count_2010, col = Type)) + geom_point(size = .3) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x -1, method.args = list(scale.est = 'Huber'), size = .5, lty = 2, se = FALSE, col = 1) + guides(col = guide_legend(title = 'Agency Type'))


summary(prior_fit)
beta_0 <- coef(prior_fit)
se_beta_0 <- vcov(prior_fit)^.5
sigma2_hat <- prior_fit$s^2

parms_prior <- list(beta_0 = beta_0, se_beta_0 = se_beta_0, sigma2_hat = sigma2_hat)
write_rds(parms_prior, "parms_prior.rds")

#hierarchcical regression analysis -----

by_state <- prior_data %>% 
  group_by(State) %>% 
  filter(n() >= 3) %>% 
  nest() 
  
state_model <- function(df){
  MASS::rlm(sqrt_count_2010 ~ sqrt_count_2008 - 1, scale.est = 'Huber', data = df, maxit = 100)
}

models <- by_state$data %>% 
  map(state_model)

state_slopes <-
  models %>% 
  map(coef) %>%
  unlist()

plot(state_slopes)


