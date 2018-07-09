# summarize calls data and make plots
library(tidyverse)
library(MASS)
fig_path <- file.path('..', '..', 'figs')
log_scale  <- read_rds(path = file.path(here::here(), 'out', 'log_scale.rds'))
n_prior  <- read_rds(path = file.path(here::here(), 'out', 'n_prior.rds'))
phones
#phones$year<- phones$year - mean(phones$year)
phones <- bind_cols(phones)
if(log_scale) phones$calls <- log(phones$calls)
plot(phones$year, phones$calls, ylab = if_else(as.logical(log_scale), 'log calls', 'calls'))

x <- phones$year- mean(phones$year)
y <- phones$calls
rl_fit <- read_rds(file.path(here::here(), 'out', 'rl_fit.rds'))
t_fit <- read_rds(file.path(here::here(), 'out', 't_fit.rds'))
nu <- 5
normal_fit <- read_rds(file.path(here::here(), 'out', 'normal_fit.rds'))
mixture_fit <- read_rds(file.path(here::here(), 'out', 'mixture_fit.rds'))
# samples from 'good' part of model
mixture_fit_good <- cbind(mixture_fit$beta0[,1],mixture_fit$beta1[,1], mixture_fit$sigma[,1])
apply(mixture_fit_good, 2, mean)
rl_fit$fit
#plot(log(rl_fit$w + 1e-20), cex = .2)


#functions needed posteriors ----

#plot posteriors
plot_post_imp <- function(rl_fit, variable = 1, xlim = NULL){
  plot(density(rl_fit$impSamps[,variable], weights = rl_fit$w), xlim = xlim)
}

# random sample from ppd for rl_fit
rpd_rl <- function(rl_fit,x, n){
  # x = scalar, 
  # y_tilde = vector of ys to eval pdd 
  # n number of samples
  w <- rl_fit$w/sum(rl_fit$w)
  n_imps <- nrow(rl_fit$impSamps)
  betas <- rl_fit$impSamps[,1:2]
  means <-c(1, x)%*%t(betas)
  sig <- sqrt(rl_fit$impSamps[,3])
  ind <- sample(1:n_imps, size = n, prob = w, replace = TRUE)
  rnorm(n, means[ind], sig[ind])
}
apply(rl_fit$impSamps, 2, function(x) mean(sample(x, 10000, prob = rl_fit$w, replace = TRUE)))


rpd_normal <- function(normal_fit,x, n){
  # x = scalar, 
  # y_tilde = vector of ys to eval pdd 
  # n number of samples
  n_imps <- nrow(normal_fit$mcmc)
  betas <- normal_fit$mcmc[,1:2]
  means <-c(1, x)%*%t(betas)
  sig <- sqrt(normal_fit$mcmc[,3])
  ind <- sample(1:n_imps, size = n, replace = TRUE)
  rnorm(n, means[ind], sig[ind])
}

rpd_mixture <- function(mixture_fit_good,x, n){
  # x = scalar, 
  # y_tilde = vector of ys to eval pdd 
  # n number of samples
  n_imps <- nrow(mixture_fit_good)
  betas <- mixture_fit_good[,1:2]
  means <-c(1, x)%*%t(betas)
  sig <- mixture_fit_good[,3]
  ind <- sample(1:n_imps, size = n, replace = TRUE)
  rnorm(n, means[ind], sig[ind])
}
mean_mix <- apply(mixture_fit_good, 2, mean)
mean_mix
means_rl <- apply(rl_fit$impSamps, 2, function(x) mean(sample(x, 100000, prob = rl_fit$w, replace = TRUE)))
library(MASS)
robust <- rlm(phones$calls[-c(1:3)]~c(phones$year-mean(phones$year))[-c(1:3)], psi = psi.bisquare)

plot(x, y)
abline(a = mean_mix[1], b = mean_mix[2])
abline(a = means_rl[1], b = means_rl[2], col= 2)
abline(a = coef(robust)[1], b =coef(robust)[2], col= 3)

rpd_t <- function(t_fit,x, n){
  # x = scalar - year to predict
  # n number of samples
  n_mcmc <- nrow(t_fit$mcmc)
  betas <- t_fit$mcmc[,1:2]
  means <-c(1, x)%*%t(betas)
  sig <- sqrt(t_fit$mcmc[,3])
  ind <- sample(1:n_mcmc, size = n, replace = TRUE)
  rt(n, nu)*sig[ind] + means[ind]
}

pd_t <- function(t_fit,x, y_tilde){
  # x = scalar - year to predict
  # y_tilde = grid to fit ppd
  n_mcmc <- nrow(t_fit$mcmc)
  betas <- t_fit$mcmc[,1:2]
  means <-c(1, x)%*%t(betas)
  sig <- sqrt(t_fit$mcmc[,3])
  sapply(y_tilde, function(y){
    mean(dt((y - means)/sig, nu)/sig) 
})
}




# ppds and credible intervals -----
x_grid <- seq(min(x), max(x), length.out = 50)
n_samps <- 1e4
set.seed(123)
# rest. likelihood ----
ppd_rl <- sapply(x_grid, function(x){
  rpd_rl(rl_fit,x, n = n_samps)
})

ppd_rl <- as_tibble(as.data.frame(cbind(t(ppd_rl), x_grid)))

ppd_rl_gather <- gather(ppd_rl, key = 'sample', value = 'value', -x_grid)
length(unique(ppd_rl_gather$sample))

rl_cred_bounds <- ppd_rl_gather %>% 
  group_by(x_grid) %>% 
  summarize(lower = quantile(value, probs = .05), upper = quantile(value, probs = .95), Model = 'restricted likelihood')

# t Model -----
ppd_t <- sapply(x_grid, function(x){
  rpd_t(t_fit,x, n = n_samps)
})

ppd_t <- as_tibble(as.data.frame(cbind(t(ppd_t), x_grid)))

ppd_t_gather <- gather(ppd_t, key = 'sample', value = 'value', -x_grid)

t_cred_bounds <- ppd_t_gather %>% 
  group_by(x_grid) %>% 
  summarize(lower = quantile(value, probs = .05), upper = quantile(value, probs = .95), Model = 't')

# normal Model on subset of the data ----
ppd_normal <- sapply(x_grid, function(x){
  rpd_normal(normal_fit,x, n = n_samps)
})

ppd_normal <- as_tibble(as.data.frame(cbind(t(ppd_normal), x_grid)))

ppd_normal_gather <- gather(ppd_normal, key = 'sample', value = 'value', -x_grid)

normal_cred_bounds <- ppd_normal_gather %>% 
  group_by(x_grid) %>% 
  summarize(lower = quantile(value, probs = .05), upper = quantile(value, probs = .95), Model = 'Normal')


#Mixture Model ---
ppd_mixture <- sapply(x_grid, function(x){
  rpd_mixture(mixture_fit_good,x, n = n_samps)
})

ppd_mixture <- as_tibble(as.data.frame(cbind(t(ppd_mixture), x_grid)))

ppd_mixture_gather <- gather(ppd_mixture, key = 'sample', value = 'value', -x_grid)

mixture_cred_bounds <- ppd_mixture_gather %>% 
  group_by(x_grid) %>% 
  summarize(lower = quantile(value, probs = .05), upper = quantile(value, probs = .95), Model = 'Mixture')



# combine credible bounds -----
cred_bounds <- bind_rows(normal_cred_bounds, t_cred_bounds, mixture_cred_bounds, rl_cred_bounds) %>% 
  gather(key = 'quantile', value = 'y', -x_grid, -Model)

phones_df <- bind_cols(phones)
phones_df$Calls <- c(rep('Used for prior', n_prior), rep('Used for fit', nrow(phones_df) - n_prior))

ggplot(cred_bounds, aes(x = x_grid + mean(phones$year), y = y, group = interaction(quantile,Model), col = Model )) + geom_smooth(se = FALSE, lwd = .5) +
  geom_point(data = phones_df, mapping = aes(x = year, y = calls, group = NULL, col = NULL, shape = Calls)) + labs(x = 'year', y = 'log calls') + theme_bw() + scale_color_brewer(palette = 'Set2') #+ geom_abline(slope =  mean_mix[2]  , intercept = mean_mix[1] - mean_mix[2]*mean(phones$year))
ggsave(file.path(fig_path, 'calls_predictive.png'))

# # lengths of credible intervals
# cred_length <- spread(cred_bounds, quantile, y) %>%  mutate(difference = upper - lower)
# 
# ggplot(filter(cred_length, Model == "restricted likelihood"), aes(x = x_grid, y = difference)) + geom_smooth(se = FALSE)
#   
# ggplot(filter(cred_length, Model == "t"), aes(x = x_grid, y = difference)) + geom_smooth(se = FALSE)    
# 
# 



# y_tilde<- seq(-150, 150, length.out = 500)
# ppd_t1 <- sapply(x_grid, function(x){
#   pd_t(t_fit,x, y_tilde)
# })
# 
# areas <- apply(ppd_t1, 2, function(y) cumsum(y[-1]*diff(y_tilde)))
# inds_lower <- apply(areas, 2, function(y) which(abs(y-.05) == min(abs(y-.05))))
# inds_upper <- apply(areas, 2, function(y) which(abs(y-.95) == min(abs(y-.95))))
# lowers_t <- y_tilde[inds_lower]
# uppers_t <- y_tilde[inds_upper]
# 
# t1_cred_bounds <- tibble(x_grid = x_grid, lower = lowers_t, upper = uppers_t, Model = 't1')