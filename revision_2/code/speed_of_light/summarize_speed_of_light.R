# summarize the speed of light fits
library(tidyverse)
fig_path <- file.path('..', '..', 'figs')
nu <- 5

# read results ---- 
fit_tukey <- readRDS('out/fit_tukey.rds')
fit_huber <- readRDS('out/fit_huber.rds')
fit_lms <- readRDS('out/fit_lms.rds')
fit_lts <- readRDS('out/fit_lts.rds')
fit_normal <- readRDS('out/fit_normal.rds')
fit_t <- readRDS('out/fit_t.rds')
fit_restricted_t <- readRDS("out/fit_restricted_t.rds")
fit_mixture <- readRDS('out/fit_mixture.rds')

apply(fit_t, 2, mean)
apply(fit_t, 2, sd)

apply(fit_restricted_t, 2, mean)
apply(fit_restricted_t, 2, sd)


extract_pdf <- function(fit, model){
beta <- tibble(beta = fit$muPost[,1], posterior = fit$muPost[,2], Model = model)
sigma2 <- tibble(sigma2 = fit$sigma2Post[,1], posterior = fit$sigma2Post[,2], Model = model)
list(beta = beta, sigma2 = sigma2)
}

pdf_tukey <- extract_pdf(fit_tukey, 'Tukey')
pdf_huber <- extract_pdf(fit_huber, 'Huber')
pdf_lms <- extract_pdf(fit_lms, 'LMS')
pdf_lts <- extract_pdf(fit_lts, 'LTS')



compute_pdf <- function(samples, model){
  pdf_beta <- density(samples[,1])
  pdf_sigma2 <- density(samples[,2])
  beta <- tibble(beta = pdf_beta$x, posterior = pdf_beta$y, Model = model)
  sigma2 <- tibble(sigma2 = pdf_sigma2$x, posterior = pdf_sigma2$y, Model = model)
  list(beta = beta, sigma2 = sigma2)
}

pdf_normal <- compute_pdf(fit_normal, 'Normal')
pdf_t <- compute_pdf(fit_t, 't')
fit_mix_good <- cbind(fit_mixture$beta, fit_mixture$sigma2)
pdf_mixture <- compute_pdf(fit_mix_good, 'Mixture')
pdf_restricted_t <- compute_pdf(fit_restricted_t, 'Restricted_t')


model_order <- c('Normal', 't', 'Mixture', 'Tukey', 'Huber', 'LMS', 'LTS', 'Restricted_t')
beta_pdfs <- bind_rows(pdf_normal$beta, pdf_t$beta, pdf_mixture$beta, 
                       pdf_tukey$beta, pdf_huber$beta, pdf_lms$beta, pdf_lts$beta,
                       pdf_restricted_t$beta)
beta_pdfs$Model <- factor(beta_pdfs$Model, levels = model_order)


sigma2_pdfs <- bind_rows(pdf_normal$sigma2, pdf_t$sigma2, pdf_mixture$sigma2, 
                         pdf_tukey$sigma2, pdf_huber$sigma2, pdf_lms$sigma2, pdf_lts$sigma2,
                         pdf_restricted_t$sigma2)
sigma2_pdfs$Model <- factor(sigma2_pdfs$Model, levels = model_order)


# plot beta ----

ggplot(beta_pdfs, aes(x = beta, y = posterior, group = Model, col = Model))  + 
  geom_line() + theme_bw() + scale_color_brewer(palette = 'Set2') + 
  xlab(expression(beta)) + ylab('Posterior')#+ geom_hline(yintercept = 0, col = 'black')

ggsave(file.path(fig_path, 'speed_of_light_beta.png'))


# plot sigma2 ----

ggplot(sigma2_pdfs, aes(x = sigma2, y = posterior, group = Model, col = Model))  + geom_line() + theme_bw() + scale_color_brewer(palette = 'Set2') + xlab(expression(sigma^2)) + ylab('Posterior') #+ scale_x_log10() #+ geom_hline(yintercept = 0, col = 'black')

ggsave(file.path(fig_path, 'speed_of_light_sigma2.png'))


# Predictive Distributions ----


predictive_t <- function(fit_t, model){
  mu <- as.numeric(fit_t[,1])
  sigma2 <- as.numeric(fit_t[,2])
  y_tilde <- seq(-50, 100, length.out = 100)
  pred_dist <- sapply(y_tilde, function(yy){
    mean(dt((yy-mu)/sqrt(sigma2), df = nu)/sqrt(sigma2))
  })
  tibble(y_tilde = y_tilde, pred_dist = pred_dist, Model = model)
}

post_pred_t <- predictive_t(fit_t, model='t')
plot(post_pred_t$y_tilde, log(post_pred_t$pred_dist), type = 'l')
sum(diff(post_pred_t$y_tilde)*post_pred_t$pred_dist[-1])
sum(diff(post_pred_t$y_tilde)*post_pred_t$y_tilde[-1]*post_pred_t$pred_dist[-1])


post_pred_restricted_t <- predictive_t(fit_restricted_t, model= "Restricted_t")
plot(post_pred_restricted_t$y_tilde, log(post_pred_restricted_t$pred_dist), type = 'l')
sum(diff(post_pred_restricted_t$y_tilde)*post_pred_restricted_t$pred_dist[-1])
sum(diff(post_pred_restricted_t$y_tilde)*post_pred_restricted_t$y_tilde[-1]*post_pred_restricted_t$pred_dist[-1])


predictive_normal <- function(fit_normal){
  mu <- as.numeric(fit_normal[,1])
  sigma2 <- as.numeric(fit_normal[,2])
  y_tilde <- seq(-50, 100, length.out = 100)
  pred_dist <- sapply(y_tilde, function(yy){
    mean(dnorm(yy, mu, sqrt(sigma2)))
  })
  tibble(y_tilde = y_tilde, pred_dist = pred_dist, Model = 'Normal')
}


post_pred_normal <- predictive_normal(fit_normal)
plot(post_pred_normal$y_tilde, log(post_pred_normal$pred_dist))
sum(diff(post_pred_normal$y_tilde)*post_pred_normal$pred_dist[-1])

sum(diff(post_pred_normal$y_tilde)*post_pred_normal$y_tilde[-1]*post_pred_normal$pred_dist[-1])



# using only the 'good' piece. i.e. the small variance component
predictive_mixture <- function(fit_mix_good){
  mu <- as.numeric(fit_mix_good[,1])
  sigma2 <- as.numeric(fit_mix_good[,2])
  y_tilde <- seq(-50, 100, length.out = 100)
  pred_dist <- sapply(y_tilde, function(yy){
    mean(dnorm(yy, mu, sqrt(sigma2)))
  })
  tibble(y_tilde = y_tilde, pred_dist = pred_dist, Model = 'Mixture')
}

post_pred_mixture <- predictive_mixture(fit_mix_good)
plot(post_pred_mixture$y_tilde, log(post_pred_mixture$pred_dist))
sum(diff(post_pred_mixture$y_tilde)*post_pred_mixture$pred_dist[-1])

sum(diff(post_pred_mixture$y_tilde)*post_pred_mixture$y_tilde[-1]*post_pred_mixture$pred_dist[-1])
mean(fit_mixture$beta)



predictive_rlm <- function(fit_tukey, model){
  y_tilde <- seq(-50, 100, length.out = 100)
  joint <- as.data.frame(fit_tukey$jointPost)
  names(joint) <- c('mu', 'sigma2', 'posterior')
  dmu <- unique(diff(unique(joint$mu)))[1]
  dsigma2 <- unique(diff(unique(joint$sigma2)))[1]
  pred_dist <- sapply(y_tilde, function(yy){
  sum(dnorm(yy, joint$mu, sqrt(joint$sigma2))*joint$posterior)*dmu*dsigma2
  })
  tibble(y_tilde = y_tilde, pred_dist = pred_dist, Model = model)
}


post_pred_tukey <- predictive_rlm(fit_tukey, 'Tukey')
post_pred_huber <- predictive_rlm(fit_huber, 'Huber')
post_pred_lts <- predictive_rlm(fit_lts, 'LTS')
post_pred_lms <- predictive_rlm(fit_lms, 'LMS')


sum(post_pred_tukey$pred_dist[-1]*diff(post_pred_tukey$y_tilde))
sum(post_pred_huber$pred_dist[-1]*diff(post_pred_huber$y_tilde))
sum(post_pred_lts$pred_dist[-1]*diff(post_pred_lts$y_tilde))
sum(post_pred_lms$pred_dist[-1]*diff(post_pred_lms$y_tilde))


sum(post_pred_tukey$pred_dist[-1]*diff(post_pred_tukey$y_tilde)*post_pred_tukey$y_tilde[-1])
sum(post_pred_lms$pred_dist[-1]*diff(post_pred_lms$y_tilde)*post_pred_lms$y_tilde[-1])




predictive_distributions <- bind_rows(post_pred_normal, 
                                      post_pred_t,
                                      post_pred_mixture, 
                                      post_pred_tukey, 
                                      post_pred_huber, 
                                      post_pred_lms, 
                                      post_pred_lts,
                                      post_pred_restricted_t)
predictive_distributions$Model <- factor(predictive_distributions$Model, levels = model_order)

predictive_distributions <- predictive_distributions %>% 
  mutate(log_density = log(pred_dist))


#xlim(c(-10, 60))

# plot posterior predictive ----
ggplot(predictive_distributions, aes(x = y_tilde, y = log_density, col = Model)) + 
  geom_line() + ylim(c(-25,-2)) + theme_bw() +
  scale_color_brewer(palette = 'Set2') + xlab(expression(y)) + 
  ylab(expression(log('predictive density'))) + 
  scale_x_continuous(limits = c(-10,60), breaks = c(-10,20, 25, 30, 35,60))
ggsave(file.path(fig_path, 'speed_of_light_predictive.png'))



# other computations ----

# compute_post_mean <- function(fit){
#   sum(diff(fit$muPost[,1])*fit$muPost[-1,1]*fit$muPost[-1,2])
# }
# 
# rlm_tukey(y)
# compute_post_mean(fit_tukey)
# rlm_huber(y)
# compute_post_mean(fit_huber)
# rlm_lts(y)
# compute_post_mean(fit_lts)
# rlm_lms(y)
# compute_post_mean(fit_lms)
# mean(fit_normal[,1])
# mean(fit_t[,1])
# 
# 
# 
# 
# 
# 
# 
