# checking calculation of kl_divergence for the classical models.

library(tidyverse)
library(MASS)


fn.gen.data <- function(sigma2) {
  p <- c(.1,.2,.3)
  n <- c(25,50, 100)
  m <- c(9, 25)
  factors <- expand.grid(p,n,m)
  yList <- list()
  factorsList <- list()
  theta <- numeric(nrow(factors)*5)
  for (nfac in 1:nrow(factors)) {
    p <- factors[nfac,1] 
    n <- factors[nfac,2]
    m <- factors[nfac,3]
    for (j in 1:5) { 
      thetai <- rnorm(1) #muTrue=0, tau2True=1
      theta[5*(nfac - 1) + j] <- thetai
      ps <- rbinom(n,1,p)
      vars <- ifelse(ps, m, 1)*sigma2
      yList[[5*(nfac - 1) + j]] <- rnorm(n, thetai, sqrt(vars))
      factorsList[[5*(nfac - 1) + j]] <- c(p,n,m)
    }
  }
  out <- list()
  out$factorsList <- factorsList
  out$yList <- yList
  out$theta <- theta
  out
}


#exact: for use with the robust regression methods
KL_twoNormals <- function(mu1,sigma2_1, mu2, sigma2_2){
  #KL(p,q) where p~N(mu_1,sigma2_1),q~N(mu_2, sigma2_2)
  sig_1 <- sqrt(sigma2_1)
  sig_2 <- sqrt(sigma2_2)
  log(sig_2/sig_1) + (sigma2_1 + (mu1 - mu2)^2)/(2*sigma2_2)-1/2  
}




sig2True <- 4 
n_sim <- 100
results_list <- lapply(1:n_sim, function(i){
YList <- fn.gen.data(sig2True)
factorsList <- YList$factorsList
factorsMat <- as.data.frame(matrix(unlist(factorsList), nrow = 90, ncol = 3, byrow = TRUE))
names(factorsMat) <- c("p", "n", "m")
thetaTrue <- YList$theta
factorsMat$theta <- thetaTrue
YList <- YList$yList
ests <- sapply(YList, function(y){
  fit <- rlm(y ~ 1,  psi = psi.huber, scale.est = 'Huber')
  c(coef(fit), fit$s)
})
factorsMat$theta_hat <- ests[1,]
factorsMat$sigma2_hat <- ests[2,]^2
kl <- KL_twoNormals(mu1 = thetaTrue, sigma2_1 = 4, mu2 = factorsMat$theta_hat, sigma2_2 = factorsMat$sigma2_hat)
factorsMat$kl <- kl
factorsMat
# tst2 <- numeric(90)
# for(j in 1:90){
#   tst2[j] <- KL_twoNormals(mu1 = thetaTrue[j], sigma2_1 = 4, mu2 = factorsMat$theta_hat[j], sigma2_2 = factorsMat$sigma2_hat[j])
# }
})


results <- bind_rows(results_list, .id = 'simulation') %>% as_tibble()


ave_results <- results %>% 
  group_by(simulation, m, n, p) %>% 
  summarize(mean_kl = mean(kl)) %>% 
  ungroup() %>% 
  group_by(m, n, p) %>% 
  summarize(mean = mean(mean_kl), sd = sd(mean_kl), num = n())

theme_set(theme_bw())
ggplot(ave_results, aes(x = as.factor(n), y = mean)) +
  geom_point() + geom_line(mapping = aes(group = 1)) + geom_errorbar(aes(ymin = mean- sd/sqrt(num), ymax = mean + sd/sqrt(num)), linetype = 1,  width = 0) +
  facet_wrap(~m+p)
ggsave('six_panel_check.png')


df_kl_mean_n <- results %>% 
  group_by(simulation, n) %>% 
  summarise(KL = mean(kl)) %>% 
  ungroup() %>% 
  group_by(n) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())

df_kl_mean_m <- results %>% 
  group_by(simulation, m) %>% 
  summarise(KL = mean(kl)) %>% 
  ungroup() %>% 
  group_by(m) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())
df_kl_mean_p <- results %>% 
  group_by(simulation, p) %>% 
  summarise(KL = mean(kl)) %>% 
  ungroup() %>% 
  group_by(p) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())




# using facet wrap
df_nmp <- (bind_rows(df_kl_mean_n,df_kl_mean_m,df_kl_mean_p))
df_nmp <- (df_nmp %>% gather(variable, value, n,m,p, na.rm = TRUE))


ggplot(df_nmp , aes(x = as.factor(value), y = mean_KL, group = 1)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = mean_KL - sd_KL/sqrt(num), ymax = mean_KL + sd_KL/sqrt(num)), linetype = 1,  width = 0) +
  facet_wrap(~ variable,  labeller = label_bquote(.(variable)), scales = "free_x")
ggsave('mnp_check.png')




