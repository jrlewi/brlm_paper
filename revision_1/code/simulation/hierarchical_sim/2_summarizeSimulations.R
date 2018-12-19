#
# Summarize the simulations ---
#
library('tidyverse')
New_results <- FALSE


sims <- 1:30
ass.vec <- c(1.25, 5, 10) #c(1.25, 2.5,5,10,20) 
scale_vec <- c(0.5,1, 2) #round(c(0.5, 1/sqrt(2),1, sqrt(2), 2),2)
sig2 <- 4

if(New_results){

dfs_all <- vector('list', length(sims))
results_df <- file.path(getwd(), 'results')
for(data_sim in sims){
# load the data, along with the factor levels, and true values of theta
data <- read_rds(file.path(getwd(),'data_sig2_4', paste0('data_', data_sim, '.rds')))
n_groups <- length(data$theta)
factor_list <- as.data.frame(do.call(rbind, data$factorsList))  
names(factor_list) <- c('p', 'n','m')

# MSE evaluation metric -----

dfs <- vector('list', length(ass.vec)*length(scale_vec)*length(sig2))
i <- 1
for(ss in sig2){
  for(as in ass.vec){
    for(sc in scale_vec){
      rds_name <- paste0("_",data_sim, "__as_", as,"__scale_as_", sc, "__sig2_", ss, ".rds")
    
    huber <- read_rds(file.path(results_df, paste0("huber", rds_name)))
    
    #computes the estimates --- 
    theta_huber_rest <- apply(huber$theta, 2, mean)
    theta_huber_rlm  <- sapply(huber$robustFits, function(x) x$coefficients)
    sigma2_huber_rest <- apply(huber$sigma2, 2, mean)
    sigma2_huber_rlm  <- sapply(huber$robustFits, function(x) summary(x)$sigma^2)
    
    huber_df <- as_tibble(cbind(factor_list, theta = rep(data$theta,2), 
                       theta_hat = c(theta_huber_rest, theta_huber_rlm), 
                       sigma2_hat = c(sigma2_huber_rest,  sigma2_huber_rlm),
                       statistic = rep('Huber', 2*n_groups), 
                       method = c(rep('restricted', n_groups), rep('rlm', n_groups)),
                       a_s = rep(as, 2*n_groups),
                       scale_as = rep(sc, 2*n_groups),
                       sigma2 = rep(ss, 2*n_groups)
                       ))
    
    
    tukey <- read_rds(file.path(results_df, paste0("tukey", rds_name)))
    
    #computes the estimates --- 
    theta_tukey_rest <- apply(tukey$theta, 2, mean)
    theta_tukey_rlm  <- sapply(tukey$robustFits, function(x) x$coefficients)
    sigma2_tukey_rest <- apply(tukey$sigma2, 2, mean)
    sigma2_tukey_rlm  <- sapply(tukey$robustFits, function(x) summary(x)$sigma^2)
    
    tukey_df <- as_tibble(cbind(factor_list, theta = rep(data$theta,2), 
                       theta_hat = c(theta_tukey_rest, theta_tukey_rlm), 
                       sigma2_hat = c(sigma2_tukey_rest,  sigma2_tukey_rlm),
                       statistic = rep('Tukey', 2*n_groups), 
                       method = c(rep('restricted', n_groups), rep('rlm', n_groups)),
                       a_s = rep(as, 2*n_groups),
                       scale_as = rep(sc, 2*n_groups),
                       sigma2 = rep(ss, 2*n_groups)))
    
    
    
    
    
    normal <- read_rds(file.path(results_df, paste0("normal", rds_name)))
    #computes the estimates --- 
    theta_normal <- apply(normal$theta, 2, mean)
    sigma2_normal <- apply(normal$sigma2, 2, mean)
    
    normal_df <-  as_tibble(cbind(factor_list, theta = data$theta, 
                       theta_hat = c(theta_normal), 
                       sigma2_hat =  sigma2_normal, 
                       statistic = rep('Normal', n_groups), 
                       method =  rep('Normal', n_groups),
                       a_s = rep(as, n_groups),
                       scale_as = rep(sc, n_groups),
                       sigma2 = rep(ss, n_groups)
    ))
    
    df <- bind_rows(huber_df, tukey_df, normal_df)
    dfs[[i]] <- df
    i <- i+1
    }
  }
}

dfs_all[[data_sim]] <- bind_rows(dfs)

}

df_estimates <- bind_rows(dfs_all, .id = 'simulation')

} else {
summarize_data_frames <- readRDS(file.path(here::here(), "summarize_data_frames.rds")) 
 df_estimates <- summarize_data_frames[[2]]
}

#sims <- 1:length(unique(df_estimates$simulation))

df_mse <- df_estimates %>% dplyr::select(-p,-n,-m) %>%  
  group_by(simulation, statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(MSE = mean((theta - theta_hat)^2)) %>% 
  ungroup() %>% 
  group_by(statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(mean_MSE = mean(MSE), sd_MSE = sd(MSE)) %>% ungroup() %>% 
  mutate(a_s = as.character(a_s), scale_as = as.character(scale_as))

labels_vals <- scale_color_discrete(
                name="Method/Statistic",
                breaks= c("restricted.Huber",
                          "restricted.Tukey",
                          "rlm.Huber",
                          "rlm.Tukey"),
                labels=c("Restricted/Huber", 
                         "Restricted/Tukey", 
                         "Rlm/Huber", 
                         "Rlm/Tukey"))  


theme_set(theme_bw())
ggplot(df_mse %>%  filter(method != 'Normal' & statistic != 'Normal'), aes(x = as.factor(a_s), y = mean_MSE, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_MSE - sd_MSE/sqrt(length(sims)), ymax = mean_MSE + sd_MSE/sqrt(length(sims))), linetype = 1, width = .1, position = position_dodge(width = .5)) + geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~scale_as,  labeller = label_bquote(c == .(scale_as))) +
  labs(x = expression(a[s]), y = 'Average MSE') + labels_vals + theme(text = element_text(family = 'Times')) 
  
#  theme(text = element_text(family = 'Times')) + guides(col = guide_legend(title="Method/Statistic")) +
#+  guides(fill=guide_legend(title="Method/Statistic"))
ggsave(file.path(getwd(), "..", "..", "..", "figs", 'mse_sim2_facet_scale.png'))



ggplot(df_mse, aes(x = as.factor(scale_as), y = mean_MSE, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_MSE - sd_MSE/sqrt(length(sims)), ymax = mean_MSE + sd_MSE/sqrt(length(sims))), linetype = 1, width = .1, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~a_s,   labeller = label_bquote(a[s] == .(a_s))) +
  labs(x = expression(c),  y = 'Average MSE') + labels_vals + theme(text = element_text(family = 'Times')) 

ggsave(file.path(getwd(), "..", "..", "..", "figs", 'mse_sim2_facet_as.png'))



# K - L divergence evaluation metric -----
# K-L Divergence E[log f(y_good_i|theta_i, sig2True)/ f(y_good_i) ], expected value taken with respect to  f(y_good_i|theta_i, sig2True)


# function to compute the predictive distribution on a grid for one group


if(New_results){
compute_pred_dist<-function(ygrid, thetaSamples, sigma2Samples){
  #ygrid: grid of y values to evaluate the pred distribution
  #thetaSamples, sigma2Samples MCMC samples from the given group
  #OR for the robust regressions, these are just the estimates of theta and sigma2 from the robust regressions
  if(length(ygrid)>1000){
    sapply(ygrid, FUN=function(x) {
      mean(dnorm(x, thetaSamples, sqrt(sigma2Samples)))
    } 
    ) } else {
      ygridN <- matrix(ygrid, nrow = length(thetaSamples), ncol = length(ygrid), byrow = TRUE)
      colMeans(dnorm(ygridN, thetaSamples, sqrt(sigma2Samples)))
    }
}



# function to compute the K-L metric, for one group

compute_KL<-function(thetaTrue,thetaSamples, sigma2Samples, sig2True, ngridpts){
  #thetaTrue: scalar of the true value of theta for the given group
  #thetaSamples: MCMC sample of the theta_i from the given group
  #OR for robust regression it is just the estimate of theta for the given group
  #sigma2Sample: MCMC sample of the sigma2_i from the given group
  #OR for robust regression it is just the estimate of sigma2 for the given group
  #sig2True: scalar of the true value of sigma2
  #ygrid: grid of y values for the quadrature approximation of the expectation: not used anymore
  #ngridpts: number of grid points
  
  ygrid<-seq(thetaTrue-5*sqrt(sig2True),thetaTrue+5*sqrt(sig2True), length.out=ngridpts)
  predDist<-compute_pred_dist(ygrid,  thetaSamples, sigma2Samples)
  trueDist<-dnorm(ygrid, thetaTrue, sqrt(sig2True))
  dif <- diff(ygrid)
  kl <- ((log(trueDist)-log(predDist))*trueDist)[-1]
  KL <- sum((kl*dif)) #[kl!=Inf])
  KL
}


#exact: for use with the robust regression methods
KL_twoNormals<-function(mu1,sigma2_1, mu2, sigma2_2){
  #KL(p,q) where p~N(mu_1,sigma2_1),q~N(mu_2, sigma2_2)
  sig_1<-sqrt(sigma2_1)
  sig_2<-sqrt(sigma2_2)
  log(sig_2/sig_1)+(sigma2_1+(mu1-mu2)^2)/(2*sigma2_2)-1/2  
}



# Function to compute the KL for each group

compute_KL_each_group<-function(ngridpts, thetaTrueVect, sig2True,thetaSamplesMat, sigma2SamplesMat){
  # thetaTrueVect: the vector of true theta values for each group
  # thetaSampleMat: nTot by nGroups matrix of MCMC samples of the theta_i's
  #OR for robust regression a vector of estimates of the theta_is
  # sigma2SamplesMat: nTot by nGroups matrix of MCMC samples of the sigma2_i's
  #OR for robust regression a vector of estimates of the sigma2_is
  if(class(thetaSamplesMat)=='matrix'){
    if(class(sigma2SamplesMat)!='matrix'){ stop('error: thetaSamplesMat and sigma2SamplesMat must both be matrices of MCMC sampples or both be vectors of robust regression estimates')}  
    #create a list out of the columns of thetaSamplesMat and sigma2SamplesMat
    #this will be fed into mapply
    thetaSamplesMat<-split(thetaSamplesMat, col(thetaSamplesMat))
    sigma2SamplesMat<-split(sigma2SamplesMat, col(sigma2SamplesMat))
  } else {
    if(class(thetaSamplesMat)!='numeric' || class(sigma2SamplesMat)!='numeric'){
      stop('error: thetaSamplesMat and sigma2SamplesMat must both be matrices of MCMC sampples or both be vectors of robust regression estimates')
    }
  } 
  mapply(compute_KL,thetaTrue=thetaTrueVect,thetaSamples=thetaSamplesMat, sigma2Samples=sigma2SamplesMat, MoreArgs=list(sig2True=sig2True, ngridpts=ngridpts))
}

# Note: some numerical issues occur for the plug in predictive distribution for the robust regression estimates. The pred dist at the tail ygrid values is numerically zero. Could compute log pred for the robsust estimates. For unified code function fn.compute.pred to handle both MCMC samples and plug in, I handled this in compute.KL by adding KL<-sum((kl*dif)[kl!=Inf]). Both ways work here, but it is a note of caution. 

sig2 <- 4
ngridpts <- 100
# factor_list <- as.data.frame(do.call(rbind, data$factorsList))  
# names(factor_list) <- c('p', 'n','m')


dfs_all <- vector('list', length(sims))
for(data_sim in sims){
data <- read_rds(file.path(getwd(),'data_sig2_4', paste0('data_', data_sim, '.rds')))
factor_list <- as.data.frame(do.call(rbind, data$factorsList))  
names(factor_list) <- c('p', 'n','m')
dfs <- vector('list', length(ass.vec)*length(scale_vec)*length(sig2))
i <- 1
for(ss in sig2){
  for(as in ass.vec){
    for(sc in scale_vec){
      
      rds_name <- paste0("_",data_sim, "__as_", as,"__scale_as_", sc, "__sig2_", ss, ".rds")
      #huber
      fit <- read_rds(file.path(results_df, paste0("huber", rds_name)))
      #compute KL 
      kl <- compute_KL_each_group(ngridpts, thetaTrueVect = data$theta, sig2True = sig2,thetaSamplesMat = fit$theta, sigma2SamplesMat = fit$sigma2)
      huber_kl <- as_tibble(cbind(factor_list, KL = kl, statistic = 'Huber', method = 'restricted', a_s = as,
                                      scale_as = sc,
                                      sigma2 = ss))
      #KL for classical fit
      ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
      shat2<-sapply(fit$robustFits, FUN=function(x) x$s^2)
      kl <-  mapply(KL_twoNormals, mu1=data$theta,sigma2_1=sig2, mu2=ests,  sigma2_2=shat2)
     
    huber_kl_rlm <-  as_tibble(cbind(factor_list, KL = kl, statistic = 'Huber', method = 'rlm', a_s = as,
                         scale_as = sc,
                         sigma2 = ss))
      #tukey
      fit <- read_rds(file.path(results_df, paste0("tukey", rds_name)))
      #compute KL 
      kl <- compute_KL_each_group(ngridpts, thetaTrueVect = data$theta, sig2True = sig2,thetaSamplesMat = fit$theta, sigma2SamplesMat = fit$sigma2)
      tukey_kl <- as_tibble(cbind(factor_list, KL = kl, statistic = 'Tukey', method = 'restricted', a_s = as,
                                      scale_as = sc,
                                      sigma2 = ss))
      #KL for classical fit
      ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
      shat2<-sapply(fit$robustFits, FUN=function(x) x$s^2)
      kl <-  mapply(KL_twoNormals, mu1=data$theta,sigma2_1=sig2, mu2=ests,  sigma2_2=shat2)
      
      tukey_kl_rlm <-  as_tibble(cbind(factor_list, KL = kl, statistic = 'Tukey', method = 'rlm', a_s = as,
                                           scale_as = sc,
                                           sigma2 = ss))
      
      #normal
     fit <- read_rds(file.path(results_df, paste0("normal", rds_name)))
      #compute KL 
      kl <- compute_KL_each_group(ngridpts, thetaTrueVect = data$theta, sig2True = sig2,thetaSamplesMat = fit$theta, sigma2SamplesMat = fit$sigma2)
     normal_kl <- as_tibble(cbind(factor_list, KL = kl, statistic = 'Normal', method = 'Normal', a_s = as,
                                      scale_as = sc,
                                      sigma2 = ss))
      
      df <- bind_rows(huber_kl,huber_kl_rlm,tukey_kl,tukey_kl_rlm, normal_kl)
      dfs[[i]] <- df
      i <- i+1
    }
  }
}
dfs_all[[data_sim]] <- bind_rows(dfs)
print(data_sim)
}
df_kl <- bind_rows(dfs_all, .id = 'simulation')
saveRDS(list(df_kl, df_estimates), file = 'summarize_data_frames.rds') 
} else {
  summarize_data_frames <- readRDS(file.path(here::here(), "summarize_data_frames.rds")) 
  df_kl <- summarize_data_frames[[1]]
}



#Note - warnings are just coercion rules - variables treated as factors without the same levels coerced to characters...no big deal here. 

#fctrs <- c('p', 'n', 'm', 'statistic', 'method', 'a_s', 'scale_as')
#df_kl <- df_kl %>% mutate_at(fctrs, funs(factor))

df_kl_mean <- df_kl %>% 
  group_by(simulation, statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(KL = mean(KL)) %>% 
  ungroup() %>% 
  group_by(statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL))

ggplot(df_kl_mean %>% filter(method != 'Normal', statistic != 'Normal') , aes(x = as.factor(a_s), y = mean_KL, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_KL - sd_KL/sqrt(length(sims)), ymax = mean_KL + sd_KL/sqrt(length(sims))), linetype = 1, width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~scale_as,   labeller = label_bquote(c == .(scale_as))) +
  labs(x = expression(a[s]),  y = 'Average KL') + theme_bw() + labels_vals + theme(text = element_text(family = 'Times'))

ggsave(file.path(getwd(), "..", "..", "..", "figs", 'kl_sim2_facet_scale.png'), width = 6, height = 4)

# tmp <- df_kl %>% filter(scale_as == 2, a_s == 5, method == 'restricted', statistic == 'Huber')
# plot(tmp$KL)
#mean(tmp$KL[-c(1:1000)])

df_kl_mean %>% filter(method == 'Normal', statistic == 'Normal') %>% ungroup() %>% 
  summarise(min = min(mean_KL), max = max(mean_KL))

library(MCMCpack)
invg <- function(x){
  a <<- 10
  cc <<- 0.5
  b <<- 4*a*cc
  dinvgamma(x, a, b)  
}
curve(invg, 0, 20, main = paste0('a = ', a, '  ', 'c = ', cc))
abline(v = 4)

invg <- function(x, a, cc){
  # a <- a_cc[1]
  # cc <- a_cc[2]
  b <- 4*a*cc
  tibble(x = x, prior = dinvgamma(x, a, b), a = a, cc = cc)
}
x <- seq(.01, 20, by = .01)

a_cc <- expand.grid(ass.vec, scale_vec)

priors <- purrr::map2(.x = a_cc[,1], .y = a_cc[,2], .f = invg, x = x) %>% bind_rows() %>% 
  mutate(a = as.factor(a))

ggplot(priors, aes(x, prior, group = a, col = a)) + 
  geom_line() + 
  facet_wrap(~cc,labeller = label_bquote(c == .(cc))) +
  theme(text = element_text(family = 'Times')) + geom_vline(xintercept = 4, lty = 2, col = 'gray', lwd = .5) + 
  labs(x = expression(sigma^2),  y = 'density')
ggsave(file.path(getwd(), "..", "..", "..", "figs", 'priors_sigma2.png'), width = 6, height = 4) 
 

ggplot(df_kl_mean %>% filter(method != 'Normal', statistic != 'Normal'), aes(x = as.factor(scale_as), y = mean_KL, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_KL - sd_KL/sqrt(length(sims)), ymax = mean_KL + sd_KL/sqrt(length(sims))), linetype = 1, width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~a_s,   labeller = label_bquote(a[s] == .(a_s))) +
  labs(x = expression(c),  y = 'Average KL') +  theme_bw() + labels_vals + theme(text = element_text(family = 'Times'))

ggsave(file.path(getwd(), "..", "..", "..", "figs", 'kl_sim2_facet_as.png'), width = 6, height = 4)

# saveRDS(list(df_kl, df_estimates), file = 'summarize_data_frames.rds')
# investigate convergence of variance estimators in mixture model ---

# Effects of n, p, or m

#note the rlm results are simply repeated for each a_s and scale_as - when grouping by n (or p, or m) - the averages for each simulation are computed over the repeated results too - but these averages are the same as if one result was used. The SE are computed as the sd of these averages - so it doesn't change the calculation if I subset the rlm results to get rid of the repeats



#n ----
df_kl_mean_n <- df_kl %>% #filter(method != 'rlm') %>% 
  group_by(simulation, n, statistic, method, sigma2, a_s, scale_as) %>% 
  summarise(KL = mean(KL)) %>% 
  ungroup() %>% 
  group_by(n, statistic, method, sigma2, a_s, scale_as) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())



# m----
df_kl_mean_m <- df_kl %>% 
  group_by(simulation, m, statistic, method, sigma2, a_s, scale_as) %>% 
  summarise(KL = mean(KL)) %>% 
  ungroup() %>% 
  group_by(m, statistic, method, sigma2,  a_s, scale_as) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())



# p -----
df_kl_mean_p <- df_kl %>% 
  group_by(simulation, p, statistic, method, sigma2,  a_s, scale_as) %>% 
  summarise(KL = mean(KL)) %>% 
  ungroup() %>% 
  group_by(p, statistic, method, sigma2,  a_s, scale_as) %>% 
  summarise(mean_KL = mean(KL), sd_KL = sd(KL), num = n())




# combine for a facet wrap plot
df_nmp <- (bind_rows(df_kl_mean_n,df_kl_mean_m,df_kl_mean_p))
df_nmp <- (df_nmp %>% gather(variable, value, n,m,p, na.rm = TRUE))

df_mnp_sub <- df_nmp %>% filter(method != 'Normal', statistic != 'Normal', scale_as == 1, a_s == 5)
ggplot(df_mnp_sub , aes(x = as.factor(value), y = mean_KL, col = interaction(method,statistic), group = interaction(method,statistic))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = mean_KL - sd_KL/sqrt(num), ymax = mean_KL + sd_KL/sqrt(num)), linetype = 1, width = 0, position = position_dodge(width = .5)) +
  facet_wrap(~ variable,   labeller = label_bquote(.(variable)), scales = "free_x") +
  labs(x = "Value",  y = 'Average KL') + theme_bw() + labels_vals + theme(text = element_text(family = 'Times'))
ggsave(file.path(getwd(), "..", "..", "..", "figs", 'kl_sim2_mnp.png'), width = 6, height = 4)



# work in Steve's office ------


ave_kl_group <- df_kl %>% filter(a_s == 5, scale_as == 1) %>%
  group_by(simulation, statistic, method, p,m,n) %>% 
  summarise(mean_kl = mean(KL)) %>% 
  ungroup() %>%  
  group_by(statistic, method, p, m, n) %>% 
  summarise(mean = mean(mean_kl), sd = sd(mean_kl))
  
ave_kl_group %>% ungroup() %>% 
  group_by(statistic, method) %>% 
  summarize(mn = mean(mean))



ggplot(ave_kl_group %>% filter(!method == 'Normal'), aes(x = as.factor(n), y = mean, group = interaction(method, statistic), col = interaction(method, statistic))) +
  geom_point(position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = mean - sd/sqrt(length(sims)), ymax = mean + sd/sqrt(length(sims))), linetype = 1, width = 0, position = position_dodge(width = .5)) +
  facet_wrap(~ m + p)



df_kl %>% filter(a_s == 10, scale_as == 1) %>%
  group_by(simulation, statistic, method, p,m,n) %>% 
  summarise(mean_kl = mean(KL)) %>% 
  ungroup() %>%  
  filter(method == 'restricted', statistic == 'Huber', m == 9, p == 0.1)







kl_fun  <- function(x){
  .5*(x^2 - 1) - log(x)
}

curve(kl_fun, from = .1, to = 2)


#Asymptotics -----
library(MASS)
gen_data <- function(n,m,p){
  sigma2 <- 4
  ps <- rbinom(n,1,p)
  vars <- ifelse(ps, m, 1)*sigma2
  rnorm(n, 0, sqrt(vars))
}
x <- gen_data(n = 1000000 ,m = 25,p = .3)
fit <- rlm(x ~ 1, psi = psi.huber, scale.est = 'Huber')
fit$s^2

#investigate behavior of KL main effects a functions of m,n,p by looking at estimates

# est_mnp_rlm <- df_estimates %>% 
#   filter(method == 'rlm', a_s == '1.25', scale_as == '0.5') %>% #git rid of repeated rlm data.frames
#   gather(variable, value, p,n,m)   
# names(est_mnp_rlm)
# 
# View(head(filter(est_mnp_rlm, variable == 'p', value == "0.1")))
#   
# View(head(filter(est_mnp_rlm, variable == 'n', value == "25")))   
# View(df_estimates %>% 
#        filter(method == 'rlm', a_s == '1.25', scale_as == '0.5'))
# 
# ggplot(est_mnp_rlm, aes(x = as.factor(value), y = theta - theta_hat, col = statistic)) +
#   geom_boxplot() +
#   facet_wrap(~ variable,   labeller = label_bquote(.(variable)), scales = "free_x") + labs(y = bquote(theta - widehat(theta)), x = 'Value') + geom_hline(yintercept = 0, lty = 2) +
#   theme(text = element_text(family = 'Times'))
# ggsave(file.path(getwd(), "..", "..", "..", "figs", 'rlm_mnp_theta.png'), width = 6, height = 4)
# 
# ggplot(est_mnp_rlm, aes(x = as.factor(value), y = sqrt(sigma2_hat), col = statistic)) +
#   geom_boxplot() +
#   facet_wrap(~ variable,   labeller = label_bquote(.(variable)), scales = "free_x") + geom_hline(yintercept = 2, lty = 2) +
#   labs(y = bquote(widehat(sigma)), x = 'Value') +
#   theme(text = element_text(family = 'Times'))
# ggsave(file.path(getwd(), "..", "..", "..", "figs", 'rlm_mnp_sigma.png'), width = 6, height = 4)
# 
# est_mnp_tukey <- df_estimates %>% 
#   filter(method == 'restricted', statistic == 'Tukey') %>%
#   gather(variable, value, p,n,m)   
# 
# tmp_df <- df_estimates %>% 
#   filter(method == 'rlm',statistic == 'Tukey', a_s == '1.25', scale_as == '0.5') %>% dplyr::select(p,n,m, sigma2_hat) %>% mutate(sigma_hat = sqrt(sigma2_hat), p = as.factor(p), n = as.factor(n), m = as.factor(m)) %>% 
#   as_tibble()
# 
# fit_mnp <- lm(sigma_hat  ~ (m + n+ p)^2, data = tmp_df )
# summary(fit_mnp) 
# 
# 
# 
# tmp_df <- df_kl %>% 
#   filter(method == 'restricted', statistic == 'Tukey') %>% dplyr::select(p,n,m, KL) %>%  
#   mutate(p = as.factor(p), n = as.factor(n), m = as.factor(m)) 
# 
# fit_mnp <- lm(log(KL)  ~ (m + n+ p)^2 + m*n*p, data = tmp_df )
# summary(fit_mnp) 
# plot(fit_mnp)
# 
# ggplot(est_mnp_tukey, aes(x = as.factor(value), y = theta - theta_hat)) +
#   geom_boxplot() + geom_hline(yintercept = 0, lty = 2) +
#   facet_wrap(~ variable,   labeller = label_bquote(.(variable)), scales = "free_x") + labs(y = bquote(theta - widehat(theta)), x = 'Value') +
#   theme(text = element_text(family = 'Times'))
# 
# ggplot(est_mnp_tukey, aes(x = as.factor(value), y = sqrt(sigma2_hat))) +
#   geom_boxplot() +
#   facet_wrap(~ variable,   labeller = label_bquote(.(variable)), scales = "free_x") + geom_hline(yintercept = 2, lty = 2) +
#   labs(y = bquote(widehat(sigma)), x = 'Value') +
#   theme(text = element_text(family = 'Times'))


