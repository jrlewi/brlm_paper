#
# Summarize the simulations ---
#
library('tidyverse')

sims <- 1:10
dfs_all <- vector('list', length(sims))
for(data_sim in sims){
# load the data, along with the factor levels, and true values of theta
data <- read_rds(file.path(getwd(),'data_sig2_4', paste0('data_', data_sim, '.rds')))
n_groups <- length(data$theta)
results_df <- file.path(getwd(), 'results')

ass.vec <- c(1.25, 5, 10) #c(1.25, 2.5,5,10,20) 
scale_vec <- c(0.5,1, 2) #round(c(0.5, 1/sqrt(2),1, sqrt(2), 2),2)
sig2 <- 4


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
    huber_df <- tibble(theta = rep(data$theta,2), 
                       theta_hat = c(theta_huber_rest, theta_huber_rlm), 
                       statistic = rep('Huber', 2*n_groups), 
                       method = c(rep('restricted', n_groups), rep('rlm', n_groups)),
                       a_s = rep(as, 2*n_groups),
                       scale_as = rep(sc, 2*n_groups),
                       sigma2 = rep(ss, 2*n_groups)
                       )
    
    
    tukey <- read_rds(file.path(results_df, paste0("tukey", rds_name)))
    
    #computes the estimates --- 
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
    #computes the estimates --- 
    theta_normal <- apply(normal$theta, 2, mean)

    normal_df <- tibble(theta = data$theta, 
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

dfs_all[[data_sim]] <- bind_rows(dfs)

}

df_estimates <- bind_rows(dfs_all, .id = 'simulation')

df_mse <- df_estimates %>% 
  group_by(simulation, statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(MSE = mean((theta - theta_hat)^2)) %>% 
  ungroup() %>% 
  group_by(statistic, method, a_s, scale_as, sigma2) %>% 
  summarise(mean_MSE = mean(MSE), sd_MSE = sd(MSE)) %>% ungroup() %>% 
  mutate(a_s = as.character(a_s), scale_as = as.character(scale_as))

labels_vals <- scale_color_discrete(name="Method/Statistic",
                                   breaks= c("restricted.Huber","rlm.Huber","restricted.Tukey","rlm.Tukey"),
                                   labels=c("Restricted/Huber", "Rlm/Huber", "Restricted/Tukey", "Rlm/Tukey"))  

theme_set(theme_bw())
ggplot(df_mse %>%  filter(method != 'Normal' & statistic != 'Normal'), aes(x = as.factor(a_s), y = mean_MSE, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_MSE - sd_MSE/sqrt(length(sims)), ymax = mean_MSE + sd_MSE/sqrt(length(sims))), linetype = 2, width = .1, position = position_dodge(width = .5)) + geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~scale_as,  labeller = label_bquote(c == .(scale_as))) +
  labs(x = expression(a[s]), y = 'Average MSE') + labels_vals + theme(text = element_text(family = 'Times')) 
  
#  theme(text = element_text(family = 'Times')) + guides(col = guide_legend(title="Method/Statistic")) +
#+  guides(fill=guide_legend(title="Method/Statistic"))
ggsave(file.path(getwd(), "..", "..", "..", "figs", 'mse_sim2_facet_scale.png'))



ggplot(df_mse, aes(x = as.factor(scale_as), y = mean_MSE, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_MSE - sd_MSE/sqrt(length(sims)), ymax = mean_MSE + sd_MSE/sqrt(length(sims))), linetype = 2, width = .1, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~a_s,   labeller = label_bquote(a[s] == .(a_s))) +
  labs(x = expression(c),  y = 'Average MSE') + labels_vals + theme(text = element_text(family = 'Times')) 

ggsave(file.path(getwd(), "..", "..", "..", "figs", 'mse_sim2_facet_as.png'))



# K - L divergence evaluation metric -----
# K-L Divergence E[log f(y_good_i|theta_i, sig2True)/ f(y_good_i) ], expected value taken with respect to  f(y_good_i|theta_i, sig2True)


# function to compute the predictive distribution on a grid for one group

compute_pred_dist<-function(ygrid, thetaSamples, sigma2Samples){
  #ygrid: grid of y values to evaluate the pred distribution
  #thetaSamples, sigma2Samples MCMC samples from the given group
  #OR for the robust regressions, these are just the estimates of theta and sigma2 from the robust regressions
  if(length(ygrid)>1000){
    sapply(ygrid, FUN=function(x) {
      mean(dnorm(x, thetaSamples, sqrt(sigma2Samples)))
    } 
    ) } else {
      ygridN<-matrix(ygrid, nrow=length(thetaSamples), ncol=length(ygrid), byrow = TRUE)
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
  
  ygrid<-seq(thetaTrue-10*sqrt(sig2True),thetaTrue+10*sqrt(sig2True), length.out=ngridpts)
  predDist<-compute_pred_dist(ygrid,  thetaSamples, sigma2Samples)
  trueDist<-dnorm(ygrid, thetaTrue, sqrt(sig2True))
  dif<-diff(ygrid)
  kl<-((log(trueDist)-log(predDist))*trueDist)[-1]
  KL<-sum((kl*dif)[kl!=Inf])
  
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
factor_list <- as.data.frame(do.call(rbind, data$factorsList))  
names(factor_list) <- c('p', 'n','m')


dfs_all <- vector('list', length(sims))
for(data_sim in sims){
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


#Note - warnings are just coercion rules - variables treated as factors wihtout the same levels coerced to characters...no big deal here. 

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



ggplot(df_kl_mean %>% filter(method != 'Normal', statistic != 'Normal'), aes(x = as.factor(scale_as), y = mean_KL, col = interaction(method,statistic), group = interaction(method,statistic))) + geom_errorbar(aes(ymin = mean_KL - sd_KL/sqrt(length(sims)), ymax = mean_KL + sd_KL/sqrt(length(sims))), linetype = 1, width = 0, position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap(~a_s,   labeller = label_bquote(a[s] == .(a_s))) +
  labs(x = expression(c),  y = 'Average KL') +  theme_bw() + labels_vals + theme(text = element_text(family = 'Times'))

ggsave(file.path(getwd(), "..", "..", "..", "figs", 'kl_sim2_facet_as.png'), width = 6, height = 4)

saveRDS(list(df_kl, df_estimates), file = 'summarize_data_frames.rds')
# investigate convergence of variance estimators in mixture model ---



# gen_data <- function(m,n,p){
#   ps <- rbinom(n,1,p)
#   vars<-ifelse(ps, m, 1)*sig2
#   rnorm(n, 0, sqrt(vars))
# }
# m <- 25; n <- 100; p <- .3
# n_sims <- 100
# s2_ests <- sapply(1:n_sims, function(sim){
#   y <- gen_data(m,n,p)
#   fit <- MASS::rlm(y~1, scale.est = "Huber")
#   fit$s^2
# })
# 
# s2_ests <- tibble(var_estimates = s2_ests, m = rep(m, n_sims), n = rep(n, n_sims), p = rep(p, n_sims))
# 
# 
# ggplot(s2_ests, aes(x = var_estimates, col = factor(m), fill = factor(n), alpha = factor(p))) + geom_histogram() + geom_vline(xintercept = sig2)
# ggsave(file.path(getwd(), "..", "..", "..", "figs", 'hist_sig2_ests.png'))
# 
# 
# 
# # look at individual simualtions....
# tmp <- df_kl %>% 
#   group_by(simulation, statistic, method, a_s, scale_as, sigma2) %>% 
#   summarise(KL = mean(KL)) 
# 
# #View(tmp %>% ungroup() %>% group_by(simulation) %>%  
# #  summarise(min = min(KL)))
# #unique(tmp$simulation)
# 
# 
# tmp <- tmp %>% filter(simulation == '10')
# # ggplot(tmp , aes(x = as.factor(scale_as), y = KL, col = interaction(method,statistic), group = interaction(method,statistic))) +
# #   geom_point(position = position_dodge(width = .5)) + geom_line(position = position_dodge(width = .5)) +
# #   facet_wrap(~a_s) + theme_bw()
# 
# 
# ggplot(tmp , aes(x = as.factor(a_s), y = KL, col = interaction(method,statistic), group = interaction(method,statistic))) +
#   geom_point(position = position_dodge(width = .5)) + geom_line(position = position_dodge(width = .5)) +
#   facet_wrap(~scale_as) + theme_bw()





