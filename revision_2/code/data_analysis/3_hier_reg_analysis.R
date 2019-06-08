# Hierarchical regression - i.e., grouped by states. 

# download brlm package if needed
# library(devtools)
# install_github('jrlewi/brlm')
library(sampling) #for stratified sampling
library(MCMCpack)
library(MASS)
library(brlm)
library(tidyverse)

reps <- 10 #50 # number of training sets
sim_number <- 3 #for batch processing
nburn <- 1e4 #set length of mcmc chains
nkeep <- 1e4
nburnt <- 1.5e4
nkeept <- 1e4 #for the t-model.
maxit <- 1000 #parameter in MASS::rlm

# aux funcitons ----

fn.compute.marginals.hierModelNormal<-function(betalsamples, sigma2lsamples, yhold,Xhold){
  
  
  fits<-function(betahats, X){X%*%betahats}
  
  betalSampsList <- lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList) <- names(yhold)
  #holdout means
  
  nMuList<-mapply(fits, betalSampsList, Xhold)
  
  #sd's across samples for each holdout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)
  
  mapply(FUN=function(yh,mean,sd){
    pdfs <- dnorm(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE))
    mns <- rowMeans(pdfs)
    sds <- apply(pdfs, 1, sd)
    cbind(mns, sds)
  }
  ,yhold,nMuList,nSigmaList)
  
}

fn.compute.marginals.hierModelTmodel <- function(betalsamples, sigma2lsamples, yhold,Xhold){
  
  fits<-function(betahats, X){X%*%betahats}
  
  betalSampsList<-lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList)<-names(yhold)
  #holdout means
  nMuList<-mapply(fits, betalSampsList, Xhold)
  
  #sd's across samples for each houldout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)
  
  mapply(FUN=function(yh,mean,sd){
    # rowMeans(tdensity(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE), nu = nu))
    pdfs <- tdensity(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE), nu = nu)
    mns <- rowMeans(pdfs)
    sds <- apply(pdfs, 1, sd)
    cbind(mns, sds)
  }
  ,yhold,nMuList,nSigmaList)
  
}

#for classical regs
get_pseudo_samps <- function(inv_chi2_samps, models_lm){ 
  list(inv_chi2_samps, models_lm) %>% 
  pmap(.f = function(inv,m){
    s <- summary(m)$sigma^2
    b <- coef(m)
    ss <- s*inv
    cov_un <- as.numeric(summary(m)$cov.unscaled)
    sds <- sqrt(cov_un*ss)
    nn <- length(sds)
    bb <- rnorm(nn,mean = b, sd = sds)
    cbind(bb,ss)
  })}



# load data and prior information ----
analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))

# parms_prior <- read_rds(file.path(here::here(), 'parms_prior.rds'))
analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>%
  group_by(State) %>% 
  filter(n() >= 25) %>% ungroup() %>% 
  mutate(State = factor(State)) %>% #filter(!State %in% c(14, 30)) %>% 
  arrange(State)

analysis_data <- analysis_data %>% 
  filter(!State %in% c('14', '30')) %>% 
  mutate(State = droplevels(State))


#defined globally.
nGroups <<- length(unique(analysis_data$State))

state_sizes <- analysis_data %>% 
  group_by(State) %>% 
  summarise(n = n())

#Set prior parameters ----
parms_prior <- read_rds(file.path(here::here(), 'hier_parms_prior.rds'))
mu0 <<- parms_prior$mu_0
Sigma0 <<- parms_prior$Sigma_0
a0 <<- parms_prior$a_0
b0 <<- parms_prior$b_0
mu_bstr <<- parms_prior$mu_bstr
psi_bstr <<- parms_prior$psi_bstr
w1 <<- parms_prior$w1
w2 <<- parms_prior$w2
a_psir <<- parms_prior$a_psir
b_psir <<- parms_prior$b_psir


trend <- sqrt_count_2012 ~ sqrt_count_2010 - 1 + 
  Associate_Count +
  Office_Employee_Count

nu <<- 5 #df for t-model
ns <- floor(nrow(analysis_data)*.5)  #c(1000, 2000) #sample size for training set
# reps <- 10 #50 # number of training sets
# sim_number <- 1 #for batch processing
# nburn <- 1e4 #set length of mcmc chains
# nkeep <- 1e4
# nburnt <- 1.5e4
# nkeept <- 1e4 #for the t-model.
# maxit <- 1000 #parameter in MASS::rlm

# reps <- 50
# nburn <- 1e3 #set length of mcmc chains
# nkeep <- 1e3
# nburnt <- 1e3
# nkeept <- 1e3 #for the t-model.
# maxit <- 1000 #parameter in MASS::rlm


#set seed 
set.seed(123 + sim_number)
N <- nrow(analysis_data)
p <<- length(mu0)
nModels <- 7
# indices to id models in the arrays
ols_ind <- 1; 
rlm_ind <- 2; 
rlm_huber_ind <- 3
norm_ind <- 4
rest_ind <- 5
rest_hub_ind <- 6
t_ind <- 7
# norm_pt_ind <- 8
# rest_pt_ind <- 9
# rest_hub_pt_ind <- 10

for(n in ns){
  
  percent_sample <- n/N
  strata_sizes <- round(percent_sample*state_sizes$n)
 
  #make sure to sample exactly n values in total
  while(sum(strata_sizes) < n){
    ind <- which.max(strata_sizes)
    strata_sizes[ind] <- strata_sizes[ind] + 1
  }

  while(sum(strata_sizes) > n){
    ind <- which.max(strata_sizes)
    strata_sizes[ind] <- strata_sizes[ind] - 1
  }
  
  
  # Set storing objects -----
  
  # holdout samples ----
  # matrices to store 
  # 1 holdout samples
  # 2 logical: open type 1 agencies
  # 3 
  y_hold <- matrix(NA,  nrow=reps, ncol=N-n) 
  y_open <- y_type1 <- matrix(NA, nrow=reps, ncol=N-n) 
  holdIndicesMatrix <- matrix(NA, nrow=reps, ncol=N-n)
  
# estimates ----

#regression estimates or posterior means for each model.
#order is beta_ls, sigma
group_estimates <- group_estimates_sds <- array(NA, c(nModels, p+1,nGroups,reps)) 
#overall Beta
estimates <- estimates_sd <- array(NA, c(nModels, p, reps))
#marginals  and predictions of each y in holdout set  and each model ------
marginals <- predictions  <- array(NA, c(nModels , reps, N-n))
marginals_sd <- array(NA, c(nModels , reps, N-n))
# M-H acceptance rates -----
acceptY <- array(NA, c(2, nGroups, reps)) #acceptance rates for augmented y's in restricted models
# convergence diagnostics ----
group_converge <- array(NA, c(nModels, p+1, nGroups, reps))
#order: Beta, bstar, mu_rho, psi_rho, rho
converge <- array(NA, c(nModels, p + 4, reps))


#auxilary functions and constants ---- 
#t density with center and scale and df = nu
  # tdensity<-function(y, mean, sigma, nu){
  #   (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
  # }
  fits <- function(betahats, X){X%*%betahats}
  
 
# simulation -----
 strt <- Sys.time()
    for(i in 1:reps){
    resample <- TRUE
    
    while(resample){  
    strat_sample <- sampling::strata(analysis_data, "State", 
                                     strata_sizes, 
                                     method = c('srswor'))
      
      holdIndices <- c(1:N)[-strat_sample$ID_unit]
      holdIndicesMatrix[i,] <- holdIndices
      train <- analysis_data[strat_sample$ID_unit,]
      train$index <- strat_sample$ID_unit
      hold <- analysis_data[holdIndices,]
      hold$index <- holdIndices
      yholdout <- hold$sqrt_count_2012
      y_hold[i,] <- yholdout
  
      #make sure sample has at least 2 unique values
      #for additional variables so regressions can be fit
      chk <- hold %>% 
        group_by(State) %>% 
        summarise(sd_1  = sd(Associate_Count), 
                  sd_2 = sd(Office_Employee_Count)) %>% select(-State)
     any_0_sd <- map(chk,.f = function(s) any(s == 0)) %>% unlist() %>% 
       any()
     if(!any_0_sd){resample <- FALSE}
    }
      
      #type 1, open agencies to predict in holdout set
      open <- hold$Count_2012 > 0 
      type1 <- hold$Type == '1'
      y_open[i,] <-open
      y_type1[i,] <-type1
      
      
#  Fit classical regressions models -------
      
#prepare data for fitting function; get ols preds and marginals on holdoutset
      #y is list of responses from each group
      #X is list of design matrices for each group
  
# OLS on training -----      
        by_state <- train %>% 
        group_by(State) %>% 
        nest() 
      nis <- apply(by_state,1, function(dd) dd$data %>% nrow())
      state_lm <- function(df){
        lm(trend, y = TRUE, x = TRUE, data = df)
      }
      
      models_lm <- by_state$data %>% 
        map(state_lm) 
      names(models_lm) <- by_state$State
      
# convert response and design matrix for training set to lists for brlm functions
      y <- models_lm %>% 
        map(.f = function(m) m$y)
      X <- models_lm %>% 
        map(.f = function(m) m$x)
      
      betaHats <- models_lm %>% 
        map(.f = function(m) as.data.frame(coef(m))) %>% 
        bind_cols() %>% 
        as.matrix()
      
      group_estimates[ols_ind, 1:p, , i] <- betaHats
    
      sigHats <- models_lm %>% 
        map(.f = function(m) summary(m)$s) %>% 
        unlist()
      group_estimates[ols_ind, p+1, , i] <- sigHats^2

      #prepare the holdout data for predictions
      #yhold is list of holdout responses from each group
      #Xhold is list of holfout design matrices for each group
      by_state_hold <- hold %>% 
        group_by(State) %>% 
        nest()
      
      models_lm_hold <- by_state_hold$data %>% 
      map(state_lm) 
      names(models_lm_hold) <- by_state_hold$State

# convert response and design matrix for holdout set
      yhold <- models_lm_hold %>% 
        map(.f = function(m) m$y)
      
      Xhold <- models_lm_hold %>% 
        map(.f = function(m) m$x)
     
#predictions on the holdoutset
      olsPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_lm, .f = function(x, y){
        predict(y, newdata = x)
        })
      
  predictions[ols_ind, i, ] <- olsPreds %>% unlist()
  
  
  
#marginals on holdoutset  
  #get psuedo posterior samples - roughly same results
  # inv_chi2_samps <- lapply(nis, function(n) (n-p)/rchisq(nkeep, df = n-p))
  # pseudo_samps <- get_pseudo_samps(inv_chi2_samps, models_lm)
  # 
  # beta_samps <- sapply(pseudo_samps, function(x) x[,1])
  # betal <- array(NA, c(p,nkeep, nGroups))
  # betal[p,,] <- beta_samps
  # sigma2_samps <- sapply(pseudo_samps, function(x) x[,2])
  # 
  # ols_marg_mn_sd <-fn.compute.marginals.hierModelNormal(betal,  sigma2_samps, yhold,Xhold)
  # ols_marg <- lapply(ols_marg_mn_sd, function(x) x[,1])
  # ols_marg_sd <- lapply(ols_marg_mn_sd, function(x) x[,2])
  # marginals[ols_ind, i,] <- ols_marg %>% unlist() 
  # marginals_sd[ols_ind, i,] <- ols_marg_sd %>% unlist() 
  
olsMarginals <-
        list(yhold, olsPreds, models_lm) %>%
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma)
        })
marginals[ols_ind, i, ] <- olsMarginals %>% unlist()  

    
  #RLM Tukey on training -----
      state_rlm <- function(df){
        rlm(trend, psi=psi.bisquare, scale.est='Huber', 
            data = df, maxit=1000)
      }
      
      models_rlm <- by_state$data %>% 
        map(state_rlm) 
      names(models_rlm) <- by_state$State
      
      betaHats <- models_rlm %>% 
        map(.f = function(m) as.data.frame(coef(m))) %>% 
        bind_cols() %>% 
        as.matrix()
      group_estimates[rlm_ind, 1:p, , i] <- betaHats
      
      sigHats <- models_rlm %>% 
        map(.f = function(m) summary(m)$sigma) %>% 
        unlist()
      group_estimates[rlm_ind, p+1, , i] <- sigHats^2
      
      rlmPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_rlm, .f = function(x, y){
          predict(y, newdata = x)
        })
      
      predictions[rlm_ind, i, ] <- rlmPreds %>% unlist()
     
    # marginals ----
      rlmMarginals <- list(yhold, rlmPreds, models_rlm) %>%
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma)
        })
      
      marginals[rlm_ind, i, ] <- rlmMarginals %>% unlist()

      
      
      
      
      
    #rlm on training: Huber ----
      
      state_rlm <- function(df){
        rlm(trend, psi=psi.huber, 
            scale.est='Huber', data = df, maxit=1000)
      }
      
      models_rlm <- by_state$data %>% 
        map(state_rlm) 
      names(models_rlm) <- by_state$State
      
      betaHats <- models_rlm %>% 
        map(.f = function(m) as.data.frame(coef(m))) %>% 
        bind_cols() %>% 
        as.matrix()
      group_estimates[rlm_huber_ind, 1:p, , i] <- betaHats
      
      sigHats <- models_rlm %>% 
        map(.f = function(m) summary(m)$sigma) %>% 
        unlist()
      group_estimates[rlm_huber_ind, p+1, , i] <- sigHats^2
      
      
      rlmPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_rlm, .f = function(x, y){
          predict(y, newdata = x)
        })
      
      predictions[rlm_huber_ind, i, ] <- rlmPreds %>% unlist()  
      
        rlmMarginals <- list(yhold, rlmPreds, models_rlm) %>%
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma)
        })

        marginals[rlm_huber_ind, i, ] <- rlmMarginals %>% unlist()


# Hierarchical Models -----    

################################################
# normal theory bayes model ----
################################################
swSq <<- 1      

#tunning parameters for MH step on bstar, mu_rho, psi_rho, and rho
      step_logbstar <- abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
      mu_rho_step <- .3
      psi_rho_step <- (a_psir/b_psir^2)^.5 #mean/sd
      rho_step <- .1
      
      # nis <- unlist(lapply(y, length), use.names=FALSE)
      sigs2 <- group_estimates[ols_ind, p+1, , i]
      step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
      if(any(is.na(step_Z))){
        sigs2 <-  1
        step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
      }
      
      
      
      nTheory <- brlm::hierNormTheoryLm(y,
                                        X,
                                        nkeep,
                                        nburn,
                                        mu0,
                                        Sigma0,
                                        a0, 
                                        b0,
                                        mu_bstr,
                                        psi_bstr,
                                        swSq = 1,
                                        w1,
                                        w2, 
                                        a_psir,
                                        b_psir,
                                        step_logbstar, 
                                        mu_rho_step, 
                                        psi_rho_step, 
                                        rho_step,
                                        step_Z)
      #betal
      #betal <- array(NA, c(p,nkeep, nGroups))
      betal <- aperm(nTheory$betal, c(1,3,2)) #format expected for marginals computation
      postMeansBetal <- apply(betal,c(1,3) , mean)
      group_estimates[norm_ind, 1:p,,i]  <- postMeansBetal
      
      postSDsBetal <- apply(betal,c(1,3) , sd)
      group_estimates_sds[norm_ind, 1:p,,i]  <- postSDsBetal
      
    
      
      mcmc_format <- 
        lapply(seq(dim(betal)[3]), function(x) betal[ , , x]) 
      mcmc_format <- do.call(rbind, mcmc_format) %>% t() %>% 
        mcmc()
      group_converge[norm_ind, 1:p, ,i] <- 
        abs(geweke.diag(mcmc_format)$z)
      
      #sigma2s
      postMeansSigma2s <- colMeans(nTheory$sigma2s)
      group_estimates[norm_ind, p+1,,i] <- postMeansSigma2s
      
      postSDsSigma2s <- apply(nTheory$sigma2s,2,sd)
      group_estimates_sds[norm_ind, p+1,,i] <- postSDsSigma2s
      
      group_converge[norm_ind, p+1, ,i] <- 
        abs(geweke.diag(mcmc(nTheory$sigma2s))$z)
      
      
      #Beta
      postMeansBETA <- colMeans(nTheory$Beta)
      estimates[norm_ind, 1:p, i] <-  postMeansBETA
      
      postSDsBeta <- apply(nTheory$Beta,2 , sd)
      estimates_sd[norm_ind, 1:p, i] <- postSDsBeta  
      
      converge[norm_ind, 1:p, i] <- 
        abs(geweke.diag(mcmc(nTheory$Beta))$z)
      
     
      #bstar converge
      converge[norm_ind, p+1, i] <- 
        abs(geweke.diag(mcmc(nTheory$bstar))$z)
      #mu_rho converge
      converge[norm_ind, p+2, i] <- 
        abs(geweke.diag(mcmc(nTheory$mu_rho))$z)
      #psi_rho_converge 
      converge[norm_ind, p+3, i] <- 
        abs(geweke.diag(mcmc(nTheory$psi_rho))$z)
      #rho_converge
      converge[norm_ind, p + 4, i] <- 
        abs(geweke.diag(mcmc(nTheory$rho))$z)
      
#predictions on holdout set
      postMeansBetalList <- 
        split(postMeansBetal, rep(1:ncol(postMeansBetal), 
                                  each = nrow(postMeansBetal)))
      # sigs <- split(postMeansSigma2s^.5, rep(1:length(postMeansSigma2s), each = 1))
      nTheoryPreds <- mapply(fits, postMeansBetalList,Xhold)
      predictions[norm_ind, i, ] <- nTheoryPreds %>% unlist()
 
  
#computing marginal likelihoods for each element in holdout sample
nTheory_marg_mn_sd <- fn.compute.marginals.hierModelNormal(betal, 
                        nTheory$sigma2s, yhold,Xhold)
nTheory_marg <- lapply(nTheory_marg_mn_sd, function(x) x[,1])
marginals[norm_ind, i,] <- nTheory_marg %>% unlist() 
nTheory_marg_sd <- lapply(nTheory_marg_mn_sd, function(x) x[,2])
marginals_sd[norm_ind, i,] <- nTheory_marg_sd %>% unlist() 


################################################
# restricted likelihood models -----
################################################
#Tukey version ----

sigs2 <- group_estimates[rlm_ind, p+1, , i]
step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
if(any(is.na(step_Z))){
  sigs2 <-  1
  step_Z <-abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
}

      restricted <- brlm::hierNormTheoryRestLm(y,
                                        X,
                                        regEst = 'Tukey',
                                        scaleEst='Huber',
                                        nkeep, 
                                        nburn,
                                        mu0,
                                        Sigma0,
                                        a0, 
                                        b0,
                                        mu_bstr,
                                        psi_bstr,
                                        swSq = 1,
                                        w1,
                                        w2, 
                                        a_psir,
                                        b_psir,
                                        maxit=maxit,
                                        step_logbstar, 
                                        mu_rho_step, 
                                        psi_rho_step, 
                                        rho_step,
                                        step_Z)

#betal
#betal <- array(NA, c(p,nkeep, nGroups))
# betal[p,,] <- t(restricted$betal) #get in format expected for marginals computation

      betal <- aperm(restricted$betal, c(1,3,2)) #format expected for marginals computation
      
      
      mcmc_format <- 
        lapply(seq(dim(betal)[3]), function(x) betal[ , , x]) 
      mcmc_format <- do.call(rbind, mcmc_format) %>% t() %>% 
        mcmc()
      group_converge[rest_ind, 1:p, ,i] <- 
        abs(geweke.diag(mcmc_format)$z)
      
      
postMeansBetal <- apply(betal,c(1,3) , mean)
group_estimates[rest_ind, 1:p,,i] <-postMeansBetal
postSDsBetal<-apply(betal,c(1,3) , sd)
group_estimates_sds[rest_ind, 1:p,,i] <- postSDsBetal

#sigma2s
postMeansSigma2s <- colMeans(restricted$sigma2s)
group_estimates[rest_ind, p+1,,i] <- postMeansSigma2s

postSDsSigma2s<-apply(restricted$sigma2s,2,sd)
group_estimates_sds[rest_ind, p+1,,i] <- postSDsSigma2s

group_converge[rest_ind, p+1, ,i] <- 
  abs(geweke.diag(mcmc(restricted$sigma2s))$z)

#Beta
postMeansBETA <- colMeans(restricted$Beta)
estimates[rest_ind, 1:p, i] <-  postMeansBETA 

postSDsBeta <- apply(restricted$Beta,2 , sd)
estimates_sd[rest_ind, 1:p, i] <- postSDsBeta

converge[rest_ind, 1:p, i] <- abs(geweke.diag(mcmc(restricted$Beta))$z)

#bstar converge
converge[rest_ind, p + 1, i] <- abs(geweke.diag(mcmc(restricted$bstar))$z)
#mu_rho converge
converge[rest_ind, p + 2, i] <- abs(geweke.diag(mcmc(restricted$mu_rho))$z)
#psi_rho_converge 
converge[rest_ind, p + 3, i] <- abs(geweke.diag(mcmc(restricted$psi_rho))$z)
#rho_converge
converge[rest_ind, p + 4, i] <- abs(geweke.diag(mcmc(restricted$rho))$z)

#Acceptance rates for new y's
acceptY[1,,i] <- restricted$yAccept 

#predictionss on holdout set
postMeansBetalList <- split(postMeansBetal, 
                        rep(1:ncol(postMeansBetal), 
                            each = nrow(postMeansBetal)))
restrictedPreds <- mapply(fits, postMeansBetalList,Xhold)

predictions[rest_ind, i, ] <- restrictedPreds %>% unlist()

#computing marginal likelihoods for each element in holdout sample
rest_marg_mn_sd <- 
  fn.compute.marginals.hierModelNormal(betal, restricted$sigma2s, yhold,Xhold)
rest_marg <- lapply(rest_marg_mn_sd, function(x) x[,1])
marginals[rest_ind, i,] <- rest_marg %>% unlist() 
rest_marg_sd <- lapply(rest_marg_mn_sd, function(x) x[,2])
marginals_sd[rest_ind, i,] <- rest_marg_sd %>% unlist() 



# rest_pt_marg <- list(yhold,  restrictedPreds, postMeansSigma2s^.5) %>%
#   pmap(.f = function(y, prd, sig){
#     dnorm(y,prd,sig) 
#   })  
# marginals[rest_pt_ind, i,] <- rest_pt_marg %>% unlist() 
# 

################################################

################################################     
#Huber version ----

sigs2 <- group_estimates[rlm_huber_ind, p+1, , i]
step_Z <-abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
if(any(is.na(step_Z))){
  sigs2 <-  1
  step_Z <-abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
}

      restricted_huber <- brlm::hierNormTheoryRestLm(y,
                                               X,
                                               regEst = 'Huber',
                                               scaleEst='Huber',
                                               nkeep, 
                                               nburn,
                                               mu0,
                                               Sigma0,
                                               a0, 
                                               b0,
                                               mu_bstr,
                                               psi_bstr,
                                               swSq = 1,
                                               w1,
                                               w2, 
                                               a_psir,
                                               b_psir,
                                               maxit=maxit,
                                               step_logbstar, 
                                               mu_rho_step, 
                                               psi_rho_step, 
                                               rho_step,
                                               step_Z)
#betal
      betal <- aperm(restricted_huber$betal, c(1,3,2)) #format expected for marginals computation
      
      
      mcmc_format <- 
        lapply(seq(dim(betal)[3]), function(x) betal[ , , x]) 
      mcmc_format <- do.call(rbind, mcmc_format) %>% t() %>% 
        mcmc()
      group_converge[rest_hub_ind, 1:p, ,i] <- 
        abs(geweke.diag(mcmc_format)$z)
      
postMeansBetal <- apply(betal,c(1,3) , mean)
group_estimates[rest_hub_ind, 1:p,,i] <-postMeansBetal
postSDsBetal<-apply(betal,c(1,3) , sd)
group_estimates_sds[rest_hub_ind, 1:p,,i] <- postSDsBetal

#sigma2s
postMeansSigma2s <- colMeans(restricted_huber$sigma2s)
group_estimates[rest_hub_ind, p+1,,i] <- postMeansSigma2s

postSDsSigma2s<-apply(restricted_huber$sigma2s,2,sd)
group_estimates_sds[rest_hub_ind, p+1,,i] <- postSDsSigma2s

group_converge[rest_hub_ind, p+1, ,i] <- 
  abs(geweke.diag(mcmc(restricted_huber$sigma2s))$z)

#Beta
postMeansBETA <- colMeans(restricted_huber$Beta)
estimates[rest_hub_ind, 1:p, i] <-  postMeansBETA 

postSDsBeta <- apply(restricted_huber$Beta,2 , sd)
estimates_sd[rest_hub_ind, 1:p, i] <- postSDsBeta

converge[rest_hub_ind, 1:p, i] <- abs(geweke.diag(mcmc(restricted_huber$Beta))$z)

#bstar converge
converge[rest_hub_ind, p + 1, i] <- abs(geweke.diag(mcmc(restricted_huber$bstar))$z)
#mu_rho converge
converge[rest_hub_ind, p + 2, i] <- abs(geweke.diag(mcmc(restricted_huber$mu_rho))$z)
#psi_rho_converge 
converge[rest_hub_ind, p + 3, i] <- abs(geweke.diag(mcmc(restricted_huber$psi_rho))$z)
#rho_converge
converge[rest_hub_ind, p + 4, i] <- abs(geweke.diag(mcmc(restricted_huber$rho))$z)

#Acceptance rates for new y's
acceptY[2,,i] <- restricted_huber$yAccept 

#predictionss on holdout set
postMeansBetalList <- split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
restricted_huberPreds <- mapply(fits, postMeansBetalList,Xhold)

predictions[rest_hub_ind, i, ] <- restricted_huberPreds %>% unlist()

#computing marginal likelihoods for each element in holdout sample
rest_marg_mn_sd <- fn.compute.marginals.hierModelNormal(betal, restricted_huber$sigma2s, yhold,Xhold)
rest_marg <- lapply(rest_marg_mn_sd, function(x) x[,1])
marginals[rest_hub_ind, i,] <- rest_marg %>% unlist() 
rest_marg_sd <- lapply(rest_marg_mn_sd, function(x) x[,2])
marginals_sd[rest_hub_ind, i,] <- rest_marg_sd %>% unlist() 

# rest_hub_pt_marg <- list(yhold,  restricted_huberPreds, postMeansSigma2s^.5) %>%
#   pmap(.f = function(y, prd, sig){
#     dnorm(y,prd,sig) 
#   })  
# marginals[rest_hub_pt_ind, i,] <- rest_hub_pt_marg  %>% unlist() 

################################################      
# t-model ----
#Note the change of prior on sigma2. In the other models the variance parameter is given  an IG(a_0,b_0) prior
#the variance for the t model is (nu/(nu-2))sigma2~IG(a_0, b_0),  implies sigma2=(nu-2)/nu var(Y)~IG(a_0,(nu-2)/nu b_0)
    
      tModel <- brlm::hier_TLm(y,
                            X,
                            nkeept,
                            nburnt,
                            mu0,
                            Sigma0,
                            a0, 
                            ((nu-2)/nu)*b0,
                            mu_bstr,
                            psi_bstr,
                            swSq = 1,
                            w1,
                            w2, 
                            a_psir,
                            b_psir,
                            nu,
                            step_logbstar,
                            mu_rho_step,
                            psi_rho_step,
                            rho_step,
                            step_Z)
  
#betal
betal <- aperm( tModel$betal, c(1,3,2)) #format expected for marginals computation

  mcmc_format <- 
    lapply(seq(dim(betal)[3]), function(x) betal[ , , x]) 
  mcmc_format <- do.call(rbind, mcmc_format) %>% t() %>% 
      mcmc()
  group_converge[t_ind, 1:p, ,i] <- 
  abs(geweke.diag(mcmc_format)$z)

postMeansBetal <- apply(betal,c(1,3) , mean)
group_estimates[t_ind, 1:p,,i] <- postMeansBetal
postSDsBetal <- apply(betal,c(1,3) , sd)
group_estimates_sds[t_ind, 1:p,,i] <- postSDsBetal



#sigma2s
postMeansSigma2s <- colMeans(tModel$sigma2s)
group_estimates[t_ind, p+1,,i] <- postMeansSigma2s

postSDsSigma2s<-apply(tModel$sigma2s,2,sd)
group_estimates_sds[t_ind, p+1,,i] <- postSDsSigma2s

group_converge[t_ind, p+1, ,i] <- abs(geweke.diag(mcmc(tModel$sigma2s))$z)



#Beta
postMeansBETA <- colMeans(tModel$Beta)
estimates[t_ind, 1:p, i] <-  postMeansBETA 

postSDsBeta <- apply(tModel$Beta,2 , sd)
estimates_sd[t_ind, 1:p, i] <- postSDsBeta

converge[t_ind, 1:p, i] <- abs(geweke.diag(mcmc(tModel$Beta))$z)

#bstar converge
converge[t_ind, p + 1, i] <- abs(geweke.diag(mcmc(tModel$bstar))$z)
#mu_rho converge
converge[t_ind, p + 2, i] <- abs(geweke.diag(mcmc(tModel$mu_rho))$z)
#psi_rho_converge 
converge[t_ind, p + 3, i] <- abs(geweke.diag(mcmc(tModel$psi_rho))$z)
#rho_converge
converge[t_ind, p + 4, i] <- abs(geweke.diag(mcmc(tModel$rho))$z)

#preds on holdout set
postMeansBetalList <- split(postMeansBetal, rep(1:ncol(postMeansBetal), 
                                                each = nrow(postMeansBetal)))
tModelPreds <- mapply(fits, postMeansBetalList,Xhold)
predictions[t_ind, i, ]  <- tModelPreds %>% unlist()

#computing marginal likelihoods for each element in holdout sample
t_marg_mn_sd <- fn.compute.marginals.hierModelTmodel(betal, tModel$sigma2s, yhold,Xhold)
t_marg <- lapply(t_marg_mn_sd, function(x) x[,1])
marginals[t_ind, i,] <- t_marg %>% unlist() 
t_marg_sd <- lapply(t_marg_mn_sd, function(x) x[,2])
marginals_sd[t_ind, i,] <- t_marg_sd %>% unlist() 


################################################
print(i)
end <- Sys.time() - strt
print(end)
    }
end <- Sys.time() - strt
print(end)  
###############################
  #save output for each n
  out <- list(y_hold = y_hold,
              y_open = y_open,
              y_type1 = y_type1,
              holdIndices = holdIndicesMatrix,
              group_estimates = group_estimates,
              group_estimates_sds = group_estimates_sds,
              estimates = estimates,
              estimates_sd = estimates_sd, 
              marginals = marginals,
              marginals_sd = marginals_sd,
              predictions =  predictions,
              acceptY = acceptY,
              group_converge = group_converge,
              converge = converge
              )
  
  write_rds(out, file.path(here::here(), paste0('hier_reg_n', n, '_sim_number_', sim_number, '.rds' )))
  
}
