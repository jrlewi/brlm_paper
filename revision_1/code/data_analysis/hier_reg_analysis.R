# Hierarchical regression - i.e., grouped by states. 

# 
# library(devtools)
# install_github('jrlewi/brlm')
library(sampling) #for stratified sampling
library(MCMCpack)
library(MASS)
library(brlm)
library(tidyverse)

#load data and prior information

analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))
parms_prior <- read_rds(file.path(here::here(), 'parms_prior.rds'))
analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>%
  group_by(State) %>% 
  filter(n() >= 10) %>% ungroup() %>% 
  mutate(State = factor(State)) %>% 
  arrange(State)

nGroups <- length(unique(analysis_data$State))

state_sizes <- analysis_data %>% 
  group_by(State) %>% 
  summarise(n = n())

#Set prior parameters ----
#place holders - needed to update
mu0 <- 1
Sigma0 <- .1
a0 <- 1
b0 <- 1
mu_bstr <- .5
psi_bstr <- 1
swSq <- 1
w1 <- 1
w2 <- 1
a_psir <- 1
b_psir <- 1


nu <- 3 #df for t-model

#n <- 1000  
ns <- c(1000) #, 2000) #sample size for training set
reps <- 2 # number of training sets

nburn <- 200 #set length of mcmc chains
nkeep <- 200
nkeept <- 400 #for the t-model.
maxit <- 400 #parameter in MASS::rlm

#set seed 
set.seed(123)
N <- nrow(analysis_data)
p <- length(mu_0)


for(n in ns){
  
  percent_sample <- n/N
  strata_sizes <- round(percent_sample*state_sizes$n, 0)
  
  
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
  # rlm regression estimates for classical robust regression
  # Posterior Means for the Bayesian regressions
  
  rlmEstimates <- array(NA, c(p+1,nGroups,reps)) #Tukey's Estimates
  rlmEstimatesHuber <- array(NA, c(p+1,nGroups,reps)) #Huber's
  olsEstimates <- array(NA, c(p+1,nGroups,reps)) #ols
  nTheoryEstimates <- array(NA, c(p+1,nGroups,reps)) #normal theory Bayes
  restrictedEstimates <- array(NA, c(p+1,nGroups,reps)) #Our method with Tukey's
  restrictedEstimatesHuber <- array(NA, c(p+1,nGroups,reps)) #Our method with Huber's
  tmodelEstimates <- array(NA, c(p+1,nGroups,reps))) #thick tailed Bayes
  
  
  #marginals of each y in holdout set  ------
  margRlm <- matrix(NA, nrow=reps, ncol=N-n)
  margRlmHuber <- matrix(NA, nrow=reps, ncol=N-n)
  margOls <- matrix(NA, nrow=reps, ncol=N-n)
  margNTheory <- matrix(NA, nrow=reps, ncol=N-n)
  margRest <- matrix(NA, nrow=reps, ncol=N-n)
  margRestHuber <- matrix(NA, nrow=reps, ncol=N-n)
  margT <- matrix(NA, nrow=reps, ncol=N-n)
  
  
  
  # Predictions of each y in holdout set  -----
  rlmPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  rlmPredHuberMat <- matrix(NA, nrow=reps, ncol=N-n)
  olsPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  nTheoryPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  restPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  restPredHuberMat <- matrix(NA, nrow=reps, ncol=N-n)
  tPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  
  
  # M-H acceptance rates -----
  acceptY <- numeric(reps) #acceptance rates for augmented y's in restricted model using Tukey's
  acceptYHuber <- numeric(reps) #acceptance rates for augmented y's in restricted model using Huber's
  acceptT <- numeric(reps) #acceptance rates for MH in the tmodel
  
  
  
  
  #auxilary functions and constants ---- 
  #t density with center and scale and df = nu
  tdensity<-function(y, mean, sigma){
    (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
  }
  
 
  # simulation -----
  system.time(  
    for(i in 1:reps){

    strat_sample <- sampling::strata(analysis_data, "State", strata_sizes, method = c('srswor'))
      
      
      #trainIndices <- sort(sample(1:N, n))
      holdIndices <- c(1:N)[-strat_sample$ID_unit]
      holdIndicesMatrix[i,] <- holdIndices
      train <- analysis_data[strat_sample$ID_unit,]
      hold <- analysis_data[holdIndices,]
      yholdout <- hold$sqrt_count_2012
      y_hold[i,] <- yholdout
     
 
      
      #type 1, open agencies to predict in holdout set
      open <- hold$Count_2012 > 0 
      type1 <- hold$Type == '1'
      y_open[i,] <-open
      y_type1[i,] <-type1
      
      
      # Set up regressions ----
      
      #prepare data for fitting function; get ols preds and marginals on holdoutset
      #y is list of responses from each group
      #X is list of design matrices for each group
  
        by_state <- train %>% 
        group_by(State) %>% 
        nest() 
      
      state_lm <- function(df){
        lm(sqrt_count_2012 ~ sqrt_count_2010 - 1, y = TRUE, x = TRUE, data = df)
      }
      
      models_lm <- by_state$data %>% 
        map(state_lm) 
      names(models_lm) <- by_state$State
      
      
      y <- models_lm %>% 
        map(.f = function(m) m$y)
      
      X <- models_lm %>% 
        map(.f = function(m) m$x)
      
      
      betaHats <- models_lm %>% 
        map(.f = function(m) coef(m)) %>% 
        unlist()
      
      sigHats <- models_lm %>% 
        map(.f = function(m) summary(m)$s) %>% 
        unlist()
      

      
      #prepare the holdout data for predictions
      #yhold is list of holdout responses from each group
      #Xhold is list of holfout design matrices for each group
      by_state_hold <- hold %>% 
        group_by(State) %>% 
        nest()
      
      models_lm_hold <- by_state_hold$data %>% 
      map(state_lm) 
      names(models_lm_hold) <- by_state_hold$State
      
      yhold <- models_lm_hold %>% 
        map(.f = function(m) m$y)
      
      Xhold <- models_lm_hold %>% 
        map(.f = function(m) m$x)
      
      olsPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_lm_hold, .f = function(x, y){
          predict(y, newdata = x)
        })
      
    
      olsMarginals <- list(yhold, olsPreds, models_lm) %>% 
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma) 
        })
          

      # olsPredListofLists[[i]][[group]]<-olsPreds
      # olsMarginalsListofLists[[i]][[group]]<-dnorm(yhold[[group]],olsPreds,summary(fit)$sigma)
        
   #Fit regressions models ------
      
    #unpooled regressions -----

    #rlm on training:Tukey -----
  
      state_rlm <- function(df){
        rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.bisquare, scale.est='Huber', data = df, maxit=1000)
      }
      
      models_rlm <- by_state$data %>% 
        map(state_rlm) 
      names(models_rlm) <- by_state_hold$State
      
      #fix this ----
      sigma2Int <- rlmfit$s^2 #to start mcmc
      rlmEstimates[,i] <- c(coef(rlmfit),rlmfit$s^2)
      rlmPreds <- Xholdout%*%coef(rlmfit)
      rlmPredMat[i,] <- rlmPreds
      margRlm[i,] <- dnorm(yholdout,mean=rlmPreds, sd=rlmfit$s)
      
    
    #rlm on training: Huber ----
      
      state_rlm_huber <- function(df){
        rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.huber, scale.est='Huber', data = df, maxit=1000)
      }
      
      models_rlm_huber <- by_state$data %>% 
        map(state_rlm_huber) 
      names(models_rlm_huber) <- by_state_hold$State
      
      #fix this ----- 
      rlmEstimatesHuber[,i]<-c(coef(rlmfitHuber),rlmfitHuber$s^2)
      rlmPredsHuber<-Xholdout%*%coef(rlmfitHuber)
      rlmPredHuberMat[i,] <- rlmPredsHuber
      margRlmHuber[i,] <- dnorm(yholdout,mean=rlmPredsHuber, sd=rlmfitHuber$s)
      
      
    #ols on training done above----
     models_lm 
      olsEstimates[,i] <- c(coef(olsFit),summary(olsFit)$sigma^2)
      olsPreds <- Xholdout%*%coef(olsFit)
      olsPredMat[i,] <- olsPreds
      margOls[i,] <- dnorm(yholdout,mean=olsPreds, sd=summary(olsFit)$sigma)
      

  #Hierarchical Models -----      
      
#normal theory bayes model ----
      
      #tunning parameters for MH step on bstar, mu_rho, psi_rho, and rho
      step_logbstar <- abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
      mu_rho_step <- .3 #(w1/(w1+w1))/sqrt(w1*w2/((w1+w2)^2*(w1+w1+1)))
      psi_rho_step <- a_psir^.5 #mean/sd
      rho_step <- .1
      
      nis <- unlist(lapply(y, length), use.names=FALSE)
      step_Z <- brlm:::fn.compute.Z(mean(sigHats^2), a0, b0)/(sqrt(nis))
      
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
    
   # fix this -----
      # #posterior means of beta
      # ##plot(nTheory$mcmc, ask=FALSE, density=FALSE)
      # postMeansNtheory <- colMeans(nTheory$mcmc)
      # nTheoryEstimates[,i] <- postMeansNtheory
      # betaNTheory <- postMeansNtheory[1:p]
      # sigma2Ntheory <- postMeansNtheory[p+1]
      # nTheoryPreds <- Xholdout%*%betaNTheory
      # nTheoryPredMat[i,] <- nTheoryPreds
      
      # get marginals f(y_h) for each element in holdout set 
      nMuMatrix <- Xholdout%*%(t(nTheory$mcmc)[1:p,])
      #sd's across samples for each houldout set
      nSigmaMat <- matrix(sqrt(rep(t(nTheory$mcmc)[p+1,],N-n)),N-n,nkeep, byrow = TRUE)
      #nSigmaMat2<-(sqrt(t(nTheory$mcmc)[p+1,]))
      #this is estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for normal model
      margNTheory[i,] <- rowMeans(dnorm(yholdout,mean = nMuMatrix, sd = nSigmaMat))
      
      

# restricted likelihood models -----
#Tukey version ----
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
                                        swSq=1,
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
      
      #fix this ----
      postMeansRest <- colMeans(restricted$mcmc)
      restrictedEstimates[,i] <- postMeansRest
      betaRest <- postMeansRest[1:p]
      sigma2Rest <- postMeansRest[p+1]
      restPreds <- Xholdout%*%betaRest
      restPredMat[i,] <- restPreds
      acceptY[i] <- mean(restricted$yAccept)
      
      
      # get marginals f(y_h) for each element in holdout set 
      #means across samples for each holdout set
      restMuMatrix <- Xholdout%*%(t(restricted$mcmc)[1:p,])
      #sd's across samples for each houldout set
      restSigmaMat <- matrix(sqrt(rep(t(restricted$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
      
      
      #estiamate (L(y_h)) of the marginal f(y_h) for each y_h in the holdout set for  restricted model
      margRest[i,] <- rowMeans(dnorm(yholdout,mean=restMuMatrix, sd=restSigmaMat))
      
      
#Huber version ----
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
                                               swSq=1,
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
      
      #fix this ----- 
      postMeansRestHuber <- colMeans(restrictedHuber$mcmc)
      restrictedEstimatesHuber[,i] <- postMeansRestHuber
      betaRestHuber <- postMeansRestHuber[1:p]
      sigma2RestHuber <- postMeansRestHuber[p+1]
      restPredsHuber <- Xholdout%*%betaRestHuber
      restPredHuberMat[i,] <- restPredsHuber
      acceptYHuber[i] <- mean(restrictedHuber$yAccept)
      
      #means across samples for each holdout set
      restMuMatrixHuber <- Xholdout%*%(t(restrictedHuber$mcmc)[1:p,])
      #sd's across samples for each houldout set
      restSigmaMatHuber <- matrix(sqrt(rep(t(restrictedHuber$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
      
      #estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for         restricted model
      margRestHuber[i,] <- rowMeans(dnorm(yholdout,mean=restMuMatrixHuber, sd=restSigmaMatHuber))
      
      
      
  # t-model ----
      #Note the change of prior on sigma2. In the other models the variance parameter is given  an IG(a_0,b_0) prior
      #the variance for the t model is (nu/(nu-2))sigma2~IG(a_0, b_0),  implies sigma2=(nu-2)/nu var(Y)~IG(a_0,(nu-2)/nu b_0)
      
      tModel <- brlm::hier_TLm(y,
                            X,
                            nkeep,
                            nburn,
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
      
      postMeansT <- colMeans(tmodel$mcmc)
      tmodelEstimates[,i] <- postMeansT
      betaT <- postMeansT[1:p]
      sigma2T <- postMeansT[p+1]
      TPreds <- Xholdout%*%betaT
      tPredMat[i,] <- TPreds
      acceptT[i] <- mean(tmodel$acceptSigma2)
      
      
      # estimate marginals.
      tMuMatrix <- Xholdout%*%(t(tmodel$mcmc)[1:p,])
      tSigmaMat <- matrix(sqrt(rep(t(tmodel$mcmc)[p+1,],N-n)),N-n,nkeept, byrow = TRUE)
      
      margT[i,] <- rowMeans(tdensity(yholdout, mean=tMuMatrix,sigma=tSigmaMat))
      print(i)
    }
  )
  

  
  out <- list(y_hold = y_hold,
              y_open = y_open,
              y_type1 = y_type1,
              acceptY = acceptY,
              acceptYHuber = acceptYHuber,
              acceptT = acceptT,
              margRlm = margRlm,
              margRlmHuber = margRlmHuber,
              margOls = margOls,
              margNTheory = margNTheory,
              margRest = margRest,
              margRestHuber = margRestHuber,
              margT = margT,
              holdIndicesMatrix = holdIndicesMatrix,
              rlmPredMat = rlmPredMat,
              rlmPredHuberMat = rlmPredHuberMat,
              olsPredMat = olsPredMat,
              nTheoryPredMat = nTheoryPredMat,
              restPredMat = restPredMat,
              restPredHuberMat = restPredHuberMat, 
              tPredMat = tPredMat,
              rlmEstimates = rlmEstimates,
              rlmEstimatesHuber = rlmEstimatesHuber,
              olsEstimates = olsEstimates,
              nTheoryEstimates = nTheoryEstimates,
              restrictedEstimates = restrictedEstimates,
              restrictedEstimatesHuber = restrictedEstimatesHuber,
              tmodelEstimates = tmodelEstimates)
  
  write_rds(out, file.path(here::here(), paste0('pooled_reg_n', n, '.rds' )))
  
}