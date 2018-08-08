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
  tmodelEstimates <- array(NA, c(p+1,nGroups,reps)) #thick tailed Bayes
  
  
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
  fits<-function(betahats, X){X%*%betahats}
  
 
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
        map2(.x = ., .y = models_lm, .f = function(x, y){
          predict(y, newdata = x)
        })
      
      olsMarginals <- list(yhold, olsPreds, models_lm) %>% 
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma) 
        })
          
        
   #Fit regressions models ------
      
    #unpooled regressions -----

    #rlm on training:Tukey -----
  
      state_rlm <- function(df){
        rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.bisquare, scale.est='Huber', data = df, maxit=1000)
      }
      
      models_rlm <- by_state$data %>% 
        map(state_rlm) 
      names(models_rlm) <- by_state$State
      
      
      rlmPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_rlm, .f = function(x, y){
          predict(y, newdata = x)
        })
      
      rlmMarginals <- list(yhold, rlmPreds, models_rlm) %>% 
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma) 
        })
      
      
    
      
    
    #rlm on training: Huber ----
      
      state_rlm <- function(df){
        rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.huber, scale.est='Huber', data = df, maxit=1000)
      }
      
      models_rlm <- by_state$data %>% 
        map(state_rlm) 
      names(models_rlm) <- by_state$State
      
      rlmPreds <- by_state_hold$data %>% 
        map2(.x = ., .y = models_rlm, .f = function(x, y){
          predict(y, newdata = x)
        })
      
      rlmMarginals_huber <- list(yhold, rlmPreds, models_rlm) %>% 
        pmap(.f = function(y, prd, m){
          dnorm(y,prd,summary(m)$sigma) 
        })
      


#Hierarchical Models -----      

##################################################    
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
      #betal
      betal <- array(NA, c(p,nkeep, nGroups))
      betal[p,,] <- t(nTheory$betal) #format expected for marginals computation
      postMeansBetal <- apply(betal,c(1,3) , mean)
      #postMeansBetal <- apply(nTheory$betal,1 , mean)
      diagnostic_betal <- geweke.diag(mcmc(t(nTheory$betal)))
      nTheorybetalConverge[] <- abs(diagnostic_betal$z)
      nTheoryGroupBetaMeans[,,i] <- postMeansBetal 
      postSDsBetal <- apply(nTheory$betal,c(1) , sd)
      nTheoryGroupBetaSDs[,,i]<-postSDsBetal
      #Beta
      postMeansBETA <- mean(nTheory$Beta) #colMeans(nTheory$Beta)
      nTheoryBetaMeans[,i] <- postMeansBETA
      postSDsBeta <- sd(nTheory$Beta) #apply(nTheory$Beta,2 , sd)
      nTheoryBetaSDs[,i] <- postSDsBeta
      nTheoryBETAconverge[,i] <- abs(geweke.diag(mcmc(nTheory$Beta))$z)
      #sigma2s
      postMeansSigma2s <- colMeans(nTheory$sigma2s)
      nTheorySigma2Means[,i] <- postMeansSigma2s
      postSDsSigma2s <- apply(nTheory$sigma2s,2,sd)
      nTheorySigma2SDs[,i] <- postSDsSigma2s
      #sigma2s converge?
      nTheorySigma2Converge[,i] <- abs(geweke.diag(mcmc(nTheory$sigma2s))$z)
      #bstar converge
      nTheorybstarConverge[i] <- abs(geweke.diag(mcmc(nTheory$bstar))$z)
      #mu_rho converge
      nTheoryMuRhoConverge[i] <- abs(geweke.diag(mcmc(nTheory$mu_rho))$z)
      #psi_rho_converge 
      nTheoryPsiRhoConverge[i] <- abs(geweke.diag(mcmc(nTheory$psi_rho))$z)
      #rho_converge
      nTheoryRhoConverge[i] <- abs(geweke.diag(mcmc(nTheory$rho))$z)
      
     
      #predictions on holdout set
      postMeansBetalList<-split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
      nTheoryPreds <- mapply(fits, postMeansBetalList,Xhold)
      nTheoryPredListofLists[[i]]<-nTheoryPreds
      
#computing marginal likelihoods for each element in holdout sample
nTheoryMarginalsListofLists[[i]]<-brlm::fn.compute.marginals.hierModelNormal(betal, nTheory$sigma2s, yhold,Xhold)
  
      
  
################################################
# restricted likelihood models -----

################################################
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

#betal
betal <- array(NA, c(p,nkeep, nGroups))
betal[p,,] <- t(restricted$betal) #get in format expected for marginals computation

diagnostic_betal <- geweke.diag(mcmc(t(restricted$betal)))
restrictedbetalConverge[] <- abs(diagnostic_betal$z)

postMeansBetal<-apply(betal,c(1,3) , mean)
restrictedGroupBetaMeans[,,i]<-postMeansBetal
postSDsBetal<-apply(betal,c(1,3) , sd)
restrictedGroupBetaSDs[,,i]<-postSDsBetal

#post means Beta
postMeansBETA<-mean(restricted$Beta) #colMeans(restricted$Beta)
restrictedBetaMeans[,i]<-postMeansBETA
postSDsBeta<-sd(restricted$Beta) #apply(restricted$Beta,2 , sd)
restrictedBetaSDs[,i]<-postSDsBeta

#Beta converge?
restrictedBETAconverge[,i]<-abs(geweke.diag(mcmc(restricted$Beta))$z)

#post means sigma2s
postMeansSigma2s<-colMeans(restricted$sigma2s)
restrictedSigma2Means[,i]<-postMeansSigma2s

#post sds sigma2s
postSDsSigma2s<-apply(restricted$sigma2s,2,sd)
restrictedSigma2SDs[,i]<-postSDsSigma2s
#sigma2s converge?
restrictedSigma2Converge[,i]<-abs(geweke.diag(mcmc(restricted$sigma2s))$z)

#bstar converge
restrictedbstarConverge[i]<-abs(geweke.diag(mcmc(restricted$bstar))$z)
#mu_rho converge
restrictedMuRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$mu_rho))$z)
#psi_rho_converge 
restrictedPsiRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$psi_rho))$z)
#rho_converge
restrictedRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$rho))$z)

#Acceptance rates
yAccept[i,]<-restricted$yAccept #just the column means

#predictionss on holdout set
postMeansBetalList <- split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
restrictedPreds <- mapply(fits, postMeansBetalList,Xhold)
restrictedPredListofLists[[i]]<-restrictedPreds

#computing marginal likelihoods for each element in holdout sample
restrictedMarginalsListofLists[[i]]<-fn.compute.marginals.hierModelNormal(betal, restricted$sigma2s, yhold,Xhold)

################################################

################################################     
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
#betal
betal <- array(NA, c(p,nkeep, nGroups))
betal[p,,] <- t(restricted_huber$betal) #get in format expected for marginals computation

diagnostic_betal <- geweke.diag(mcmc(t(restricted_huber$betal)))
restricted_huberbetalConverge[] <- abs(diagnostic_betal$z)

postMeansBetal<-apply(betal,c(1,3) , mean)
restricted_huberGroupBetaMeans[,,i]<-postMeansBetal
postSDsBetal<-apply(betal,c(1,3) , sd)
restricted_huberGroupBetaSDs[,,i]<-postSDsBetal

#post means Beta
postMeansBETA<-mean(restricted_huber$Beta) #colMeans(restricted_huber$Beta)
restricted_huberBetaMeans[,i]<-postMeansBETA
postSDsBeta<-sd(restricted_huber$Beta) #apply(restricted_huber$Beta,2 , sd)
restricted_huberBetaSDs[,i]<-postSDsBeta

#Beta converge?
restricted_huberBETAconverge[,i]<-abs(geweke.diag(mcmc(restricted_huber$Beta))$z)

#post means sigma2s
postMeansSigma2s<-colMeans(restricted_huber$sigma2s)
restricted_huberSigma2Means[,i]<-postMeansSigma2s

#post sds sigma2s
postSDsSigma2s<-apply(restricted_huber$sigma2s,2,sd)
restricted_huberSigma2SDs[,i]<-postSDsSigma2s
#sigma2s converge?
restricted_huberSigma2Converge[,i]<-abs(geweke.diag(mcmc(restricted_huber$sigma2s))$z)

#bstar converge
restricted_huberbstarConverge[i]<-abs(geweke.diag(mcmc(restricted_huber$bstar))$z)
#mu_rho converge
restricted_huberMuRhoConverge[i]<-abs(geweke.diag(mcmc(restricted_huber$mu_rho))$z)
#psi_rho_converge 
restricted_huberPsiRhoConverge[i]<-abs(geweke.diag(mcmc(restricted_huber$psi_rho))$z)
#rho_converge
restricted_huberRhoConverge[i]<-abs(geweke.diag(mcmc(restricted_huber$rho))$z)

#Acceptance rates
yAccept[i,]<-restricted_huber$yAccept #just the column means

#predictionss on holdout set
postMeansBetalList <- split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
restricted_huberPreds <- mapply(fits, postMeansBetalList,Xhold)
restricted_huberPredListofLists[[i]]<-restricted_huberPreds

#computing marginal likelihoods for each element in holdout sample
restricted_huberMarginalsListofLists[[i]]<-fn.compute.marginals.hierModelNormal(betal, restricted_huber$sigma2s, yhold,Xhold)  

      
################################################      
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
  
#betal
betal <- array(NA, c(p,nkeep, nGroups))
betal[p,,] <- t(tModel$betal) #get in format expected for marginals computation

diagnostic_betal <- geweke.diag(mcmc(t(tModel$betal)))
tModelbetalConverge[] <- abs(diagnostic_betal$z)

postMeansBetal<-apply(betal,c(1,3) , mean)
tModelGroupBetaMeans[,,i]<-postMeansBetal
postSDsBetal<-apply(betal,c(1,3) , sd)
tModelGroupBetaSDs[,,i]<-postSDsBetal

#post means Beta
postMeansBETA <- mean(tModel$Beta) #colMeans(tModel$Beta)
tModelBetaMeans[,i] <- postMeansBETA
postSDsBeta <- sd(tModel$Beta) #apply(tModel$Beta,2 , sd)
tModelBetaSDs[,i] <- postSDsBeta

#Beta converge?
tModelBETAconverge[,i]<-abs(geweke.diag(mcmc(tModel$Beta))$z)

#post means sigma2s
postMeansSigma2s <- colMeans(tModel$sigma2s)
tModelSigma2Means[,i]<-postMeansSigma2s

#post sds sigma2s
postSDsSigma2s<-apply(tModel$sigma2s,2,sd)
tModelSigma2SDs[,i]<-postSDsSigma2s
#sigma2s converge?
tModelSigma2Converge[,i]<-abs(geweke.diag(mcmc(tModel$sigma2s))$z)

#bstar converge
tModelbstarConverge[i]<-abs(geweke.diag(mcmc(tModel$bstar))$z)
#mu_rho converge
tModelMuRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$mu_rho))$z)
#psi_rho_converge 
tModelPsiRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$psi_rho))$z)
#rho_converge
tModelRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$rho))$z)

#preds on holdout set
postMeansBetalList <- split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
tModelPreds <- mapply(fits, postMeansBetalList,Xhold)
tModelPredListofLists[[i]]<-tModelPreds

#computing marginal likelihoods for each element in holdout sample
tModelMarginalsListofLists[[i]] <- fn.compute.marginals.hierModelTmodel(betal, tModel$sigma2s, yhold,Xhold)
################################################

    }
  )
  
###############################
  #fix this - what do I need to save for each run?
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
  
  write_rds(out, file.path(here::here(), paste0('hier_reg_n', n, '.rds' )))
  
}