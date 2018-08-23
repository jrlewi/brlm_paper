# Regressions within a single state 

# download brlm package if needed
# library(devtools)
# install_github('brlm', 'jrlewi')
library(MCMCpack)
library(MASS)
library(brlm)
library(tidyverse)


strt <- Sys.time()

states <- c(2, 15, 27, 36)
ns <- c(25, 50) 
v_inflate <- c(100)

sims <-1:50
reps <-  length(sims)# number of training sets

nburn <- 1e4 #set length of mcmc chains
nkeep <- 1e4
nkeept <- 1e4
nburnt <- 2e4 #for the t-model.
nu <- 5
#set seed 
set.seed(min(sims))

for(State_keep in states){
  
#prior -----
prior_data <- read_rds(file.path(here::here(), 'data', 'prior_data.rds'))
# cnts <- prior_data %>%  group_by(State) %>% summarize(n = n())
# View(cnts)
prior_data <- prior_data %>% 
  mutate(sqrt_count_2008 = sqrt(Count_2008), sqrt_count_2010 = sqrt(Count_2010))  %>% filter(State == State_keep, Type == 1, Count_2010 > 0)

# ggplot(prior_data) + geom_point(aes(x = sqrt_count_2008, y = sqrt_count_2010)) + xlim(c(0,100)) +ylim(c(0,100))

# pooled regression analysis ----
prior_fit <- MASS::rlm(sqrt_count_2010 ~ sqrt_count_2008 - 1 , scale.est = 'Huber', data =  prior_data, maxit = 100)

# theme_set(theme_bw(base_family = 'Times'))
# ggplot(prior_data, aes(x = sqrt_count_2008, y = sqrt_count_2010)) + geom_point(size = 1) + xlim(c(0,150)) +  ylim(c(0, 130)) + stat_smooth(method = "rlm", formula = y ~ x - 1, method.args = list(scale.est = 'Huber', maxit = 100), size = .5, lty = 2, se = FALSE, col = 1) + guides(col = guide_legend(title = 'Agency Type'))

#load analysis data -------

analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))
# cnts <- analysis_data %>%  group_by(State) %>% summarize(n = n()) 
# 
# View(cnts)

analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>% filter(State == State_keep)

ggplot(analysis_data) + geom_point(aes(x = sqrt_count_2010, y = sqrt_count_2012)) + xlim(c(0,100))

#MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1 , scale.est = 'Huber', data =  analysis_data, maxit = 100)

N <- nrow(analysis_data)
p <- length(coef(prior_fit))

for(n_percent in ns){ # percent to use as training set. 
  
  converge_matrix <- array(NA, c(4, 2, length(v_inflate), length(sims))) 
  
  n <- floor(n_percent*N/100)
  # Set storing objects -----
  
  # holdout samples ----
  # matrices to store 
  # 1 holdout samples
  # 2 logical: open type 1 agencies
  # 3 
  y_hold <- matrix(NA,  nrow=reps, ncol=N-n) #mse's for predictions on entire holdoutset
  y_open <- y_type1 <- matrix(NA, nrow=reps, ncol=N-n) #on just the open type 1 agencies
  holdIndicesMatrix<-matrix(NA, nrow=reps, ncol=N-n)
  
  
  
  # estimates ----
  # rlm regression estimates for classical robust regression
  # Posterior Means for the Bayesian regressions
  
  rlmEstimates <- matrix(NA, nrow=p+1, ncol=reps) #Tukey's Estimates
  rlmEstimatesHuber <- matrix(NA, nrow=p+1, ncol=reps) #Huber's
  olsEstimates <- matrix(NA, nrow=p+1, ncol=reps) #ols
  nTheoryEstimates <- array(NA, c(p+1, reps, length(v_inflate))) #normal theory Bayes
  restrictedEstimates <-  array(NA, c(p+1, reps, length(v_inflate))) #Our method with Tukey's
  restrictedEstimatesHuber <-  array(NA, c(p+1, reps, length(v_inflate))) #Our method with Huber's
  tmodelEstimates <- array(NA, c(p+1, reps, length(v_inflate))) #thick tailed Bayes
  
  
  #marginals of each y in holdout set  ------
  margRlm <- matrix(NA, nrow=reps, ncol=N-n)
  margRlmHuber <- matrix(NA, nrow=reps, ncol=N-n)
  margOls <- matrix(NA, nrow=reps, ncol=N-n)
  margNTheory <- array(NA,c(reps,length(v_inflate),N-n))
  margRest <- array(NA,c(reps,length(v_inflate),N-n))
  margRestHuber <- array(NA,c(reps,length(v_inflate),N-n))
  margT <- array(NA,c(reps,length(v_inflate),N-n))
  

  # Predictions of each y in holdout set  -----
  rlmPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  rlmPredHuberMat <- matrix(NA, nrow=reps, ncol=N-n)
  olsPredMat <- matrix(NA, nrow=reps, ncol=N-n)
  nTheoryPredMat <- array(NA,c(reps,length(v_inflate),N-n))
  restPredMat <- array(NA, c(reps,length(v_inflate),N-n))
  restPredHuberMat <- array(NA, c(reps,length(v_inflate),N-n))
  tPredMat <- array(NA, c(reps,length(v_inflate),N-n))
  
  
  # M-H acceptance rates -----
  acceptY <- array(NA, c(reps,length(v_inflate))) #acceptance rates for augmented y's in restricted model using Tukey's
  acceptYHuber <- array(NA, c(reps,length(v_inflate))) #acceptance rates for augmented y's in restricted model using Huber's
  acceptT <- array(NA, c(reps,length(v_inflate))) #acceptance rates for MH in the tmodel
  
  
  
  
  #auxilary functions and constants ---- 
  #t density with center and scale and df = nu
  tdensity<-function(y, mean, sigma, nu){
    (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
  }
  
  
  # simulation -----
 # system.time(  
    for(i in sims){
      
      trainIndices <- sort(sample(1:N, n))
      holdIndices <- c(1:N)[-trainIndices]
      holdIndicesMatrix[i,] <- holdIndices
      train <- analysis_data[trainIndices,]
      hold <- analysis_data[holdIndices,]
      yholdout <- hold$sqrt_count_2012
      y_hold[i,] <- yholdout
      
      
      #type 1, open agencies to predict in holdout set
      open <- hold$Count_2012 > 0 
      type1 <- hold$Type == '1'
      y_open[i,] <-open
      y_type1[i,] <-type1
      
      
      #regresssions ----
      y <- train$sqrt_count_2012 
      
      #get model matrix from holdoutset
      #fit1 is on the holdoutset
      fit1 <- lm(sqrt_count_2012 ~ sqrt_count_2010 - 1,data = hold)
      Xholdout <- model.matrix(fit1)
      
      
      #rlm on training:Tukey -----
      rlmfit <- rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.bisquare, scale.est='Huber',data = train, maxit=1000)
      X <- model.matrix(rlmfit) #model matrix for all regressions
      # p <- ncol(X) #number of regression coefs
      sigma2Int <- rlmfit$s^2 #to start mcmc
      rlmEstimates[,i] <- c(coef(rlmfit),rlmfit$s^2)
      rlmPreds <- Xholdout%*%coef(rlmfit)
      #rlmPreds <- predict(rlmfit, newdata = data.frame(sqrt_count_2010 = Xholdout))
      rlmPredMat[i,] <- rlmPreds
      margRlm[i,] <- dnorm(yholdout,mean=rlmPreds, sd=rlmfit$s)
     
      
      #rlm on training: Huber ----
      rlmfitHuber <- rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, psi=psi.huber, scale.est='Huber',data = train, maxit=1000)
      
      rlmEstimatesHuber[,i]<-c(coef(rlmfitHuber),rlmfitHuber$s^2)
      rlmPredsHuber<-Xholdout%*%coef(rlmfitHuber)
      rlmPredHuberMat[i,] <- rlmPredsHuber
      margRlmHuber[i,] <- dnorm(yholdout,mean=rlmPredsHuber, sd=rlmfitHuber$s)
      
      
      #ols on training ----
      olsFit <- lm(sqrt_count_2012 ~ sqrt_count_2010 - 1,data=train)
      olsEstimates[,i] <- c(coef(olsFit),summary(olsFit)$sigma^2)
      olsPreds <- Xholdout%*%coef(olsFit)
      olsPredMat[i,] <- olsPreds
      margOls[i,] <- dnorm(yholdout,mean=olsPreds, sd=summary(olsFit)$sigma)
      
#bayes models

for(beta_var_inflate in v_inflate){
v_ind <- which(beta_var_inflate == v_inflate)   
        #summary(prior_fit)
        beta_0 <- coef(prior_fit)
        se_beta_0 <- vcov(prior_fit)^.5
        var_scalar <- floor(nrow(prior_data)*beta_var_inflate/100)
        var_beta_0 <- var_scalar*se_beta_0^2 
        sigma2_hat <- prior_fit$s^2
        
        a_0 <- 5
        b_0 <- sigma2_hat*(a_0 - 1)      
      
            
      #normal theory bayes model ----
      
      nTheory <- brlm::bayesLm(y, X , 
                               mu0 = beta_0, 
                               Sigma0 = var_beta_0, 
                               a0 = a_0, 
                               b0 = b_0, 
                               sigma2Int = sigma2Int, 
                               nkeep = nkeep, 
                               nburn= nburn)
      #posterior means of beta
      ##plot(nTheory$mcmc, ask=FALSE, density=FALSE)
      postMeansNtheory <- colMeans(nTheory$mcmc)
      nTheoryEstimates[,i, v_ind] <- postMeansNtheory
      betaNTheory <- postMeansNtheory[1:p]
      sigma2Ntheory <- postMeansNtheory[p+1]
      nTheoryPreds <- Xholdout%*%betaNTheory
      nTheoryPredMat[i, v_ind,] <- nTheoryPreds
     
      # get marginals f(y_h) for each element in holdout set 
      nMuMatrix <- Xholdout%*%(t(nTheory$mcmc)[1:p,])
      #sd's across samples for each houldout set
      nSigmaMat <- matrix(sqrt(rep(t(nTheory$mcmc)[p+1,],N-n)),N-n,nkeep, byrow = TRUE)
      #nSigmaMat2<-(sqrt(t(nTheory$mcmc)[p+1,]))
      #this is estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for normal model
      margNTheory[i,v_ind,] <- rowMeans(dnorm(yholdout,mean = nMuMatrix, sd = nSigmaMat))
      
      converge_matrix[1,, v_ind, i] <- geweke.diag(nTheory$mcmc)$z
      
        
      # restricted likelihood -----
      
      
      #Tukey version ----
      restricted <- brlm::restrictedBayesLm(y, X, 
                                            regEst='Tukey', 
                                            scaleEst='Huber',
                                            mu0 = beta_0, 
                                            Sigma0 = var_beta_0, 
                                            a0 = a_0, 
                                            b0 = b_0, 
                                            sigma2Int = sigma2Int, 
                                            nkeep = nkeep, 
                                            nburn= nburn, 
                                            maxit=1000)
      #plot(restricted$mcmc, ask=FALSE, density=FALSE)
      postMeansRest <- colMeans(restricted$mcmc)
      restrictedEstimates[,i,v_ind] <- postMeansRest
      betaRest <- postMeansRest[1:p]
      sigma2Rest <- postMeansRest[p+1]
      restPreds <- Xholdout%*%betaRest
      restPredMat[i,v_ind,] <- restPreds
      acceptY[i,v_ind] <- mean(restricted$yAccept)
      
      
      # get marginals f(y_h) for each element in holdout set 
      #means across samples for each holdout set
      restMuMatrix <- Xholdout%*%(t(restricted$mcmc)[1:p,])
      #sd's across samples for each houldout set
      restSigmaMat <- matrix(sqrt(rep(t(restricted$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
      
      
      #estiamate (L(y_h)) of the marginal f(y_h) for each y_h in the holdout set for  restricted model
      margRest[i,v_ind,] <- rowMeans(dnorm(yholdout,mean=restMuMatrix, sd=restSigmaMat))
      # plot(margRest[i,v_ind,], margRlm[i,])
      # abline(0,1)
      # mean(margRest[i,v_ind,]); mean(margRlm[i,])
      converge_matrix[2,, v_ind, i] <- geweke.diag(restricted$mcmc)$z
      
      
      #Huber version ----
      restrictedHuber <- brlm::restrictedBayesLm(y, X,
                                                 regEst='Huber',
                                                 scaleEst='Huber',
                                                 mu0 = beta_0, 
                                                 Sigma0 = var_beta_0, 
                                                 a0 = a_0, 
                                                 b0 = b_0, 
                                                 sigma2Int = sigma2Int, 
                                                 nkeep = nkeep, 
                                                 nburn= nburn, 
                                                 maxit=1000)
      
      
      postMeansRestHuber <- colMeans(restrictedHuber$mcmc)
      restrictedEstimatesHuber[,i,v_ind] <- postMeansRestHuber
      betaRestHuber <- postMeansRestHuber[1:p]
      sigma2RestHuber <- postMeansRestHuber[p+1]
      restPredsHuber <- Xholdout%*%betaRestHuber
      restPredHuberMat[i,v_ind,] <- restPredsHuber
      acceptYHuber[i,v_ind] <- mean(restrictedHuber$yAccept)
      
      #means across samples for each holdout set
      restMuMatrixHuber <- Xholdout%*%(t(restrictedHuber$mcmc)[1:p,])
      #sd's across samples for each houldout set
      restSigmaMatHuber <- matrix(sqrt(rep(t(restrictedHuber$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
      
      #estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for         restricted model
      margRestHuber[i,v_ind,] <- rowMeans(dnorm(yholdout,mean=restMuMatrixHuber, sd=restSigmaMatHuber))
      
      
      converge_matrix[3, ,v_ind, i] <- geweke.diag(restrictedHuber$mcmc)$z
      
      
      # t-model ----
      
      #Note the change of prior on sigma2. In the other models the variance parameter is given  an IG(a_0,b_0) prior
      #the variance for the t model is (nu/(nu-2))sigma2~IG(a_0, b_0),  implies sigma2=(nu-2)/nu var(Y)~IG(a_0,(nu-2)/nu b_0)
      
      tmodel <- brlm::bayesTdistLm(y, X,
                                   mu0 = beta_0, 
                                   Sigma0 = var_beta_0, 
                                   a0 = a_0, 
                                   b0 = ((nu-2)/nu)*b_0,
                                   parInit = NULL,
                                   nu = nu, 
                                   nkeep = nkeept, 
                                   nburn = nburnt,
                                   rwTune = NULL)  
      postMeansT <- colMeans(tmodel$mcmc)
      tmodelEstimates[,i,v_ind] <- postMeansT
      betaT <- postMeansT[1:p]
      sigma2T <- postMeansT[p+1]
      TPreds <- Xholdout%*%betaT
      tPredMat[i,v_ind,] <- TPreds
      acceptT[i,v_ind] <- mean(tmodel$acceptSigma2)
      
      
      # estimate marginals.
      tMuMatrix <- Xholdout%*%(t(tmodel$mcmc)[1:p,])
      tSigmaMat <- matrix(sqrt(rep(t(tmodel$mcmc)[p+1,],N-n)),N-n,nkeept, byrow = TRUE)
      
      margT[i,v_ind,] <- rowMeans(tdensity(yholdout, mean=tMuMatrix,sigma=tSigmaMat, nu = nu))
      
      converge_matrix[4, ,v_ind, i] <- geweke.diag(tmodel$mcmc)$z
  
      # v_ind <- v_ind + 1
}
print(i)
}
# )  
 
  
  out <- list(y_hold = y_hold,
              y_open = y_open,
              y_type1 = y_type1,
              acceptY = acceptY,
              acceptYHuber = acceptYHuber,
              acceptT = acceptT,
              converge_matrix = converge_matrix,
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
  
write_rds(out, file.path(here::here(), paste0('single_reg_state_', State_keep, '_n_', n_percent, '.rds' )))
} 
print(paste0('finish state ', State_keep))
}

end <- Sys.time()
end - strt
