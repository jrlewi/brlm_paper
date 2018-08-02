# Pooled regression - i.e., single regression without accounting for variation by state


# library(devtools)
# install_github('brlm', 'jrlewi')
library(MCMCpack)
library(MASS)
library(brlm)
library(tidyverse)

#load data and prior information

analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))
parms_prior <- read_rds(file.path(here::here(), 'parms_prior.rds'))
analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012))
# 
#  fit_2 <- MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, data = analysis_data) #filter(analysis_data, Type == '1', Count_2012 > 0))
#  fit_22 <- MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, data = filter(analysis_data, Type == '1', Count_2012 > 0))
#  



beta_0 <- parms_prior$beta_0
se_beta_0 <- parms_prior$se_beta_0
var_beta_0 <- se_beta_0^2 #inflating prior info
sigma2_hat_0 <- parms_prior$sigma2_hat


sigma2_hat_0
fit_2$s^2
fit_22$s^2

xx <-seq(.94, .96, by = .0005)
plot(xx, dnorm(xx, beta_0, var_beta_0^.5), type = 'l')
# abline(v = coef(fit_2))
# abline(v = coef(fit_22), col = 2)
nu <- 3 #df for t-model

#n <- 1000  
ns <- c(25, 100)# , 1000) #sample size for training set
reps <- 10 # number of training sets

nburn <- 2000 #set length of mcmc chains
nkeep <- 2000
nkeept <- 4000 #for the t-model.

 # fit_2 <- MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, data = filter(analysis_data, Type == '1', Count_2012 > 0))
 # fit_2 <- MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, data = analysis_data)

#set seed 
set.seed(123)
N <- nrow(analysis_data)
p <- length(beta_0)



#prior parameters for sigma2 ----
a_0 <- 5
b_0 <- sigma2_hat_0*(a_0 - 1)

b_0/(a_0-1) # prior mean for sigma^2
b_0^2/((b_0-1)^2*(b_0-2)) #var 

for(n in ns){

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
nTheoryEstimates <- matrix(NA, nrow=p+1, ncol=reps) #normal theory Bayes
restrictedEstimates <- matrix(NA, nrow=p+1, ncol=reps) #Our method with Tukey's
restrictedEstimatesHuber <- matrix(NA, nrow=p+1, ncol=reps) #Our method with Huber's
tmodelEstimates <- matrix(NA, nrow=p+1, ncol=reps) #thick tailed Bayes


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
    p <- ncol(X) #number of regression coefs
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
    nTheoryEstimates[,i] <- postMeansNtheory
    betaNTheory <- postMeansNtheory[1:p]
    sigma2Ntheory <- postMeansNtheory[p+1]
    nTheoryPreds <- Xholdout%*%betaNTheory
    nTheoryPredMat[i,] <- nTheoryPreds
    
    # get marginals f(y_h) for each element in holdout set 
    nMuMatrix <- Xholdout%*%(t(nTheory$mcmc)[1:p,])
    #sd's across samples for each houldout set
    nSigmaMat <- matrix(sqrt(rep(t(nTheory$mcmc)[p+1,],N-n)),N-n,nkeep, byrow = TRUE)
    #nSigmaMat2<-(sqrt(t(nTheory$mcmc)[p+1,]))
    #this is estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for normal model
    margNTheory[i,] <- rowMeans(dnorm(yholdout,mean = nMuMatrix, sd = nSigmaMat))
    
    
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
    
tmodel <- brlm::bayesTdistLm(y, X,
                           mu0 = beta_0, 
                           Sigma0 = var_beta_0, 
                           a0 = a_0, 
                           b0 = ((nu-2)/nu)*b_0,
                           parInit = NULL,
                           nu = nu, 
                           nkeep = nkeept, 
                           nburn = nburn,
                           rwTune = NULL)  
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

# save output
# rm(list=setdiff(ls(),c('y_hold','yOpenMat','y_open_type1','acceptY','acceptYHuber','acceptT', 'n','rlmEstimates',"rlmEstimatesHuber",'nTheoryEstimates','restrictedEstimates',"restrictedEstimatesHuber",'tmodelEstimates','olsEstimates', 'mu0Star', 'Sigma0Star','margNTheory','margRest',"margRestHuber",'margT','margRlm',"margRlmHuber",'margOls','holdIndicesMatrix',"rlmPredMat", 'rlmPredHuberMat','olsPredMat','nTheoryPredMat','restPredMat','restPredHuberMat','tPredMat', 'reps', 'run')))


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