#########################################
#John Lewis
#Simulation
#data was prepared on NW laptop: 2008-2010 for prior constuction and 2010-2012 for data analysis

#Models described by bayesLM.R (requires linearRegressionHuberAndProposal2.R for the restricted methods)

#Priors constructed using nwdataPaper1PriorContruction.R and results saved in nwdataPaper1PriorConstructionWorkSpace.RData
#rlm model on 2008-2010 data used to set prior parameters for beta and sigma2
#relavent transformation of these parameters made when centering and scaling: see nwdataPaper1PriorContruction.R
#prior parameters suffixed by 'star' to indicate these are after transformation




#Task: one simulation consists of taking a sample of size n for the training set 
#Fit the rlm,ols, normal theory Bayes, restricted Bayes (with Tukeys and another with Hubers location), and t-distribution Bayes model: 6 models in total

#suffixes: Rest uses Tukey's regression estimates, RestHuber uses Huber's regression estimates

#prior parameters mu0Star, Sigma0StarUnscaled, a0Star, b0Star construction in nwdataPaper1PriorConstuction

#########################################
#n's from here (25,50,100,200,500,1000,2000)
n<-1000 #sample size for the training data
run<-4 #for large n; run multiple runs in parallel. Set the seed using the index number of seeds
seeds<-c(3,4,5,6)
reps<-25 #number of repetitions of the simulation. If simulations are split into 4 runs and we want a total on 100 reps, use reps=25
nburn<-2000 #set length of mcmc chains
nkeep<-20000
nkeept<-40000 #for the t-model.


#number of samples for the t-model a bit different to allow for longer chain
#load the MCMC algorithm functions

source('../../linearRegressionHuberAndProposal2.R') #restricted sampling functions for y
source('../../bayesLM.R')

#load prior info; this also loads the analysis set called analysisSetCs (the centered/scaled analysis set)

load("../workSpaces/nwdataPaper1PriorConstructionWorkSpace.RData")


#set seed 
set.seed(seeds[run])

N<-nrow(analysisSetCs)
p<-length(mu0Star)



#prior parameters
Sigma0Star
mu0Star

a0Star
b0Star
b0Star/(a0Star-1) #mean
b0Star^2/((a0Star-1)^2*(a0Star-2)) #var 

# matrices to store
# 1 holdout samples
# 2 logical: open agencies 
# 3 logical: open type 1 agencies
yholdoutMat<-matrix(NA,  nrow=reps, ncol=N-n) #mse's for predictions on entire holdoutset
yOpenMat<-matrix(NA, nrow=reps, ncol=N-n) #on just the open agencies
yOpenType1Mat<-matrix(NA, nrow=reps, ncol=N-n) #on just the open type 1 agencies


#######################
#Estimates:
  #rlm regression estimates for classical robust regression
  #Posterior Means for the Bayesian regressions

rlmEstimates<-matrix(NA, nrow=p+1, ncol=reps) #Tukey's Estimates
rlmEstimatesHuber<-matrix(NA, nrow=p+1, ncol=reps) #Huber's
olsEstimates<-matrix(NA, nrow=p+1, ncol=reps) #ols
nTheoryEstimates<-matrix(NA, nrow=p+1, ncol=reps) #normal theory Bayes
restrictedEstimates<-matrix(NA, nrow=p+1, ncol=reps) #Our method with Tukey's
restrictedEstimatesHuber<-matrix(NA, nrow=p+1, ncol=reps) #Our method with Huber's
tmodelEstimates<-matrix(NA, nrow=p+1, ncol=reps) #thick tailed Bayes


########################
#Matrix to save the holdout indices accross reps
#######################
holdIndicesMatrix<-matrix(NA, nrow=reps, ncol=N-n)

###############################
#Matrices to collect the marginals of each y in holdout set 
margRlm<-matrix(NA, nrow=reps, ncol=N-n)
margRlmHuber<-matrix(NA, nrow=reps, ncol=N-n)
margOls<-matrix(NA, nrow=reps, ncol=N-n)
margNTheory<-matrix(NA, nrow=reps, ncol=N-n)
margRest<-matrix(NA, nrow=reps, ncol=N-n)
margRestHuber<-matrix(NA, nrow=reps, ncol=N-n)
margT<-matrix(NA, nrow=reps, ncol=N-n)



###############################
#Matrices to collect the Predictions of each y in holdout set 
rlmPredMat<-matrix(NA, nrow=reps, ncol=N-n)
rlmPredHuberMat<-matrix(NA, nrow=reps, ncol=N-n)
olsPredMat<-matrix(NA, nrow=reps, ncol=N-n)
nTheoryPredMat<-matrix(NA, nrow=reps, ncol=N-n)
restPredMat<-matrix(NA, nrow=reps, ncol=N-n)
restPredHuberMat<-matrix(NA, nrow=reps, ncol=N-n)
tPredMat<-matrix(NA, nrow=reps, ncol=N-n)


########################
#AcceptanceRate vectors
acceptY<-numeric(reps) #acceptance rates for augmented y's in restricted model using Tukey's
acceptYHuber<-numeric(reps) #acceptance rates for augmented y's in restricted model using Huber's
acceptT<-numeric(reps) #acceptance rates for MH in the tmodel
#############################



##################################
#Define the t density df is fixed at nu
#################################
tdensity<-function(y, mean, sigma){
  (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
}

#find the min sqrt_Count2012: this corressponds to zero policies (and hence a closed agency)
zeroOnCsScale<-min(analysisSetCs$sqrt_Count2012)

#-----------------------
#begin the simulation
#-----------------------
system.time(  
  for(i in 1:reps){
    if(n>=1000){
    #define the training and holdout sets
    #make sure at least 5 per state
    satisfy<-FALSE
    while(!satisfy){  
    temp<-sort(sample(1:N, n))
    anyLess<-sum(table(analysisSetCs[temp, "Primary_Agency_State"])<5)>0
    if(!anyLess){satisfy<-TRUE}
    }
    trainIndices<-temp
    } else {
      trainIndices<-sort(sample(1:N, n))
    }
    
    
    holdIndices<-c(1:N)[-trainIndices]
    holdIndicesMatrix[i,]<-holdIndices
    trainSet<-analysisSetCs[trainIndices,]
    #table(trainSet$Primary_Agency_State)
    holdoutSet<-analysisSetCs[holdIndices,]
    yholdout<-holdoutSet$sqrt_Count2012
    yholdoutMat[i,]<-yholdout
    
    
    ##################
    #Finding agencies in holdout set that are open
    ##################
    openAgencies<-holdoutSet$sqrt_Count2012>zeroOnCsScale
    yOpenMat[i,]<-openAgencies
    ##################
    #Finding type 1 agencies that are open (in holdout set)
    ####################
    openType1<-holdoutSet$sqrt_Count2012>zeroOnCsScale & holdoutSet$Agency_Type=='Type 1'
    yOpenType1Mat[i,]<-openType1
    ######################
    
    ##########
    #response
    ##########
    y<-trainSet$sqrt_Count2012 
    #get model matrix from holdoutset
    #fit1 is on the holdoutset
    fit1<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=holdoutSet)
    Xholdout<-model.matrix(fit1)
 
    ############
    #models
    ############
    #rlm on training:Tukey
    rlmfit<-rlm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure, psi=psi.bisquare, scale.est='Huber',data=trainSet, maxit=1000)
    X<-model.matrix(rlmfit)
    p<-ncol(X)
    sigma2Int<-rlmfit$s^2
    rlmEstimates[,i]<-c(coef(rlmfit),rlmfit$s^2)
    rlmPreds<-Xholdout%*%coef(rlmfit)
    rlmPredMat[i,]<-rlmPreds
    
    margRlm[i,]<-dnorm(yholdout,mean=rlmPreds, sd=rlmfit$s)
    
    
    
    #rlm on training: Huber
    rlmfitHuber<-rlm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure, psi=psi.huber, scale.est='Huber',data=trainSet, maxit=1000)
    X<-model.matrix(rlmfitHuber)
    p<-ncol(X)
    sigma2Int<-rlmfitHuber$s^2
    rlmEstimatesHuber[,i]<-c(coef(rlmfitHuber),rlmfitHuber$s^2)
    rlmPredsHuber<-Xholdout%*%coef(rlmfitHuber)
    rlmPredHuberMat[i,]<-rlmPredsHuber
    margRlmHuber[i,]<-dnorm(yholdout,mean=rlmPredsHuber, sd=rlmfitHuber$s)
  
    
    #ols on training
    olsFit<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=trainSet)
    olsEstimates[,i]<-c(coef(olsFit),summary(olsFit)$sigma^2)
    olsPreds<-Xholdout%*%coef(olsFit)
    olsPredMat[i,]<-olsPreds
    margOls[i,]<-dnorm(yholdout,mean=olsPreds, sd=summary(olsFit)$sigma)
    
    #---------
    #normal theory bayes
    #---------
    #My function bayesLM fits same model as: MCMCfit<-MCMCregress(y~X[,2]+X[,3]+X[,4], b0=mu0Star, B0=round(solve(Sigma0Star),10),c0=2*a0Star, d0=2*b0Star)...MCMCregress is much faster bc it uses C++
    nTheory<-bayesLm(y, X,mu0Star, Sigma0Star, a0Star, b0Star,sigma2Int, nkeep, nburn)
    #posterior means of beta
    ##plot(nTheory$mcmc, ask=FALSE, density=FALSE)
    postMeansNtheory<-colMeans(nTheory$mcmc)
    nTheoryEstimates[,i]<-postMeansNtheory
    betaNTheory<-postMeansNtheory[1:p]
    sigma2Ntheory<-postMeansNtheory[p+1]
    nTheoryPreds<-Xholdout%*%betaNTheory
    nTheoryPredMat[i,]<-nTheoryPreds
    
    #means across samples for each holdout set
    nMuMatrix<-Xholdout%*%(t(nTheory$mcmc)[1:p,])
    #sd's across samples for each houldout set
    nSigmaMat<-matrix(sqrt(rep(t(nTheory$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
    #nSigmaMat2<-(sqrt(t(nTheory$mcmc)[p+1,]))
    #this is estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for normal model
    margNTheory[i,]<-rowMeans(dnorm(yholdout,mean=nMuMatrix, sd=nSigmaMat))
    

    #---------
    #restricted bayes
    #---------
    
    #---------
    #Tukey
    #---------
    restricted<-restrictedBayesLm(y, X,regEst='Tukey',scaleEst='Huber',mu0Star, Sigma0Star, a0Star, b0Star,sigma2Int, nkeep, nburn, maxit=1000)
    # Sys.time()-st
    #posterior means of beta
    ##plot(restricted$mcmc, ask=FALSE, density=FALSE)
    postMeansRest<-colMeans(restricted$mcmc)
    restrictedEstimates[,i]<-postMeansRest
    betaRest<-postMeansRest[1:p]
    sigma2Rest<-postMeansRest[p+1]
    restPreds<-Xholdout%*%betaRest
    restPredMat[i,]<-restPreds
    acceptY[i]<-mean(restricted$yAccept)
    
    
    #means across samples for each holdout set
    restMuMatrix<-Xholdout%*%(t(restricted$mcmc)[1:p,])
    #sd's across samples for each houldout set
    restSigmaMat<-matrix(sqrt(rep(t(restricted$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
    
    
    #estiamate (L(y_h)) of the marginal f(y_h) for each y_h in the holdout set for       restricted model
    margRest[i,]<-rowMeans(dnorm(yholdout,mean=restMuMatrix, sd=restSigmaMat))
    
    ######################################################

    #---------
    #Huber
    #---------
    restrictedHuber<-restrictedBayesLm(y, X,regEst='Huber',scaleEst='Huber',mu0Star, Sigma0Star, a0Star, b0Star,sigma2Int, nkeep, nburn, maxit=1000)
    # Sys.time()-st
    #posterior means of beta
    ##plot( restrictedHuber$mcmc, ask=FALSE, density=FALSE)
    postMeansRestHuber<-colMeans(restrictedHuber$mcmc)
    restrictedEstimatesHuber[,i]<-postMeansRestHuber
    betaRestHuber<-postMeansRestHuber[1:p]
    sigma2RestHuber<-postMeansRestHuber[p+1]
    restPredsHuber<-Xholdout%*%betaRestHuber
    restPredHuberMat[i,]<-restPredsHuber
    acceptYHuber[i]<-mean(restrictedHuber$yAccept)
    
    #means across samples for each holdout set
    restMuMatrixHuber<-Xholdout%*%(t(restrictedHuber$mcmc)[1:p,])
    #sd's across samples for each houldout set
    restSigmaMatHuber<-matrix(sqrt(rep(t(restrictedHuber$mcmc)[p+1,],N-n)),N-n,nkeep, byrow=TRUE)
    
    #estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for         restricted model
    margRestHuber[i,]<-rowMeans(dnorm(yholdout,mean=restMuMatrixHuber, sd=restSigmaMatHuber))
      
    #---------
    #Thick Tailed Bayesian Model
    #---------
  
    #Note: change the prior on sigma2: in the other models the variance parameter is given   an IG(a0star,b0star) prior
    
    #the variance here is (nu/(nu-2))sigma2~IG(a0star, b0Star),  implies sigma2=(nu-2)/nu var(Y)~IG(a0star,(nu-2)/nu b0star)
    
    tmodel<-bayesTdistLm2(y, X,mu0Star, Sigma0Star, a0Star,((nu-2)/nu)*b0Star,parInit=NULL,nu, nkeep=nkeept, nburn=nburn,rwTune=NULL)  
    postMeansT<-colMeans(tmodel$mcmc)
    tmodelEstimates[,i]<-postMeansT
    betaT<-postMeansT[1:p]
    sigma2T<-postMeansT[p+1]
    TPreds<-Xholdout%*%betaT
    tPredMat[i,]<-TPreds
    acceptT[i]<-tmodel$accept
    
    
    tMuMatrix<-Xholdout%*%(t(tmodel$mcmc)[1:p,])
    tSigmaMat<-matrix(sqrt(rep(t(tmodel$mcmc)[p+1,],N-n)),N-n,nkeept, byrow=TRUE)
    
    ########################
    #estiamate L(y_h) of the marginal f(y_h) for each y_h in the holdout set for t model 
    #######################
    margT[i,]<-rowMeans(tdensity(yholdout, mean=tMuMatrix,sigma=tSigmaMat))
    print(i)
    
  }
  
    )


##########################
#save the mse matrices
#########################
rm(list=setdiff(ls(),c('yholdoutMat','yOpenMat','yOpenType1Mat','acceptY','acceptYHuber','acceptT', 'n','rlmEstimates',"rlmEstimatesHuber",'nTheoryEstimates','restrictedEstimates',"restrictedEstimatesHuber",'tmodelEstimates','olsEstimates', 'mu0Star', 'Sigma0Star','margNTheory','margRest',"margRestHuber",'margT','margRlm',"margRlmHuber",'margOls','holdIndicesMatrix',"rlmPredMat", 'rlmPredHuberMat','olsPredMat','nTheoryPredMat','restPredMat','restPredHuberMat','tPredMat', 'reps', 'run')))




assign(paste('yholdoutMat',n, sep=''), yholdoutMat)
assign(paste('yOpenMat',n, sep=''), yOpenMat)
assign(paste('yOpenType1Mat',n, sep=''),yOpenType1Mat)
rm(yholdoutMat,yOpenMat,yOpenType1Mat)

assign(paste('acceptY',n, sep=''),acceptY)
assign(paste('acceptYHuber',n, sep=''), acceptYHuber)
assign(paste('acceptT',n, sep=''),acceptT)
rm(acceptY,acceptYHuber,acceptT)


assign(paste('margRlm',n, sep=''), margRlm)
assign(paste('margRlmHuber',n, sep=''), margRlmHuber)
assign(paste('margOls',n, sep=''),margOls)
assign(paste('margNTheory',n, sep=''), margNTheory)
assign(paste('margRest',n, sep=''), margRest)
assign(paste('margRestHuber',n, sep=''), margRestHuber)
assign(paste('margT',n, sep=''),margT)
assign(paste('holdIndicesMatrix',n, sep=''), holdIndicesMatrix)
rm(margRlm, margRlmHuber, margOls, margNTheory, margRest, margRestHuber, margT, holdIndicesMatrix)



assign(paste('rlmPredMat',n, sep=''), rlmPredMat)
assign(paste('rlmPredHuberMat',n, sep=''), rlmPredHuberMat)
assign(paste('olsPredMat',n, sep=''),olsPredMat)
assign(paste('nTheoryPredMat',n, sep=''),nTheoryPredMat)
assign(paste('restPredMat',n, sep=''), restPredMat)
assign(paste('restPredHuberMat',n, sep=''),restPredHuberMat)
assign(paste('tPredMat',n, sep=''),tPredMat)
rm(rlmPredMat, rlmPredHuberMat, olsPredMat, nTheoryPredMat, restPredMat, restPredHuberMat, tPredMat)



assign(paste('rlmEstimates',n, sep=''), rlmEstimates)
assign(paste('rlmEstimatesHuber',n, sep=''), rlmEstimatesHuber)
assign(paste('olsEstimates',n, sep=''),olsEstimates)
assign(paste('nTheoryEstimates',n, sep=''),nTheoryEstimates)
assign(paste('restrictedEstimates',n, sep=''), restrictedEstimates)
assign(paste('restrictedEstimatesHuber',n, sep=''),restrictedEstimatesHuber)
assign(paste('tmodelEstimates',n, sep=''),tmodelEstimates)
rm(rlmEstimates, rlmEstimatesHuber, olsEstimates, nTheoryEstimates, tmodelEstimates,restrictedEstimates, restrictedEstimatesHuber)

########################
#save the workspace
########################
if(n %in% c(1000,1500, 2000,2500)){
save.image(paste('../workSpaces/wsNwAnalysis_n',n,'run', run, '.RData', sep=''))
} else {
  save.image(paste('../workSpaces/wsNwAnalysis_n',n, '.RData', sep=''))
}
