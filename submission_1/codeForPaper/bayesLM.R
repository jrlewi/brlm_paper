#############################
#Code for the standard linear model fits:
#Functions include data models of:

# Normal
# restricted with normal prior as base
# t-distributed data with nu degrees of freedom
###############################################

#priors are the same for beta and sigma2 throughout

#############################
library(MCMCpack)
#library(msm) 
library(Matrix)
library(mcmc)
library(MASS)
#######
#Model: standard linear regression model
#######

#y_i=beta*x_i+e_i

#y_i|b,sigma2~N(x_i'beta,sigma^2)

#e_i~N(0, sigma^2)

#b~N(mu0,Sigma0)

#sigma^2~IG(a0,b0)

#mu0,Sigma0, a0, b0 are fixed

######################
#Function for MCMC
######################
#matches MCMCregress(y~X[,2]+X[,3]+X[,4]+..., b0=mu0, B0=round(solve(Sigma0),10),c0=2*a0, d0=2*b0, mcmc=, burin) from MCMCpack
#MCMCregress is much faster becuase it uses C

bayesLm<-function(y, X,mu0, Sigma0, a0, b0,sigma2Int, nkeep=1e4, nburn=1e3){
  #y is the response
  #X is the design Matrix
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  #sigma2Int is the initial value for sigma2
  #nkeep: length of final MCMC chain
  #nburn: length for burnin
  
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  #set inital values to current values
  sigma2Cur<-sigma2Int
  Sigma0Inv<-solve(Sigma0)
  
  #################################
  #Sampling Functions for Gibbs Sampler
  #################################
  #[sigma^2|---]=IG(a0+n/2, b0+.5*t(Y-Xbeta)(Y-Xbeta))
  #beta|----]=N(mun, Lambdan)
  #mun=Lambdan%*%(t(X)%*%Y+Sigma0^-1mu0)/sigma^2
  #Lambdan=sigma^2(t(X%*%X)+Sigma0^-1)^-1; (t(X%*%X)+Sigma0^-1)^-1=varCov_0
  
  #a0 b0 must be specified
  sampleSigma2<-function()
  {
    an<-a0+.5*n
    resid<-y-X%*%betaCur
    bn<-as.numeric(b0+.5*crossprod(resid,resid))
    return(rinvgamma(1,an,bn))
  }

  
  sampleBeta<-function(){
    varCov.n<-solve(t(X)%*%X/sigma2Cur+Sigma0Inv)
    mu.n<-varCov.n%*%(t(X)%*%y/sigma2Cur+Sigma0Inv%*%mu0)
    return(mvrnorm(1,mu.n,varCov.n))
  }
################################################
#Start the gibbs sampler
################################################  
  for(i in 1:total){
    #step one: [beta|sigma2,y]~N()
    betaCur<-sampleBeta()
    #step 2 [sigma2|---]~IG
    sigma2Cur<-sampleSigma2()
    betaSamples[i,]<-betaCur
    sigma2Samples[i]<-sigma2Cur
  }
  
  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
  colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out
}



################################
#Model: Restricted model with normal Theory Model above as the base
################################

######################
#Function for for restricted MCMC
######################


restrictedBayesLm<-function(y, X,regEst='Huber',scaleEst='Huber',mu0, Sigma0, a0, b0,sigma2Int, nkeep=1e4, nburn=1e3, maxit=400){
  #y is the response
  #X is the design Matrix
  #regEst: regression estimates used: options are 'Huber' and 'Tukey'
  #scaleEst='Huber' only option is Huber's proposal two
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  # sigma2Int is the initial value for sigma2
  #nkeep: length of final MCMC chain
  #nburn: length for burnin
  #maxit maximum iterations for computing robust regression estimates
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  Sigma0Inv<-solve(Sigma0)
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  yAccept<-numeric(total)
  #set inital values to current values
  sigma2Cur<-sigma2Int
  
  
  
  Q<-qr.Q(qr(X))
  projMatrix<-diag(n)-tcrossprod(Q,Q) #Q%*%t(Q)
  
  ############################
  #define the psi and chi functions
  ############################
  if(regEst=='Huber') {
    psi<-get('psi.huber') #internal
    fn.psi<-get('fn.psi.huber')
    
  } else { 
    if(regEst=='Tukey'){
      psi<-get('psi.bisquare') #internal
      fn.psi<-get('fn.psi.bisquare')
    } else {stop("only set up for Huber or Tukey regression estimates")}}
  
  if(scaleEst!='Huber'){
    stop('scale estimate only set up for Hubers Prop2 ')
  }
  fn.chi<-fn.chi.prop2
  ################################

  
  #run the robust regression
  robust<-rlm(X,y,psi=psi, scale.est=scaleEst, maxit=maxit)
  #condition on these estimates
  l1obs<-robust$coefficients
  s1obs<-robust$s
  
  #################################
  #Sampling Functions for Gibbs Sampler. Those unique to the restricted sampling are in the source code
  #################################
  #[sigma^2|---]=IG(a0+n/2, b0+.5*t(Y-Xbeta)(Y-Xbeta))
  #beta|----]=N(mun, Lambdan)
  #mun=Lambdan%*%(t(X)%*%Y+Sigma0^-1mu0)/sigma^2
  #Lambdan=sigma^2(t(X%*%X)+Sigma0^-1)^-1; (t(X%*%X)+Sigma0^-1)^-1=varCov_0
  
  #a0 b0 must be specified
  sampleSigma2<-function(data)
  {
    an<-a0+.5*n
    resid<-data-X%*%betaCur
    bn<-as.numeric(b0+.5*crossprod(resid,resid))
    return(rinvgamma(1,an,bn))
  }
  
  sampleBeta<-function(data){
    varCov.n<-solve(t(X)%*%X/sigma2Cur+Sigma0Inv)
    mu.n<-varCov.n%*%(t(X)%*%data/sigma2Cur+Sigma0Inv%*%mu0)
    return(mvrnorm(1,mu.n,varCov.n))
  }
  #choose starting y data
  y.prop<-rnorm(n)
  y.curr<-fn.comp.ystst(y.prop,X,l1obs,s1obs,psi,scaleEst,maxit)
  log.prop.den.curr <-log.prop.den(y.curr,X, projMatrix, l1obs, s1obs,fn.psi, fn.chi,n,p)
  ################################################
  #Start the gibbs sampler
  ################################################  
  for(i in 1:total){
    #step one: [beta|sigma2,y]~N()
    betaCur<-sampleBeta(y.curr)
    #step 2 [sigma2|---]~IG
    sigma2Cur<-sampleSigma2(y.curr)
    betaSamples[i,]<-betaCur
    sigma2Samples[i]<-sigma2Cur
    #step 3: augmented data 
    ySample<-fn.one.rep.y(y.curr,betaCur,sqrt(sigma2Cur),l1obs, s1obs,X, log.prop.den.curr, projMatrix,fn.psi,fn.chi, psi,scaleEst,maxit)
    y.curr<-ySample[[1]]
    yAccept[i]<-ySample[[2]]
    log.prop.den.curr<-ySample[[3]]
    }
 
  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
  colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out$yAccept<-yAccept
  out$robust<-robust
  out$coef<-l1obs
  out$s<-s1obs
  
  out
}



#######
#Model: linear regression model with t distributed errors, no hierarchy
#######

#y_i=b_j*x_ij+e_i

#y_i|b,sigma2~t_v(x_i'b,sigma^2)

#same priors

#b_i~N(mu0,Sigma0)

#sigma^2~IG(a0,b0)

#mu0,Sigma0, a0, b0 are fixed


######################
#Function for MCMC for T
######################
#method 1:
#t-data model parameterized as follows
#see Gelman page 303 
# y_i~N(t(x_i)beta,v_i)
# v_i~inv-chi^2(nu,sigma2) #nu fixed
# 
# beta~N(mu0,Sigma_0)
# sigma2~IG(a0,b0)

bayesTdistLm<-function(y, X,mu0, Sigma0 , a0, b0,parInit=NULL,nu, nkeep=1e4, nburn=1e3,rwTune=NULL){
  #y is the response
  #X is the design Matrix
  # parInit: initial values for c(beta,sigma2)
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  #rwTune=MH rando  m walk for sigma2 st.dev 
 
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  Sigma0Inv<-solve(Sigma0)

if(is.null(parInit)){
  fit1<-rlm(X,y, maxit=400)
  parInit<-c(coef(fit1),summary(fit1)$sigma^2)
}
betaCur<-parInit[1:p]
sigma2Cur<-parInit[p+1]

if(is.null(rwTune)){
  fit1<-rlm(X,y, maxit=400)
  rwTune<-summary(fit1)$sigma/5
}

  #Vectors to save samples
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  acceptSigma2<-numeric(total)

############
#Sampling functions
############
#for V
sampleV<-function(){
    as.vector(((y-X%*%betaCur)^2+nu*sigma2Cur)/rchisq(n,nu+1))
  }
#for beta
sampleBeta<-function()
{ 
  vInv<-diag(vCur^-1)
  Var<-solve(t(X)%*%vInv%*%X+Sigma0Inv)
  Mean<-Var%*%(t(X)%*%vInv%*%y+Sigma0Inv%*%mu0)
  return(mvrnorm(1,mu=Mean, Sigma=Var))
}

#Sampling for sigma2
#randomWalkMH for sigma2
#log posterior of [sigma2|everything else]
lpSigma2<-function(par){
    if(par>0){
      return((-a0-1+0.5*n*nu)*log(par)-b0/par-.5*nu*par*sum(1/vCur))}
    else return(-Inf)
  }
mhStepSigma2<-function(){
sigma2Prop<-rnorm(1, sigma2Cur, sd=rwTune)
logMHrat<-lpSigma2(sigma2Prop)-lpSigma2(sigma2Cur)
aprob<-min(1, exp(logMHrat))
if(runif(1)<aprob){
  return(c(sigma2Prop,1))
} else{
  return(c(sigma2Cur,0))
}
 }
############################

############################
#Start the Gibbs Sampler
############################

for(i in 1:total){
  vCur<-sampleV()
  betaCur<-sampleBeta()
  outSigma2<-mhStepSigma2()
  sigma2Cur<-outSigma2[1]
  betaSamples[i,]<-betaCur
  sigma2Samples[i]<-sigma2Cur
  acceptSigma2[i]<-outSigma2[2]
}
#output
colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out$acceptSigma2<-acceptSigma2
  out
}


#Second way to sample from the posterior using the t-distribution
#may be a bit faster, but more correlation
bayesTdistLm2<-function(y, X,mu0, Sigma0, a0, b0,parInit=NULL,nu, nkeep=1e4, nburn=1e3,rwTune=NULL){
  #y is the response
  #X is the design Matrix
  # parInit is the initial values for c(beta,sigma2)
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  #rwTune=MH random walk for all the params 
  
  
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  Sigma0Inv<-solve(Sigma0)
  
  if(is.null(parInit)){
    fit1<-rlm(X,y, maxit=400)
    parInit<-c(coef(fit1),summary(fit1)$sigma^2)
  }
  if(is.null(rwTune)){
    #fit1<-rlm(X,y)
    fit1<-lm(y~X[,-1])
    rwTune<-c(diag(summary(fit1)$sigma^2*solve(t(X)%*%X))^.5/5,summary(fit1)$sigma/100)
  }

#log posterior:
  lp<-function(par){if(par[p+1]>0){
   return(-n*log(sqrt(par[p+1]))-.5*(nu+1)*sum(log(1+(1/nu)*((y-X%*%par[1:p])^2/par[p+1])))-.5*(t(par[1:p]-mu0)%*%Sigma0Inv%*%(par[1:p]-mu0))-(a0+1)*log(par[p+1])-b0/par[p+1])} else{
      return(-Inf)
    }
  }
  outMH<-metrop(lp, initial=parInit, nbatch=total, scale=rwTune)
  #out$accept
  #plot(as.mcmc(out$batch),ask=FALSE,density=FALSE)
#   #output
#
  betaSamples<-outMH$batch[((nburn+1):total),1:p]
  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples,outMH$batch[((nburn+1):total),p+1]))
 colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
 fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
 out<-list()
 out$mcmc<-mcmcBetaSigma2
 out$fitted.values<-fittedValues
 out$accept<-outMH$accept
out
}
