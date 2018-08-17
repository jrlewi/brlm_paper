#
# Simulation ------- 
#

#requires sourcing linearRegressionHuberAndProposal2.R #path may need to be changed
# Adjust source for Steve's directories at work
source('linearRegressionHuberAndProposal2.R')
# source('//trad/profs$/maceachern.1/Desktop/LEWIS.LEE/codeForPaper/codeForPaper/linearRegressionHuberAndProposal2.R')


#
# Complete Model ----
# theta_i \sim N(mu, tau^2)
# Y_ij | \theta_i \sim N(\theta_i, sigma_i^2)  
# with a fairly flat distribution on mu, tau^2 and sigma_i^2.
# mu \propto 1, \tau^2 \propto 1/\tau^2 = IG(0,0), sigma_i^2 \~ IG(ass, bss)
# i runs from 1 to 90
#

#
# Data Generation ----
# see fn.gen.data
# theta_i~N(0,1) i=1,...,90
# Y_ij | \theta_i \sim (1-p) N(theta_i, sigma2)+ p N(theta_i,m*sigma2), j=1,...,n_i

#90 groups formed by factorial stucture of the 3 factors
# n_i=25, 50, 100
# p= .1,.2,.3
# m= 9, 25
# each combination observed in five groups
# sigma2 also will vary: .5, 1, 4, 10

#
# 
# Sampling functions---
#
#

# [tau^2 | rest ] ----

# tau^2 \sim IG(a,b) in functional form with a=b=0
fn.sample.tau2<-function(theta, mu, a=0, b=0){
  #theta is the vector of theta_i
  K<-length(theta)
  anew<-0.5*K+a
  bnew<-0.5*sum((theta-mu)^2)+b
  1/rgamma(1, shape=anew, rate=bnew)
}

# [mu | rest ] ----

fn.sample.mu<-function(theta, tau2){
  K<-length(theta)
  rnorm(1, mean(theta),sqrt(tau2/K))
}

# [sigma2_i| rest ] ----
# sigma2_i \sim IG(a,b) in functional form with a=b=0

# fn.sample.sigma2.i<-function(Yi, thetai, a=0,b=0){
fn.sample.sigma2.i<-function(Yi, thetai, ass , bss){  # SteveMod
  # Yi: data vector from ith group
  # thetai: ith theta
  n<-length(Yi)
#  anew<-0.5*n+a
#  bnew<-0.5*sum((Yi-thetai)^2)+b
  anew<-0.5*n+ass # SteveMod
  bnew<-0.5*sum((Yi-thetai)^2)+bss # SteveMod
  1/rgamma(1, shape=anew, rate=bnew)
}


# [theta_i| rest ]
fn.sample.theta.i<-function(Yi,sigma2i, mu, tau2){
  # Yi: data vector from ith group
  # sigma2i: variance parameter from ith group
  n<-length(Yi)
  muNew<-(tau2*sum(Yi)+mu*sigma2i)/(sigma2i+n*tau2)
  varNew<-sigma2i*tau2/(sigma2i+n*tau2)
  rnorm(1, muNew, sqrt(varNew))
}
  

#
#
# End Sampling Functions---
#
#

# One Rep complete MCMC ----
# One rep of complete data MCMC of model written above under ' Complete Model'


# fn.one.rep.complete.MCMC<-function(YList, theta, sigma2=NULL, mu, tau2=NULL,a=0,b=0){
fn.one.rep.complete.MCMC<-function(YList, theta, sigma2=NULL, mu, tau2=NULL,a=0,b=0, ass, bss){ # SteveMod
# YList: list consisting of a vector of data from each group
# theta: vector current values of each theta_i
# sigma2: vector current values of each sigma2_i
# lengths of YList, theta, and sigma2 should be equal
# mu: current value of mu
# tau2: current value of tau2
  #a,b: prior parameters for inverse gamma on tau^2
# Note: A change in the ordering would require a change in the required inputs. Inputs not needed for initiation are defaulted to NULL
  
  #sample tau2
  tau2<-fn.sample.tau2(theta, mu,a=a,b=b)
  
  #sample mu
  mu<-fn.sample.mu(theta, tau2)
  
  #sample each sigma2_i
#  sigma2<-mapply(FUN=fn.sample.sigma2.i, YList, theta)
  sigma2<-mapply(FUN=fn.sample.sigma2.i, YList, theta, ass, bss) # SteveMod
  #sample each theta_i
  theta<-mapply(FUN=fn.sample.theta.i, YList,sigma2, mu, tau2)
  out<-list()
  out$theta<-theta
  out$sigma2<-sigma2
  out$mu<-mu
  out$tau2<-tau2
  out
 
  }  
  
    
# MCMC function for Complete model ----
# Fits the model given above
#fn.complete.MCMC<-function(YList, nkeep, nburn, theta, sigma2=NULL, mu, tau2=NULL, printEach=1000,a=0,b=0){
fn.complete.MCMC<-function(YList, nkeep, nburn, theta, sigma2=NULL, mu, tau2=NULL, printEach=1000,a=0,b=0,ass,bss){ #SteveMod
  # YList: list consisting of a vector of data from each group
  # theta: vector of initial values of each theta_i
  # sigma2: vector of intitial values of each sigma2_i
  # lengths of YList, theta, and sigma2 should be equal
  # mu: initial value of mu
  # tau2: initial value of tau2
  # a,b: prior parameters for inverse gamma on tau^2
  # Note: A change in the ordering would require a change in the required inputs. Inputs not needed for initiation are defaulted to NULL
  #printEach: print progess every printEach interations
  
  nGroups<-length(YList)
  n<-sapply(YList, length)
  nTot<-nkeep+nburn
  #storage vectors and matrices
  tau2.sample<-numeric(nTot)
  mu.sample<-numeric(nTot)
  theta.sample<-matrix(NA, nrow=nTot, ncol=nGroups)
  sigma2.sample<-matrix(NA, nrow=nTot, ncol=nGroups)
  
  for(it in 1:nTot){
#    samp<-fn.one.rep.complete.MCMC(YList, theta, sigma2, mu, tau2,a=a,b=b)
samp<-fn.one.rep.complete.MCMC(YList, theta, sigma2, mu, tau2,a=a,b=b, ass=ass, bss=bss) # SteveMod
    theta<-samp$theta
    sigma2<-samp$sigma2
    mu<-samp$mu
    tau2<-samp$tau2
    
    theta.sample[it,]<-theta
    sigma2.sample[it,]<-sigma2
    mu.sample[it]<-mu
    tau2.sample[it]<-tau2
    
    out<-list()
    keep<-(nburn+1):nTot
    out$theta<-theta.sample[keep,]
    out$sigma2<-sigma2.sample[keep,]
    out$mu<-mu.sample[keep]
    out$tau2<-tau2.sample[keep]
    if(it %% printEach==0) print(paste('iteration', it, sep=' '))
  }
out  
}
  
#
#
# MCMC function for Incomplete (restricted) model ----
#
#


fn.Incomplete.MCMC<-function(YList,
                             regEst='Huber',
                             scaleEst='Huber',
                             nkeep, 
                             nburn, 
                             theta, 
                             sigma2=NULL, 
                             mu, 
                             tau2=NULL, 
                             printEach=1000,
                             maxit=400,
                             a=0,
#                             b=0){
                             b=0, ass, bss){  # SteveMod
  
  nGroups<-length(YList)
  n<-sapply(YList, length)
  X<- lapply(YList, FUN=function(y) matrix(rep(1, length(y)), ncol=1)) #design matrix for each group. Here just a vector of ones the length of the group size
  nTot<-nkeep+nburn
  
  #storage vectors and matrices
  tau2.sample<-numeric(nTot)
  mu.sample<-numeric(nTot)
  theta.sample<-matrix(NA, nrow=nTot, ncol=nGroups)
  sigma2.sample<-matrix(NA, nrow=nTot, ncol=nGroups)
  yAccept<-matrix(NA,nTot,nGroups)
  logprop.curr<-matrix(NA,nTot,nGroups)
  yMhRatios<-matrix(NA,nTot,nGroups)
  
  
  #
  # define the psi function
  #
  if(regEst=='Huber') {
    psi<-get('psi.huber') #internal
    fn.psi<-get('fn.psi.huber')
    
  } else { 
    if(regEst=='Tukey'){
      psi<-get('psi.bisquare') #internal
      fn.psi<-get('fn.psi.bisquare')
    } else {stop("only set up for Huber or Tukey regression estimates")}}
  
  if(scaleEst!='Huber'){
    stop('scale estimate only set up for Hubers Prop2')
  }
  fn.chi<-fn.chi.prop2
  
  robustList<-list();length(robustList)<-nGroups
  thetaHatObsList<-list();length(thetaHatObsList)<-nGroups
  sigHatObsList<-list();length(sigHatObsList)<-nGroups
  
  #fit the robust linear model to each group separately
  #These are the statistics being conditioned on.
  for(groups in 1:nGroups){
    robust<-rlm(X[[groups]],YList[[groups]],psi=psi, scale.est=scaleEst, maxit=maxit)
    robustList[[groups]]<-robust
    thetaHatObsList[[groups]]<-robust$coef
    sigHatObsList[[groups]]<-robust$s
  }
  
 #
 # set up the projection matrices
 #    
  projList=NULL #projection onto deviation space for each group
  
  if(is.null(projList)) {
                        projList<-list(); length(projList)<-nGroups
                         QtList<-list(); length(QtList)<-nGroups
                        for(i in 1:nGroups){
                           Q<-qr.Q(qr(X[[i]])) #columns of Q form an orthonormal basis for the column space of X. In this case it is just the one vector normalized to length 1. (maybe in the opposite direction...i.e. negative)
                           #tcrossprod(Q,Q) is XX^T in the paper, the projection onto the column space of X.
                           projList[[i]]<-diag(n[i])-tcrossprod(Q,Q)
                           #projList these are the WW^T for each group in the paper
                           QtList[[i]]<-t(Q) #save for computation. These are the U^T's is the paper (where U forms orthonormal basis for column space of X)
                         }
  }
  
 #
 #choose starting value for yi
 #
  log.prop.den.curr<-numeric(nGroups)
  yCurr<-list(); length(yCurr)<-nGroups
  for(i in 1:nGroups){
    yCurr[[i]]<-numeric(n[i])
  }
  
  for(i in 1:nGroups){
    y.prop <- rnorm(n[i])
    Xi<-X[[i]]
    thetaHatObsi<-thetaHatObsList[[i]] #observed stats in each group
    sigHatObsi<-sigHatObsList[[i]]
    proji<-projList[[i]]
    Qti<-QtList[[i]]
    yCurri <- fn.comp.ystst(y.prop,Xi,l1obs=thetaHatObsi,s1obs=sigHatObsi,psi,scale.est=scaleEst,maxit)
    yCurri<-yCurri
    yCurr[[i]]<-yCurri
    log.prop.den.curr[i]<-log.prop.den2(yCurri,Xi, proji,Qti,thetaHatObsi, sigHatObsi,fn.psi, fn.chi,n[i],p=ncol(Xi))
  }
  
# Begin MCMC 
  for(it in 1:nTot){
    
    #sample from complete model
#    samp<-fn.one.rep.complete.MCMC(yCurr, theta, sigma2, mu, tau2,a=a,b=b, ass=ass, bss=bss)
    samp<-fn.one.rep.complete.MCMC(yCurr, theta, sigma2, mu, tau2,a=a,b=b, ass=ass, bss=bss) # SteveMod
    theta<-samp$theta
    sigma2<-samp$sigma2
    mu<-samp$mu
    tau2<-samp$tau2
    
  
    # [y|rest + robust statistics] step; each group updated separately
    for(gp in 1:nGroups){
      yicurr<-yCurr[[gp]]
      Xi<-X[[gp]]
      thetaCuri<-theta[gp]
      sigma2Curi<-sigma2[gp]
      thetaHatObsi<-thetaHatObsList[[gp]]
      sigHatObsi<-sigHatObsList[[gp]]
      log.prop.den.curri<-log.prop.den.curr[gp]
      proj<-projList[[gp]]
      Qt<-QtList[[gp]]
      ySample<-fn.one.rep.y2(yicurr,thetaCuri,sqrt(sigma2Curi),thetaHatObsi, sigHatObsi,Xi, log.prop.den.curri, proj,Qt,fn.psi,fn.chi, psi,scaleEst,maxit)
      yCurr[[gp]]<-ySample[[1]]
      yAccept[it,gp]<-ySample[[2]]
      log.prop.den.curr[gp]<-ySample[[3]]
      logprop.curr[it,gp]<-ySample[[3]]
      yMhRatios[it,gp]<-ySample[[4]]
    }
  
    
    theta.sample[it,]<-theta
    sigma2.sample[it,]<-sigma2
    mu.sample[it]<-mu
    tau2.sample[it]<-tau2
    if(it %% printEach==0) print(paste('iteration', it, sep=' '))
    
  }    
    
    out<-list()
    keep<-(nburn+1):nTot
    out$theta<-theta.sample[keep,]
    out$sigma2<-sigma2.sample[keep,]
    out$mu<-mu.sample[keep]
    out$tau2<-tau2.sample[keep]
    out$robustFits<-robustList
    out$yAccept<-yAccept
    out$logprop.curr<-logprop.curr
    out$yMhRatios<-yMhRatios
    out  
}



#
# Data for the Simulation----
#

#Function generates the simulation data and saves the factors and true thetas

fn.gen.data<-function(sigma2){
  p<-c(.1,.2,.3)
  n<-c(25,50, 100)
  m<-c(9, 25)
  factors<-expand.grid(p,n,m)
  yList<-list()
  factorsList<-list()
  theta<-numeric(nrow(factors)*5)
  for(nfac in 1:nrow(factors)){
    p<-factors[nfac,1]
    n<-factors[nfac,2]
    m<-factors[nfac,3]
    for(j in 1:5){
      thetai<-rnorm(1) #muTrue=0, tau2True=1
      theta[5*(nfac-1)+j]<-thetai
      ps<-rbinom(n,1,p)
      vars<-ifelse(ps, m, 1)*sigma2
      yList[[5*(nfac-1)+j]]<-rnorm(n, thetai, sqrt(vars))
      factorsList[[5*(nfac-1)+j]]<-c(p,n,m)
  }
  }
  out<-list()
  out$factorsList<-factorsList
  out$yList<-yList
  out$theta<-theta
  out
}


# tst<-fn.gen.data(1)
# length(tst$yList[[31]])
# tst$factorsList
# set.seed(1)
# tst1<-fn.gen.data(1)
# 
# set.seed(1)
# tst2<-fn.gen.data(10)
# 
# all.equal(tst1$yList[[1]], tst2$yList[[1]])
# head(tst1$yList[[1]])
# head(tst2$yList[[1]])
# 
#
#
#
#not used----
#
#
#
#
# # Priors ----
# 
# # tau^2 prior ----
# 
# fn.prior.tau2<-function(tau2, log=TRUE){
#   den<-(1/tau2)*(tau2>0)
#   if(log) {
#     log(den)
#   } else{
#     den
#   }
# }
# 
# # mu prior ----
# 
# fn.prior.mu<-function(mu, log=TRUE){
#   den<-mu>0
#   if(log) {
#     log(den)
#   } else{
#     den
#   }
# }
# 
# 
# # theta_i prior ----
# 
# fn.prior.theta.i<-function(thetai, mu, tau2, log=TRUE){
#   dnorm(thetai, mu, sqrt(tau2), log=log)
# }
# 
# 
# 
# # sigma2_i prior ----
# fn.prior.sigma2.i<-function(sigma2i, log=TRUE){
#   den<-(1/sigma2i)*(sigma2i>0)
#   if(log) {
#     log(den)
#   } else{
#     den
#   }
# }
