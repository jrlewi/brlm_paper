# 
#one.rep  function and overall function
# ---------------------------


source("nwHierNormalModelSamplingFunctions.R")

#function for one time through gibbs sampler
#all tunning constants should be defined outside of function
fn.hier.one.rep<-function(y,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, #first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho
                          ){ #y, X are list of group level responses
  #fn.one.rep.beta.l loops through for each beta.l
  #[beta_i|-]
  sigma2<-fn.compute.sigma2(Z,a0,b0)
  betalMat<-fn.one.rep.beta.l(y, X,XtX,Beta,sigma2,bstar,betalMat)
  #[Z|-] i.e. the sigma2_is
  #fn.one.rep.Z loops through for all the z_i
  Z<-fn.one.rep.Z(Z, y,X, betalMat, rho,step=step_Z)
  sigma2<-fn.compute.sigma2(Z,a0,b0)
  #[Beta|-]
  Beta<-fn.sample.Beta(betalMat,bstar)
  
 # [mu_bstr|-]
 # mu_bstr<-fn.sample.mubstr(mu_bstr,psi_bstr, bstar,step_mubstr) 
  
  #[psi_bstr|-]
#  psi_bstr<-fn.sample.psibstr(psi_bstr,mu_bstr, bstar)
  
  #[bstar|-]
  quad1<-fn.compute.quad1(Beta, mu0)
  quad2<-fn.compute.quad2(Beta, betalMat)
  bstar<-fn.sample.bstar(bstar,v1,v2,quad1,quad2,K=nGroups,p) 
  #[mu_rho|-]
  mu_rho<-fn.sample.mu_rho(mu_rho,psi_rho, rho)
  
  #[psi_rho|-]
  psi_rho<-fn.sample.psi_rho(psi_rho,mu_rho, rho)
  
  #[rho|-]
  rho<-fn.sample.rho(rho, mu_rho, psi_rho, Z)
  
  out<-list()
  out$Beta<-Beta
  out$betalMat<-betalMat
  out$Z<-Z
  out$sigma2<-sigma2
  #out$mu_bstr<-mu_bstr
  #out$psi_bstr<-psi_bstr
  out$bstar<-bstar
  out$mu_rho<-mu_rho
  out$psi_rho<-psi_rho
  out$rho<-rho
  out
}
#fitting function
#y is a list of the responses
#X is the design Matrix-a list of the design matrices Xi
#Xi is the design matrix for each group
# sigma2Int is the vector of sigma2i initial values
#fixed hyperparameters:
#mu0: prior mean of each beta
#Sigma0 is the 'variance' matrix of beta b_i~N(mu0,b^*Sigma0)
#a0, b0 prior parameters for sigma2
#mu_bstr (originally had a prior)--mean of the beta prior for b^*
#psi_bstr (originally had a prior)--precision of the beta prior for b^*
      #in detail b*~beta(v1,v2), mu_bstr=v1/(v1+v2), psi_bstr=v1+v2
#swSq=1: keep as one
#w1, w2, a_psir, b_psir: parameters definining prior for rho
  #in detail: 
  #rho~beta(mean=mu_rho, precision=psi_rho)
  #mu_rho~beta(w1,w2) where mu_rho is the mean of rho.
  #psi_rho~gamma(a_psir, b_psir)

heirNormTheoryLm<-function(y,
                           X,
                           nkeep=1e4, nburn=1e3,
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
                           b_psir
                           ){
  mu0<<-mu0
  Sigma0<<-Sigma0
  a0<<-a0
  b0<<-b0
  #alpha_mustr<<-alpha_mustr
  #beta_mustr<<-beta_mustr
  #a_psib<<-a_psib
  #b_psib<<-b_psib
  v12<<-fn.compute.ab(mu_bstr,psi_bstr)
  v1<<-v12[1]
  v2<<-v12[2]
  swSq<<-swSq
  w1<<-w1; w2<<-w2
  a_psir<<-a_psir
  b_psir<<-b_psir
  ############
  XtX<-lapply(X, FUN=function(X) t(X)%*%X)
  p<<-length(mu0) #the number of reg. coefficients per group
  pTot<<-length(X)*p
 #n<<-length(unlist(y))
  nGroups<<-pTot/p
  Sigma0Inv<<-solve(Sigma0)
  total<-nkeep+nburn
#initialize outputs
  #list of betaSamples: each element of list is betaSamples from corresponding group
  betaGroupSamples<-array(NA, c(p,nGroups, total))
  BetaSamples<-matrix(NA,total,p)
  sigma2Samples<-matrix(NA,total,nGroups)
  Zsamples<-matrix(NA,total,nGroups)
  #mu_bstrSamples<-numeric(total)
  #psi_bstrSamples<-numeric(total)
  bstarSamples<-numeric(total)
  mu_rhoSamples<-numeric(total)
  psi_rhoSamples<-numeric(total)
  rhoSamples<-numeric(total)
  tstbst<-matrix(NA, total, 4)


  #initial values
  #starting values of parameters
  #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
  #psi_bstr<-rgamma(1, a_psib,b_psib)  
  #v<-fn.compute.ab(mu_bstr,psi_bstr)
  bstar<-mu_bstr #rbeta(1, v[1],v[2]) qbeta(.5,v[1],v[2]) #
  Beta<-mvrnorm(1,mu=mu0, Sigma=(1-bstar)*Sigma0)
  betalMat<-matrix(Beta, p,nGroups)
  mu_rho<-rbeta(1,w1,w2)
  psi_rho<-rgamma(1, a_psir,b_psir)
  ab<-fn.compute.ab(mu_rho,psi_rho)
#   rho<-rbeta(1,ab[1],ab[2])
  rho<-runif(1, .1,.9)
  Sigma_rho<-(1-rho)*diag(nGroups)+matrix(1, nGroups,nGroups)*rho
  
  Z<-mvrnorm(1, rep(0,nGroups), Sigma_rho)
  
  
  
  for(iter in 1:total){
    samp<-fn.hier.one.rep(y,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, #first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho)
    #update temp values
    Beta<-samp$Beta
    betalMat<-samp$betalMat
    Z<-samp$Z
    sigma2<-samp$sigma2
    #mu_bstr<-samp$mu_bstr
    #psi_bstr<-samp$psi_bstr
    bstar<-samp$bstar
    mu_rho<-samp$mu_rho
    psi_rho<-samp$psi_rho
    rho<-samp$rho
    #update outputs
    BetaSamples[iter,]<-Beta
    betaGroupSamples[,,iter]<-betalMat
    sigma2Samples[iter,]<-sigma2
    Zsamples[iter,]<-Z
    #mu_bstrSamples[iter]<-mu_bstr
    #psi_bstrSamples[iter]<-psi_bstr
    bstarSamples[iter]<-bstar
    mu_rhoSamples[iter]<-mu_rho
    psi_rhoSamples[iter]<-psi_rho
    rhoSamples[iter]<-rho
    } 
  out<-list()
  out$Beta<-BetaSamples[-c(1:nburn),]
  out$betal<-betaGroupSamples[,,-c(1:nburn)]
  out$sigma2s<-sigma2Samples[-c(1:nburn),]
  out$Z<-Zsamples[-c(1:nburn),]
  #out$mu_bstr<-mu_bstrSamples[-c(1:nburn)]
  #out$psi_bstr<-psi_bstrSamples[-c(1:nburn)]
  out$bstar<-bstarSamples[-c(1:nburn)]
  out$mu_rho<-mu_rhoSamples[-c(1:nburn)]
  out$psi_rho<-psi_rhoSamples[-c(1:nburn)]
  out$rho<-rhoSamples[-c(1:nburn)]
  hypers<-c(a0,b0,swSq,w1,w2,a_psir,b_psir, mu_bstr,psi_bstr)
names(hypers)<-c("a0","b0","swSq",'w1',"w2","a_psir","b_psir", "mu_bstr","psi_bstr")
  out$hypers<-hypers
  out
} 





