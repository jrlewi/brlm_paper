#############################
#Fitting Hierarchical Bayes Basic Normal Theory Linear Model
#To the Nationwide data set
#############################
library(MCMCpack)
library(msm) #for truncated normal
library(Matrix)
library(MASS)
library(plyr)
#######
#Model
#######

#------------------------------
#Code to fit the following Hierarchical linear model to the 
#NW data set
#------------------------------

#[beta| mu0, Sigma0, b^*]~N(mu0, a*Sigma0) #mu0, Sigma0- fixed
#[b_i|beta, Sigma0, b^*]~iid N(beta,b^*Sigma0) 

#b=b^*swSq, where swSq is fixed (swSq=sum w_i^2, from previous data set)

#Variance is partitioned between beta and (b_i, i=1,..,K) to connect to the non-hierarchical version of this model
# this reasoning implies a+b^*=1

#[b^*|v1, v2]~beta(v1, v2) 

#priors are placed on mu_b*=v1/(v1+v2) and psi_b*=v1+v2 in the following way
#mu_b*~beta(alpha_mu*, beta_mu*) with parameters fixed so the mean is small....see notes and setting of prior parameters
#psi_b*~gamma(a_psib, b_psi_b) 
#that is, the mean of b^*|v1,v2 is mu_b* and the variance is mu_b*(1-mu_b*)/(psi_b*+1)


#[sigma^2_i |mu_rho, K_rho, rho]~InG(a0,b0) a0, b0 fixed from previous data set (use the same values as the non-hiearchical fit). But sigma^2_i's not independent

# connect the sigma^2_i's via the formulation
# Z~N_K(0, Sigma_rho) with sigma^2_i=H^-1(phi(z_i)) where H is the cdf of IG(a0, b0) and Sigma_rho=rho*J+(1-rho)I

#to allow for data dependence in the correlation of the sigma^2_i

#[rho| mu_rho, psi_rho]~beta(alpha_rho, beta_rho) where alpha_rho and beta_rho are functions of mu_rho, psi_rho, so that the mean of the beta for rho is mu_rho and the variance is mu_rho(1-mu_rho)/(psi_rho+1)

#[mu_rho] ~beta(w1,w2), [psi_rho]~Gamma(a_psir, b_psir), fixed parameters


#finally, for the K groups i=1,...,K
# y_i=b_i*x_ij+e_ij j=1,..., n_i (n_i is the number of obs in the ith group)

# y_ij|b_i,sigma_i^2~N(b_i*x_ij,sigma_i^2)

# e_ij~(iid) N(0, sigma_i^2) 

# mu_0 and Sigma0 are fixed a priori, so are a0 and b0
# NOTE: y and X are lists of the individual responses and design matrices for each group
#Model: obtain samples from [mu_b*, psi_b*, b*,beta, beta_i, i=1,...,K, Z,  mu_rho, psi_rho,rho|Y]
#strategy: write sampling functions for each full conditional, then combine into one function for one repetition
#fixed HyperParameters are
#mu0, Sigma0,Sigma0Inv, a0,b0,alpha_mustr, beta_mustr,a_psib, b_psib
#w1, w2,a_psir, b_psir
# other fixed values are swSq (sum of squared weights), p, nGroups, J_nGroups (nGroups by nGroups matrix of ones), I_nGroups (nGroups by nGroups Identity)


#fix mu_bstr and psi_bstr

# #--------------
# #1. [mu_bstr|-]
# #function to compute shape parameters of beta from the mean mu<-a/(a+b) and psi=a+b (precision)
fn.compute.ab<-function(mu, psi){
  a<-mu*psi
  b<-psi-a
  c(a,b)
}
# 
# fn.loglike.mubstr<-function(mu_bstr, psi_bstr, bstar){
#   if(mu_bstr>0 && mu_bstr<1){
#   ab<-fn.compute.ab(mu_bstr, psi_bstr)
#   v1<-ab[1]
#   v2<-ab[2]
#   -log(beta(v1,v2))+(v1-1)*log(bstar)+(v2-1)*log(1-bstar)+(alpha_mustr-1)*log(mu_bstr)+(beta_mustr-1)*log(1-mu_bstr)
# } else {
#   -Inf
# }
# }
# 
# 
# fn.proposal.mubstr<-function(mubstr, step_mubstr){ #alpha_mustr,beta_mustr
#   rnorm(1, mubstr, step_mubstr)
# }
# 
# # 
# fn.proposal.mubstr<-function(mubstrCur, step_mubstr){ #alpha_mustr,beta_mustr
#   logmub<-rnorm(1, log(mubstrCur), step_mubstr)
#   exp(logmub)
# }
# 
# 
# fn.logpdf.prop.mubstr<-function(mu_bstrProp){
#  -log(mu_bstrProp)
# } #
# 
# 
# # fn.logpdf.prop.mubstr<-function(mu_bstrProp,mu_bstrCur,step_mubstr){
# #   dnorm(mu_bstrProp,mu_bstrProp, step_mubstr, log=TRUE)
# #   } #
# 
# fn.sample.mubstr<-function(mu_bstrCur,psi_bstr, bstar,step_mubstr) {
#   mu_bstrProp<-fn.proposal.mubstr(mu_bstrCur, step_mubstr)
#   logRat1<-fn.loglike.mubstr(mu_bstrProp, psi_bstr, bstar)-fn.loglike.mubstr(mu_bstrCur, psi_bstr, bstar)
#   logRat2<-fn.logpdf.prop.mubstr(mu_bstrCur)-fn.logpdf.prop.mubstr(mu_bstrProp)
#   mhRat4<-min(exp(logRat1+logRat2),1)
#   if(runif(1)<mhRat4){
#     return(mu_bstrProp)
#   } else {
#     return(mu_bstrCur)
#   } 
# }
# 
# #-------------------------
# 
# #--------------
# #2.[psi_bstr|--] a_psib, b_psib
# fn.loglike.psibstr<-function(psi_bstr,mu_bstr, bstar){
#   if(psi_bstr>0){
#   ab<-fn.compute.ab(mu_bstr, psi_bstr)
#   v1<-ab[1]; v2<-ab[2]
#   -log(beta(v1,v2))+(v1-1)*log(bstar)+(v2-1)*log(1-bstar)+(a_psib-1)*log(psi_bstr)-b_psib*psi_bstr
# } else {
#   -Inf
# }
#   }
# 
# #proposing from prior like this will cancel part of the log.like; but include it for now for easy debug
# # fn.proposal.psibstr<-function(){ #a_psib, b_psib
# #   rgamma(1,a_psib, rate=b_psib)
# # }
# 
# fn.proposal.psibstr<-function(psi_bstrCur, psi_bstr_step){ #a_psib, b_psib
#   rnorm(1, psi_bstrCur,psi_bstr_step)
# }
# 
# 
# # fn.logpdf.prop.psibstr<-function(psi_bstr){
# #   (a_psib-1)*log(psi_bstr)-b_psib*psi_bstr
# # } #
# 
# fn.sample.psibstr<-function(psi_bstrCur,mu_bstr, bstar) {
#   #psi_bstrProp<-fn.proposal.psibstr()
#   psi_bstrProp<-fn.proposal.psibstr(psi_bstrCur, psi_bstr_step)
#   logRat1<-fn.loglike.psibstr(psi_bstrProp,mu_bstr, bstar)-fn.loglike.psibstr(psi_bstrCur,mu_bstr, bstar)
#   logRat2<-0 #fn.logpdf.prop.psibstr(psi_bstrCur)-fn.logpdf.prop.psibstr(psi_bstrProp)
#   mhRat5<-min(exp(logRat1+logRat2),1)
#   if(runif(1)<mhRat5){
#     return(psi_bstrProp)
#   } else {
#     return(psi_bstrCur)
#   } 
# }
# 

  
# ---------------------------
#3. [b^*|-]
#involves MH step
#quad1 is the current value of (beta-mu0)'Sigma0Inv(beta-mu0)
#quad 2 is the current value of sum_i^nGroups (beta_i-beta) Sigma0Inv (beta_i-beta)

fn.compute.quad1<-function(Beta,mu0){
  t(Beta-mu0)%*%Sigma0Inv%*%(Beta-mu0)
}
fn.compute.quad2<-function(Beta, betalMat){
  diff<-betalMat-Beta
  sum(apply(diff,MARGIN=2,FUN=function(x) t(x)%*%Sigma0Inv%*%x))
}

#with mu_bstr, psi_bstr fixed....replace these imputs with v1 and v2
fn.loglike.bstar<-function(bstar,v1,v2, quad1, quad2, K,p)
{
  if(bstar>0 && bstar<1){
    #ab<-fn.compute.ab(mu_bstr, psi_bstr)
    #v1<-ab[1]; v2<-ab[2]
    (v1-.5*K*p-1)*log(bstar)+(v2-.5*p-1)*log(1-bstar)-.5*quad1/(1-bstar)-.5*swSq*quad2/bstar #-log(beta(v1,v2))
    } else { 
         -Inf}
}


# 
# fn.proposal.bstar<-function(bstarCur, step_bstar){
#   logb<-rnorm(1, log(bstarCur),step_bstar)# sd=step)
#   exp(logb)
# }


# fn.proposal.bstar<-function(bstarCur, step_bstar){
#  rnorm(1, bstarCur,step_bstar)# sd=step)
# }
#for proposing on log scale
fn.proposal.logbstar<-function(logbstarCur, step_logbstar){
  rnorm(1, logbstarCur,step_logbstar)# sd=step)
}

fn.logpdf.prop.logbstar<-function(bstarProp){
  -log(bstarProp)
} #dnorm(log(bstarProp),log(bstarCur),step_logbstar, log=TRUE)-log(bstarProp)



# fn.sample.bstar<-function(bstarCur,mu_bstr,psi_bstr,quad1,quad2, K) {
#   bstarProp<-fn.proposal.bstar(bstarCur, step_bstar)
#   logRat1<-fn.loglike.bstar(bstarProp,mu_bstr,psi_bstr, quad1, quad2, K)-fn.loglike.bstar(bstarCur,mu_bstr,psi_bstr, quad1, quad2,K)
#   logRat2<-0 #fn.logpdf.prop.bstar(bstarCur,bstarProp, step_bstar)-fn.logpdf.prop.bstar(bstarProp,bstarCur, step_bstar)
#   mhRat6<-min(exp(logRat1+logRat2),1)
#   if(runif(1)<mhRat6){
#     return(bstarProp)
#   } else {
#     return(bstarCur)
#   } 
# }
# 
#for proposing on log scale
fn.sample.bstar<-function(bstarCur,v1,v2,quad1,quad2, K,p) {
  logbstarCur<-log(bstarCur)
  logbstarProp<-fn.proposal.logbstar(logbstarCur, step_logbstar)
  bstarProp<-exp(logbstarProp)
  logRat1<-fn.loglike.bstar(bstarProp, v1,v2, quad1, quad2, K,p)-fn.loglike.bstar(bstarCur, v1,v2, quad1, quad2,K,p)
  logRat2<-fn.logpdf.prop.logbstar(bstarCur)-fn.logpdf.prop.logbstar(bstarProp)
  mhRat6<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat6){
    return(bstarProp)
  } else {
    return(bstarCur)
  } 
}








# ---------------------------

# ---------------------------

#4. [Beta|-]
#betalMat, matrix with columns equal to current beta_l's from each group
#bstar: current value of bstar
#Sigma0-fixed 
#nGroups-fixed; this is the number of groups K

fn.sample.Beta<-function(betalMat,bstar){
  a<-1-bstar 
  b<-bstar/swSq
  prec<-(1/a+nGroups/b)
  SigmaB<-(Sigma0)/prec
  muB<-(rowSums(betalMat)/b+mu0/a)/prec
  mvrnorm(1, muB, SigmaB)
}

# ---------------------------

# ---------------------------


#loop through for each group
#5. [beta_l|-] l=1,..,nGroups
#leave yl and Xl as inputs for easy transition to restricted model.
# these are the responses and design matrix for the lth group
#XtX.l t(X)%*%X for lth group
#sigma2l: the sigma^2_l for the l group
#b=b^*/swSq
#Sigma0Inv: fixed
fn.sample.beta.l<-function(yl, Xl,XtX.l,Beta,sigma2l,bstar){
  b<-bstar/swSq
  Sigma.betal<-solve(XtX.l/(sigma2l)+Sigma0Inv/b)
  mu.l<-Sigma.betal%*%(crossprod(Xl,yl)/(sigma2l)+Sigma0Inv%*%Beta/b)
  mvrnorm(1, mu.l, Sigma.betal)
}

fn.one.rep.beta.l<-function(y, X,XtX,Beta,sigma2,bstar,betalMat){
  for(l in 1:ncol(betalMat)){
    yl<-y[[l]]
    Xl<-X[[l]]
    XtX.l<-XtX[[l]]
    sigma2l<-sigma2[l]
    beta.l<-fn.sample.beta.l(yl, Xl,XtX.l,Beta,sigma2l,bstar)
    betalMat[,l]<-beta.l
  }
betalMat
}

#---------------------------------------------

#---------------------------------------------
#[mu_rho|-]

fn.loglik.mu_rho<-function(mu_rho,psi_rho, rho){
#  if(mu_rho>0.001 && mu_rho<.999){
#  if((mu_rho>0.001 && mu_rho<.999) && psi_rho>0){
  if(mu_rho>0 && mu_rho<1){
    ab<-fn.compute.ab(mu_rho,psi_rho)
    alpha_rho<-ab[1]
    beta_rho<-ab[2]
    lik<-(w1-1)*log(mu_rho)+(w2-1)*log(1-mu_rho)+(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-log(beta(alpha_rho,beta_rho))     
  } else {
    lik<--Inf
  }
  return(lik)
}

# fn.proposal.mu_rho<-function(){
#   runif(1)
#   # rbeta(1, w1,w2)
#   #rnorm(1,mu_rhoCur, mu_rho_step)
# }

fn.proposal.mu_rho<-function(mu_rhoCur, mu_rho_step){
  rnorm(1,mu_rhoCur, mu_rho_step)
}

# fn.proposal.mu_rho<-function(mu_rhoCur, mu_rho_step){
#   logmur<-rnorm(1,log(mu_rhoCur), mu_rho_step)
#   exp(logmur)
# }
# 
# fn.logpdf.prop.mu_rho<-function(mu_rhoProp, mu_rhoCur,mu_rho_step){
#   dnorm(log(mu_rhoProp),log(mu_rhoCur),mu_rho_step, log=TRUE)-log(mu_rhoProp)
# }


# fn.logpdf.prop.mu_rho<-function(mu_rho){
#   (w1-1)*log(mu_rho)+(w2-1)*log(1-mu_rho)
# }

# alpha_rho and beta_rho used in multiple steps
fn.sample.mu_rho<-function(mu_rhoCur,psi_rho, rho){
  mu_rhoProp<-fn.proposal.mu_rho(mu_rhoCur, mu_rho_step)
  #mu_rhoProp<-fn.proposal.mu_rho()
  logRat1<-fn.loglik.mu_rho(mu_rhoProp,psi_rho, rho)-fn.loglik.mu_rho(mu_rhoCur,psi_rho, rho)
  logRat2<-0 #fn.logpdf.prop.mu_rho(mu_rhoCur,mu_rhoProp,mu_rho_step)-fn.logpdf.prop.mu_rho(mu_rhoProp, mu_rhoCur,mu_rho_step)
  mhRat7<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat7){
    return(mu_rhoProp)
  } else {
    return(mu_rhoCur)
  } 
}
#---------------------------------------------

#---------------------------------------------
# [psi_rho|-]

fn.loglik.psi_rho<-function(psi_rho,mu_rho, rho){
     #if((mu_rho>0 && mu_rho<1) && psi_rho>0){(mu_rho>0.001 && mu_rho<.999) &&
  if(psi_rho>0){
    ab<-fn.compute.ab(mu_rho,psi_rho)
    alpha_rho<-ab[1]
    beta_rho<-ab[2]
   lik<-(a_psir-1)*log(psi_rho)-b_psir*psi_rho+(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-log(beta(alpha_rho,beta_rho))        
  } else {
    lik<--Inf
  }
return(lik)
}

# fn.proposal.psirho<-function(){ #a_psir, b_psir
#   rgamma(1,a_psir, rate=b_psir)
# }

fn.proposal.psirho<-function(psi_rhoCur, psi_rho_step){ #a_psir, b_psir
  rnorm(1,psi_rhoCur, psi_rho_step)
}



fn.logpdf.prop.psirho<-function(psi_rho){
  (a_psir-1)*log(psi_rho)-b_psir*psi_rho
} #

#fn.pdf.prop.mu_rho<-function()
# alpha_rho and beta_rho used in multiple steps....computer these once
fn.sample.psi_rho<-function(psi_rhoCur,mu_rho, rho){
  #psi_rhoProp<-fn.proposal.psirho()
  psi_rhoProp<-fn.proposal.psirho(psi_rhoCur, psi_rho_step)
  #print(psi_rhoProp)
   #print(c(mu_rho,rho))
  logRat1<-fn.loglik.psi_rho(psi_rhoProp,mu_rho, rho)-fn.loglik.psi_rho(psi_rhoCur,mu_rho, rho)
  logRat2<-0 #fn.logpdf.prop.psirho(psi_rhoCur)-fn.logpdf.prop.psirho(psi_rhoProp)
  mhRat1<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat1){
    return(psi_rhoProp)
  } else {
    return(psi_rhoCur)
  } 
}

#---------------------------------------------
#[rho|-]


#J_nGroups must be defined
#I_nGroups must be defined

#function to effeciently caclulate Sigma_rhoInv using its special form and a well known matrix inverse result (from the fix-rank kriging chapter 8 of the spatial handbook book); note; the K variable added for the sampling of the z_i's later
fn.compute.SigmaRhoInv<-function(rho, K){
  #nGroups<-nrow(I_nGroups)
  diag(K)/(1-rho)-matrix(1, K,K)/((1-rho)^2*(K/(1-rho)+1/rho))
}

# 
# fn.compute.SigmaRhoInv<-function(rho, K=nGroups){
#   #nGroups<-nrow(I_nGroups)
#   (diag(K)-matrix(1, K,K)*(rho/(1+(K-1)*rho)))/(1-rho)
# }
# 


# fn.compute.SigmaRhoInv<-function(rho, K=nGroups){
#   #nGroups<-nrow(I_nGroups)
#  solve(rho*matrix(1, K,K)+(1-rho)*diag(K))
# }


# I_nGroups<-diag(10)
# J_nGroups<-matrix(1, 10,10)
# rho<-.4
#  Sigma_rho<-(1-rho)*I_nGroups +rho*J_nGroups
# # 
# all.equal(Sigma_rho%*%fn.compute.SigmaRhoInv(rho, 10), I_nGroups)
# all.equal(fn.compute.SigmaRhoInv(rho,10)%*%Sigma_rho, I_nGroups)


#function to effeciently caclulate log(det(Sigma_rho)) using its special form and Sylvester's determinant theorem on Wiki
#det(I_n+AB)=det(I_m+BA)
fn.compute.logDetSigmaRho<-function(rho, K){
  if(rho>0 & rho <1){
  K*log(1-rho)+log(1+K*rho/(1-rho))} else {
    if(rho==0) {log(1)}
    else {if(rho==1){log(0)
    }
    }
}
}


# fn.compute.logDetSigmaRho<-function(rho, K=nGroups){
#   log(det(rho*matrix(1, K,K)+(1-rho)*diag(K)))
# }

#all.equal(log(det(Sigma_rho)),fn.compute.logDetSigmaRho(rho, 10))

# curve(exp(fn.compute.logDetSigmaRho(x)))

#note: this assumes the proposal for pho is beta(mu_rho, K_rho)


#alpha_rho, beta_rho computed within the sampling function
# fn.proposal.rho<-function(alpha_rho, beta_rho){
#   #rbeta(1, alpha_rho, beta_rho)#, alpha_rho, beta_rho)
#   #rtnorm(1, mean=rhoCur, sd=rho_step, lower=0, upper=1)
#   runif(1)
# }

# fn.proposal.rho<-function(){
#   runif(1)
# }

# fn.proposal.rho<-function(rhoCur, rho_step){
#   rtnorm(1, mean=rhoCur, sd=rho_step, lower=0, upper=1)
# }
# log.pdf.prop.rho<-function(rhoProp,rhoCur,rho_step){
# dtnorm(rhoProp, rhoCur, rho_step, lower=0, upper=1,log=TRUE)
# }

fn.proposal.rho<-function(rhoCur, rho_step){
  rnorm(1, rhoCur, rho_step)
}



# fn.proposal.rho<-function(rhoCur, rho_step){
#   logr<-rnorm(1,log(rhoCur), rho_step)
#   exp(logr)
# }
# 
# fn.logpdf.prop.rho<-function(rhoProp, rhoCur,rho_step){
#   dnorm(log(rhoProp),log(rhoCur),rho_step, log=TRUE)-log(rhoProp)
# }



# log.pdf.prop.rho<-function(rho,alpha_rho,beta_rho){
# (alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)
#  # dtnorm(rhoProp, rhoCur, rho_step, lower=0, upper=1,log=TRUE)
# }

#alpha_rho and beta_rho computed inside sampling function
fn.loglike.rho<-function(rho, Z,alpha_rho,beta_rho){
  if(rho<1 && rho>0){
   #if(rho<0.999 && rho>0.001){
    #Sigma_rho<-rho*J_nGroups+(1-rho)*I_nGroups #do i need this?
    K<-length(Z)
    logdetSigma_rho<-fn.compute.logDetSigmaRho(rho,K)
    Sigma_rhoInv<-fn.compute.SigmaRhoInv(rho,K)
    lik<-(alpha_rho-1)*log(rho)+(beta_rho-1)*log(1-rho)-.5*logdetSigma_rho-.5*t(Z)%*%Sigma_rhoInv%*%Z
    } else {
    lik<--Inf
  }
  lik
}

fn.sample.rho<-function(rhoCur, mu_rho, psi_rho, Z){
  ab<-fn.compute.ab(mu_rho,psi_rho)
  alpha_rho<-ab[1]
  beta_rho<-ab[2]
  #rhoProp<-fn.proposal.rho() 
  #rhoProp<-fn.proposal.rho(alpha_rho, beta_rho) 
  rhoProp<-fn.proposal.rho(rhoCur, rho_step)
 # print(c(rhoCur, rhoProp))
  logRat1<-fn.loglike.rho(rhoProp, Z,alpha_rho,beta_rho)-fn.loglike.rho(rhoCur, Z,alpha_rho,beta_rho)
  logRat2<-0 #fn.logpdf.prop.rho(rhoCur,rhoProp,rho_step)-fn.logpdf.prop.rho(rhoProp,rhoCur,rho_step)
  mhRat2<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat2){
    return(rhoProp)
  } else {
    return(rhoCur)
  } 
}

#---------------------------------------------


# ---------------------------
#
#[Zl|-]: sample one at a time
#function to compute simga2_is from a z vector (or a single z) from nGroups dimensional Normal distribution with marginal means 0 and vars 1
# a0,b0 fixed
fn.compute.sigma2<-function(Z,a0,b0){
  X<-qgamma(pnorm(Z),shape=a0, scale=1/b0) #rate=b0
  1/X
}

# pnorm(0)
# curve(dinvgamma(x,shape=.05,scale=.05), from=0,1 )
# 1/qgamma(.5,shape=.05, scale=1/.05)
# fn.compute.sigma2(0)


fn.compute.Z2<-function(sigma2,a0,b0){
  X<- 1-pgamma(1/sigma2, a0, 1/b0)
  qnorm(X)
}

fn.compute.Z<-function(sigma2,a0,b0){
  X<-integrate(function(x) {dinvgamma(x, shape=a0, scale = b0)},0, sigma2)$value #rate=b0
 qnorm(X)

}
# 
# curve(fn.compute.Z(x, a0, b0), 0, 20,n=1000, ylim=c(-10,5))
# x<-seq(.0000001, 20, length.out=100)
# z1<-numeric(length(x))
# for( i in 1:length(x)){
# z1[i]<-fn.compute.Z(x[i], a0, b0)}
# lines(x, z2, col=4)
# 
# 
# z2<-numeric(length(x))
# for( i in 1:length(x)){
#   z2[i]<-fn.compute.Z2(x[i], a0, b0)}
# lines(x, z2, col=4)
# all.equal(z1,z2)
# plot(z1,z2)


fn.hat<-function(Xl, betal) {Xl%*%betal}
#fn.diffs<-function(yl,ylhat) {yl-ylhat}

#function for the likelihood of one z. 

fn.logLike.zl<-function(zl,zminus, yl, fitsl, condsd,condMeanMult){
  #zl is the compenent I am updating, zminus is all other components; condsd is the condition sd for the conditional normal (zl|zminus): this depends only on current 
  #condMeanMult%*%zminus is the conditional mean for zl|zminus
  condMean<-condMeanMult%*%zminus
  sigma2l<-fn.compute.sigma2(zl,a0,b0)
  dnorm(zl,condMean,condsd,log=TRUE)+sum(dnorm(yl, fitsl, sqrt(sigma2l),log=TRUE))
}



fn.proposal.zl<-function(zlCur, step=1){
  rnorm(1, zlCur,step)
}

fn.logpdf.proposal.Z<-function(zlProp,zlCur,step){
  dnorm(zlProp, zlCur,step,log=TRUE)
  }

fn.sample.Zl<-function(zlCur,zminus,yl, fitsl, condsd,condMeanMult, step=1){
  zlProp<-fn.proposal.zl(zlCur, step)
  logRat1<-fn.logLike.zl(zlProp,zminus, yl, fitsl, condsd,condMeanMult)-fn.logLike.zl(zlCur,zminus, yl, fitsl, condsd,condMeanMult)
  logRat2<-0 #random walk: symmetric
  mhRat3<-min(exp(logRat1+logRat2),1)
  if(runif(1)<mhRat3){
    return(zlProp)
  } else {
    return(zlCur)
  } 
}

#one rep:sampling all of the z's
fn.one.rep.Z<-function(Zcur, y,X, betalMat, rho,step=rep(1,length(Zcur))){
  Sigma_rhoInvMinus<-fn.compute.SigmaRhoInv(rho, K=length(Zcur)-1)    #inverse (1-pho)I+rhoJ for K-1 by K-1 dim matrices
  rho_vec<-rep(rho, length(Zcur)-1)
  condsd<-sqrt(1-t(rho_vec)%*%Sigma_rhoInvMinus%*%rho_vec)
  condMeanMult<-t(rho_vec)%*%Sigma_rhoInvMinus
  for(j in 1:length(Zcur)){
    zlCur<-Zcur[j]
    zminus<-Zcur[-j]
    yl<-y[[j]]
    Xl<-X[[j]]
    betal<-betalMat[,j]
    fitsl<-fn.hat(Xl, betal)
    zlCur<-fn.sample.Zl(zlCur,zminus,yl, fitsl, condsd,condMeanMult,step[j])
    Zcur[j]<-zlCur
  }
  Zcur
} 
  
  








