#---------------
#Functions required to propose a new data vector using MH for a Bayesian linear regression model where conditioning is on estimators of Beta and sigma. 

#Estimators currently supported are M-estimators
#Estimators of Beta: currently supports estimators defined using either Huber's or Tukey's Loss.
#Estimatros of sigma: currently only supports the estimator defined using Huber's Proposal 2

#---------------
#More Details:
#The code will develop the sampling function to sample from the following full conditional in the gibbs step

#[Y|T(Y)=T(y), theta]
#where y is the observed data and the model is

#Y~N(BX, sigma^2I)
#priors
#B~ flexible prior
#sigma^2~ flexible prior
#theta=vector of all parameters including B and sigma^2
#flexible priors on all other parameters in theta

#sampling from [Y|T(y), theta] is independent of priors in the model because all priors cancel in MH ratio. See details in paper.  



#---------------
#Functions needed for the calculations of the gradients of the l
#---------------
library(MASS)

#---------------
#Defining the psi function: 2 choices: Huber's and Tukey's
#---------------

#note: psi.huber is available in the MASS library: I write my own because for deriv=0 MASS's function gives psi(x)/x. I need psi(x). 

fn.psi.huber<-function(u,k=1.345, deriv=0){
  if (!deriv) 
    return(pmax(-k, pmin(k,u)))
  abs(u) <= k
  }

#---------------
#Second option for psi Tukey's bisquare
#---------------
fn.psi.bisquare<-function(u, c = 4.685, deriv = 0)
{
  if (!deriv) 
    return(u*(1 - pmin(1, abs(u/c))^2)^2)
  t <- (u/c)^2
  ifelse(t < 1, (1 - t) * (1 - 5 * t), 0)
}



#---------------
#defining the scale estimate: for now: Huber's proposal 2
#---------------

#b is the b(k) from Huber's book
fn.chi.prop2<-function(u, k=1.345, deriv=0, b=.71016){
  if (!deriv) 
    return(fn.psi.huber(u,k, 0)^2-b)
  2*fn.psi.huber(u,k,1)*fn.psi.huber(u,k)
}


#---------------
#Function to calculate the gradients
#---------------

#the linear equation to solve is of the form A*deriv=b
#function should only input y's that satisfy the conditioning statistic values (l1obs and s1obs)


#for simplicity this is set up only to handle the default tunning constants given by fn.psi.huber, fn.psi.bisquare  and fn.psi.prop2
#k is the tunning constant for psi, k2 is tunning constant for chi;for now this is the same as k

#y: data vector; should satisfy observed statistics l1obs and s1obs
#X: design matrix
#fn.psi: custom psi function defining M-estimator for location returns psi(x) for deriv=0 and psi'(x) for deriv=1  ; default is fn.psi.huber 
#fn.chi: custom chi function defining M-estimator for scale returning chi(x)  for deriv=0 and  chi'(x) for deriv=1  ; default is fn.chi.prop2

#returns: p+1 by n matrix containing the partial derivatives of each statistic wrt the data. First p rows are for the location statistics. Last row is for the scale.
fn.grads<-function(y,X,l1obs, s1obs, fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2){
  p<-ncol(X)
  resids<-y-X%*%l1obs
  stResids<-resids/s1obs
  fn.psi <- get('fn.psi', mode = "function")
  fn.chi<- get('fn.chi', mode = "function")
  psiPrime<-fn.psi(stResids, deriv=1)
  chiPrime<-fn.chi(stResids, deriv=1)
  A<-matrix(NA,p+1,p+1)
  psiColX<-apply(X, 2, FUN=function(x) x*psiPrime)
  A[1:p,1:p]<-s1obs*t(X)%*%psiColX
  A[1:p,p+1]<-t(X)%*%(psiPrime*resids)
  A[p+1,1:p]<-s1obs*t(X)%*%chiPrime
  A[p+1,p+1]<-sum(chiPrime*resids)  
  prime1<-psiColX #for the slopes
  prime2<-c(chiPrime) #for the scale
  rhs<-s1obs*t(cbind(prime1,prime2))
  tryCatch({return(solve(A,rhs))}
           , error=function(err){print(paste('The Error:', err))
                                 return(matrix(rep(0, (p+1)*length(y)), p+1,length(y)))
           }
  )
}


#----------------
#takes the yst and transfroms it to the ystst matching the statistics
#----------------
#yst: data vector to be transformed to ystst satisfying the observed summary statistics
#X design matrix
#l1obs, s1obs: observed statsitics
#psi: psi function returning psi(x)/x or psi'(x). default is psi.huber from MASS
#scale.est: defines the scale estimator. Currently supports only 'Huber'
#maxit: number of iterations to find T(yst)

#returns: ystst- a data vector such that T(ystst)=T(y_obs)
fn.comp.ystst<-function(yst,X,l1obs,s1obs,psi=psi.huber, scale.est='Huber', maxit=400){
  fit.yst<-rlm(X,yst,psi=psi, scale.est=scale.est, maxit=maxit)
  l1yst<-fit.yst$coefficients #(fit.yst)
  s1yst<-fit.yst$s 
  res<-yst-X%*%l1yst
  ystst<-s1obs*res/s1yst+X%*%l1obs
  return(ystst)
}


#---------------
#computes the radius in deviation space
#-------------

#ystst- a data vector such that T(ystst)=T(y_obs)
#X design matrix

#returns: the radius in deviation space
fn.radius<- function(ystst,X) 
{
  Bls<-solve(t(X)%*%X)%*%t(X)%*%ystst
  resid<-ystst-X%*%Bls
  rad <-sqrt(sum(resid^2))
  return(rad)
}


#second, equivalent way to calculate radius using proj matrix
#ystst- a data vector such that T(ystst)=T(y_obs)
#proj: the projection matrix onto the deviation space (aka the least squares residual space)

#returns: the radius in deviation space
fn.radius2<-function(ystst,proj){
  resid<-proj%*%ystst
  rad <-sqrt(sum(resid^2))
  return(rad)
}

#ystst- a data vector such that T(ystst)=T(y_obs)
#X: design matrix
#proj: the projection matrix onto the deviation space (aka the least squares residual space)
#l1obs, s1obs: observed statsitics
#fn.psi: custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber 
#fn.chi: custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2

#Returns: 'Attenuations' of the transformation to ystst. These attenuations make up part of the Jacobian

fn.attenuation <- function(ystst,X, proj,l1obs,s1obs,fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2) 
{
  p<-ncol(X)
  n<-nrow(X)
  #W<-t(Vt)[,(p+1):n]
  grads.ystst<-fn.grads(ystst,X, l1obs, s1obs,fn.psi, fn.chi)
  grad.y.B<-grads.ystst[1:p,]
  grad.y.s<-grads.ystst[(p+1),]
  unit.a<-grad.y.s/sqrt(sum(grad.y.s*grad.y.s)) #normalize
  zstst<-proj%*%ystst
  unit.b<- zstst/sqrt(sum(zstst*zstst))
  atten1 <-abs(sum(unit.a * unit.b))
  qr.A<-qr(cbind(grad.y.s,t(grad.y.B)))
  Qmat<-qr.Q(qr.A)
  tmpPrj<-((diag(n)-proj)%*%Qmat)[,-1]
  tmpPrj2<-t(tmpPrj)%*%tmpPrj 
  atten2<-sqrt(abs(det(tmpPrj2)))
  return(atten1*atten2)
}

#added 4/3/2014
#updated, faster version of fn.attenuation. 
#Additional input is  Qt: Q transpose where Q is the orthonormalized X. That is, the columns of Q form an orthnormal basis for the column space of X


fn.attenuation2 <- function(ystst,X,proj, Qt,l1obs,s1obs,fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2)
{
  #Qt is just the t(Q) in Q%*%t(Q)+W%*%t(W)=I; C(X)=C(Q); Q orthonormalized X
  p<-ncol(X)
  n<-nrow(X)
  grads.ystst<-fn.grads(ystst,X, l1obs, s1obs,fn.psi, fn.chi)
  grad.y.B<-grads.ystst[1:p,]
  grad.y.s<-grads.ystst[(p+1),]
  unit.a<-grad.y.s/sqrt(sum(grad.y.s*grad.y.s)) #normalize
  zstst<-proj%*%ystst
  unit.b<- zstst/sqrt(sum(zstst*zstst))
  atten1 <-abs(sum(unit.a * unit.b))
  if(!is.na(atten1)){
  if(p==1){ #patch to allow this function to work when p=1
    qr.A<-qr(cbind(grad.y.s,grad.y.B))
  } else {  
  qr.A<-qr(cbind(grad.y.s,t(grad.y.B)))
  }
  Qmat<-qr.Q(qr.A)
  tmpPrj<-Qt%*%Qmat
  sings<-svd(tmpPrj)$d
  atten2<-prod(sings)
  return(atten1*atten2)} else {
    return(Inf) #return inf-this will cause log.prop.den2=Inf so the proposed ystst will be rejected
  }
  
}
#edited with if else statement to accomodate the tryCatch added to fn.grad 


#-------------
#Likelihood Function for full data
#------------


#ystst- a data vector such that T(ystst)=T(y_obs)
#X: design matrix
#Beta: regression coeff
#sig: sigma, standard deviation of residuals
#returns: the log likelihood 
log.fn.eval.lik <- function(ystst, X,Beta, sig)
{
  lik <- sum(dnorm(ystst,mean=X%*%Beta,sd=sig, log=TRUE))
  return(lik)
}


#-------------------------
#function to caluclate the proposal dens on the log scale 
#------------------------

#y.prop: proposed data vector satisfying T(y.prop)=T(y_obs)
#X: design matrix
#l1obs, s1obs: observed statsitics
#fn.psi: custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber 
#fn.chi: custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2
#n: length(y_obs)=length(y.prop)
#p: number of regression coefficients

#returns: the proposal density for y.prop on the log scale (assuming original vector sampled uniformly on unit sphere)
log.prop.den <-function(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi,n,p){log(fn.attenuation(y.prop,X,proj, l1obs,s1obs,fn.psi, fn.chi))-(n-p-1)*log(fn.radius(y.prop,X))}



#Function performs the MH for the data vector. 
#y.curr: current data vector
#Beta: vector of reg coefficient parameters
#sig: sigma-the standard deviation parameter
#l1obs: conditioning statistic-estimate of Beta from observed data
#s1obs: conditioning statistic-estimate of sigma from observed data
#X: regression design matrix (i)
#log.prop.den.curr: value of the proposal density for the current data vector (free of other parameters so it can be saved from iteration to iteration to save computation time)
#proj: the projection matrix onto the deviation space (aka the least squares residual space)
#fn.psi: custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber 
#fn.chi: custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2

#psi, psi function returning psi(x)/x. Must match fn.spi
#scaleEst: 'Huber', must match fn.chi. 
    #inputs should be updated to get rid of redundancy. 
#maxit: max iterations for the estimation of the summary statistics


# returns: list
# y.curr: the new (if accepted) or same (if rejected) data vector
# a: indicator of acceptance of proposed value 1 if accepted, 0 if rejected
# log.prop.den.curr: proposal density on log scale of the returned data vector
# rat: MH ratio.
# log.rat1: log ratio of proposal densities. 
fn.one.rep.y<-function(y.curr,Beta,sig,l1obs, s1obs,X, log.prop.den.curr, proj,fn.psi,fn.chi,psi, scaleEst,maxit=400)
{
  n<-length(y.curr)
  p<-ncol(X)
  # propose new y vector
  y.prop <- rnorm(n) 
  y.prop <- fn.comp.ystst(y.prop,X,l1obs, s1obs, psi,scaleEst,maxit)  
  # compute acceptance ratio in steps
  log.rat1 <- log.fn.eval.lik(y.prop,X,Beta, sig)-log.fn.eval.lik(y.curr,X,Beta, sig)
  
  log.prop.den.prop<-log.prop.den(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi, n,p)
  
  log.rat2 <- log.prop.den.curr - log.prop.den.prop
  ratPre<-exp(log.rat1+log.rat2)
  rat <- min(c(1,ratPre))
  tmp <- runif(1) 
  if (tmp < rat) {y.curr <- y.prop; a<-1;log.prop.den.curr<-log.prop.den.prop} else {
    a<-0
  } 
  return(list(y.curr,a,log.prop.den.curr,rat, log.rat1))
}

#---------------------------------
#slightly faster. must include Qt as an input where columns of Q form orthonormal basis for X #also uses fn.radius2...which might be faster than fn.radius too
#---------------------------------

#updated, faster version of log.prop.den
#Additional input is  Qt (see fn.attenuation2 )


log.prop.den2 <-function(y.prop,X, proj,Qt, l1obs, s1obs,fn.psi,fn.chi,n,p)
  {
  log(fn.attenuation2(y.prop,X,proj,Qt, l1obs,s1obs,fn.psi, fn.chi))-(n-p-1)*log(fn.radius2(y.prop,proj))

}

#updated, faster version of fn.one.rep.y
#Additional input is  Qt (see fn.attenuation2)

fn.one.rep.y2<-function(y.curr,Beta,sig,l1obs, s1obs,X, log.prop.den.curr, proj, Qt,fn.psi,fn.chi,psi, scaleEst,maxit=400)
{
  n<-length(y.curr)
  p<-ncol(X)
  # propose new y vector
  y.prop <- rnorm(n)
  y.prop <- fn.comp.ystst(y.prop,X,l1obs, s1obs, psi,scaleEst,maxit)  
  # compute acceptance ratio in steps
  log.rat1 <- log.fn.eval.lik(y.prop,X,Beta, sig)-log.fn.eval.lik(y.curr,X,Beta, sig)
  
  #log.prop.den.curr <-log(fn.attenuation(y.curr,X, proj=proj, l1obs,s1obs))-(n-p-1)*log(fn.radius(y.curr, X))
  #log.prop.den.prop<-log.prop.den(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi, n,p)
  log.prop.den.prop<-log.prop.den2(y.prop,X, proj,Qt, l1obs, s1obs,fn.psi,fn.chi, n,p)
  
  log.rat2 <- log.prop.den.curr - log.prop.den.prop
  ratPre<-exp(log.rat1+log.rat2)
  #print(rat)
  rat <- min(c(1,ratPre))
  tmp <- runif(1) 
  if (tmp < rat) {y.curr <- y.prop; a<-1;log.prop.den.curr<-log.prop.den.prop} else {
    a<-0
  }  
  #  update if needed
  return(list(y.curr,a,log.prop.den.curr,rat, log.rat1))
}



