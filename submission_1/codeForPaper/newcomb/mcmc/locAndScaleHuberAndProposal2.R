####################################
#John Lewis August 1st 2013
#Beginings of code to run correct mcmc algorithm when conditing on robust statistics
#Focus here is on Huber's loss and Huber's proposal 2 to form the estimating equations for the location and scale in the location and scale setting
#will compare results to grid point estimation to verify the match
######################################

#The code will develop the sampling function to sample from the following full conditional in the gibbs step

#[Y|T(y), theta]
#where the model is

#y_i~iid N(mu, sigma^2) i=1,2,3,...,n
#priors
#mu~
#sigma^2~
#note that since conditional on mu and sigma^2, Y is independent of any other hyperparameters, sampling from [Y|T(y), theta] should be independent of any other hyperparameters set in the model

##################################
#As laid out, calculations are done in two steps

#1) occurs in the deviation space: a stretch and an attenuation
#2) is the attenuation from the deviation space to the original space

#calculations involve the gradients of the location and scale statistics wrt the each data value


#####################################
#Functions needed for the calculations of the gradients of the location and scale statistics
######################################

library(MASS)

#note: psi.huber is available in the mass library: I write my own.
fn.psi.huber<-function(u,k=1.345){
  pmax(-k, pmin(k,u))
}

# x<-seq(-4,4, by=.001)
# plot(x,fn.psi.huber(x), type='l', lwd=2)
# abline(h=c(1.345,-1.345))
# abline(v=c(1.345,-1.345))

fn.psi.huber.prime<-function(u, k=1.345){
  abs(u)<=k
}

# system.time(plot(x,fn.psi.huber.prime(x), type='l',f lwd=2))
# system.time(a<-fn.psi.huber.prime(x))
# abline(v=c(1.345,-1.345))


fn.chi.prop2.prime<-function(u, k=1.345){
  2*fn.psi.huber.prime(u,k)*fn.psi.huber(u,k)
}

#the linear equation to solve is A*deriv=b
#grad is the deriv of location wrt y_i and deriv of scale wrt y_i
#function should only imput y's that satisfy l1obs and s1obs
fn.grads<-function(y, l1obs, s1obs, k=1.345){
  prime1<-fn.psi.huber.prime((y-l1obs)/s1obs)
  A11<-s1obs*sum(prime1)
  A12<-sum(prime1*(y-l1obs))
  prime2<-fn.chi.prop2.prime((y-l1obs)/s1obs)
  A21<-s1obs*sum(prime2)
  A22<-sum(prime2*(y-l1obs))
  Ainv<-solve(matrix(c(A11,A12,A21,A22),2,2,byrow=TRUE))
  rhs<-s1obs*rbind(prime1,prime2)
  return(Ainv%*%rhs)
  }

# n<-20
# y.curr<-rnorm(n)
# tst<-rlm(y.curr~1,scale.est='Huber')
# ts<-fn.grads(y.curr, coef(tst), tst$s, k=1.345)
# sum(ts[2,])
# 
# set.seed(1)
# x<-rnorm(n)
# system.time(fit<-rlm(x~1, scale.est='Huber'))
# l1obs<-coef(fit)
# s1obs<-fit$s
# 
# system.time(fn.grads(x, l1obs, s1obs, k=1.345))
# tst<-fn.grads(x, l1obs, s1obs, k=1.345)
# sum(tst[2,])
# l1obs<- coef(tst);s1obs<-tst$s
# n<-length(y.curr)
# # fn.grads(y.curr,as.matrix(rep(1,n),n,1),l1obs, s1obs, k=1.345, psi=fn.psi.huber, chi=fn.chi.prop2)
# # fn.grads2(y.curr,l1obs, s1obs, k=1.345)
# 


#computes y** from y*
fn.compute.ystst<-function(yst,l1obs,s1obs){
  fit.yst<-rlm(yst~1, scale.est='Huber', maxit=400)
  l1yst<-coef(fit.yst)
  s1yst<-fit.yst$s
  res<-yst-l1yst
  ystst<-s1obs*res/s1yst+l1obs
  return(ystst)
}


#computes the radius in z space
fn.radius<- function(ystst) 
{
  rad <-sqrt(sum((ystst-mean(ystst))^2))
  return(rad)
}

#construct the transformation matrix
V<-contr.helmert(n, contrasts = TRUE, sparse =FALSE)
V<-cbind(rep(1,n),V)
column.norms<-apply(V, 2, FUN=function(x){sqrt(sum(x^2))})
V<-sweep(V, 2, STATS=column.norms, FUN='/')
Vt<-t(V)

fn.comp.v<-function(ystst,Vt=Vt){
  v<-Vt%*%ystst
  return(v)
}

fn.radius<- function(ystst) 
{
  rad <-sqrt(sum((ystst-mean(ystst))^2))
  return(rad)
}

fn.radius2<-function(vstst){
  rad <-sqrt(sum(vstst[-1]^2))
  return(rad)
}
# vstst<-fn.comp.v(ystst,Vt)
# fn.radius2(vstst)
# fn.radius(ystst)
# yst<-rnorm(n)
# fityst<-rlm(yst~1, scale.est='Huber')
# s1obs<-fityst$s
# l1obs<-coef(fityst)
# ystst<-fn.compute.ystst(yst, l1obs,s1obs)
# tst<-rlm(ystst~1, scale.est='Huber')
# coef(tst)
# l1obs
# tst$s
# s1obs
# Vt%*%grad.y.l
# Vt%*%grad.y.s

# x->ystst
#computes and returns the product of the two attenuations attenuation
fn.attenuation <- function(ystst, Vt=Vt, l1obs,s1obs)
{
  grads.ystst<-fn.grads(ystst, l1obs, s1obs, k=1.345)
  grad.y.l<-grads.ystst[1,]
  grad.y.s<-grads.ystst[2,]
  
  unit.a<-fn.comp.v(ystst,Vt)[-1] 
  unit.a<-unit.a/sqrt(sum(unit.a*unit.a)) #normalize
  
  grad.v.s<-(Vt%*%grad.y.s)[-1]
  unit.b<-grad.v.s/sqrt(sum(grad.v.s*grad.v.s))
  atten1 <-abs(sum(unit.a * unit.b)) 
  #atten1 <-sqrt(1-sum(unit.a * unit.b)^2) 
  qr.A<-qr(cbind(grad.y.s,grad.y.l))
  Qmat<-qr.Q(qr.A)
  unit1<-Qmat[,2]
  unit2<-Vt[1,] 
  atten2<-abs(sum(unit1 * unit2))
  #atten2<-sqrt(1-sum(unit1 * unit2)^2)
  return(atten1*atten2)
  }



# testing the wedge product version
# qr.A<-qr(cbind(grad.y.s,grad.y.l, diag(n)[,-c(1:2)]))
# p<-1
# Qmat2<-qr.Q(qr.A)[,-c(1:(1+1))]
# B<-W%*%t(W)%*%Qmat2
# sqrt(det(t(B)%*%B))
# grad.v.s<-(Vt%*%grad.y.s)[-c(1:p)]
# unit.b<-grad.v.s/sqrt(sum(grad.v.s*grad.v.s))
# Estars<-crossprod(t(Vt),Qmat)[-c(1:p),]
# abs(det(cbind(unit.b,Estars)))
# 

#wrote this on 11/4/2013: this caluclate the second attenuation by attempting to caluculate the derivatives
#of y**=z**+l-l(z**) directly:
#I could not get this to work
#################################################

fn.attenuation2 <- function(ystst, Vt=Vt, l1obs,s1obs)
{
  grads.ystst<-fn.grads(ystst, l1obs, s1obs, k=1.345)
  grad.y.l<-grads.ystst[1,]
  grad.y.s<-grads.ystst[2,]
  
  unit.a<-fn.comp.v(ystst,Vt)[-1] 
  unit.a<-unit.a/sqrt(sum(unit.a*unit.a)) #normalize
  grad.v.s<-(Vt%*%grad.y.s)[-1]
  unit.b<-grad.v.s/sqrt(sum(grad.v.s*grad.v.s))
  atten1 <-abs(sum(unit.a * unit.b)) 
  
  grad.v.l<-(Vt%*%grad.y.l)
  
  W<-t(Vt)[,2:length(ystst)]
  zstst<-t(W%*%t(W)%*%ystst)
  l1obsZstst<-coef(rlm(t(zstst)~1, scale.est='Huber', maxit=400))
  fn.grads.z<-fn.grads(zstst,  l1obsZstst, s1obs, k=1.345)
  
  Jac<-W%*%t(W)-matrix(rep(fn.grads.z[1,], length(ystst)),length(ystst),length(ystst), byrow=TRUE)
  atten2<-abs(det(solve(Jac)))
  return(atten1*atten2)
}
# 
#################################################


# grad.y.s.norm<-grad.y.s/sqrt(sum(grad.y.s^2))
# grad.y.s.norm%*%Vt[1,]
# grad.y.s.norm%*%unit1
# Vt[1,]%*%unit1
# 
# # ytest<-ystst-mean(ystst)
# ytest<-ytest/sqrt(sum(ytest^2))
# plot(round(unit.a,5),round(ytest,5))
# abline(0,1)
#computed on the log scale
log.fn.eval.lik <- function(ystst, mu, sig)
{
  lik <- sum(dnorm(ystst,mean=mu,sd=sig, log=TRUE))
  return(lik)
}

#fn.eval.lik(y.prop,mu,sig)-fn.eval.lik(y.curr,mu,sig)





#####################################

fn.one.rep.y<- function(y.curr,mu,sig,l1obs, s1obs,n=length(y.curr))
{
  # propose new y vector
  y.prop <- rnorm(n)
  y.prop <- fn.compute.ystst(y.prop,l1obs, s1obs)  
  # compute acceptance ratio
  log.rat1 <- log.fn.eval.lik(y.prop,mu,sig)-log.fn.eval.lik(y.curr,mu,sig)
  
  log.prop.den.curr <-log(fn.attenuation(y.curr, Vt=Vt, l1obs,s1obs))-(n-2)*log(fn.radius(y.curr))
  log.prop.den.prop <-log(fn.attenuation(y.prop,Vt=Vt, l1obs,s1obs))-(n-2)*log(fn.radius(y.prop)) 
  log.rat2 <- log.prop.den.curr - log.prop.den.prop
  ratPre<-exp(log.rat1+log.rat2)
  #print(rat)
  rat <- min(c(1,ratPre))
  tmp <- runif(1) 
  if (tmp < rat) {y.curr <- y.prop; a<-1} else {
    a<-0
  } 
  #  update if needed
  return(list(y.curr,a))
}



# ######################################
# # mu=0
# # sig=1
#fn.one.rep.y<- function(y.curr,mu,sig,l1obs, s1obs,n=length(y.curr))
#{
 # # propose new y vector
 # y.prop <- rnorm(n)
 # y.prop <- fn.compute.ystst(y.prop,l1obs, s1obs)  
 # # compute acceptance ratio
 # log.rat1 <- fn.eval.lik(y.prop,mu,sig)-fn.eval.lik(y.curr,mu,sig)
 # 
 # log.prop.den.curr <--(n-2)*log(fn.radius(y.curr))
 # log.prop.den.prop <--(n-2)*log(fn.radius(y.prop)) 
 # log.rat2 <- log.prop.den.curr-log.prop.den.prop
 # ratPre<-exp(log.rat1+log.rat2)
 # #print(rat)
 # rat <- min(c(1,ratPre))
 # tmp <- runif(1) 
 # if (tmp < rat) {y.curr <- y.prop; a<-1} else {
 #   a<-0
 # } 
 # #  update if needed
 # return(list(y.curr,a))
#}



# #computed on the original scale
# fn.eval.lik <- function(ystst, mu, sig)
# {
#   lik <- prod(dnorm(ystst,mean=mu,sd=sig))
#   return(lik)
# }
# 
# 
# fn.one.rep.y<- function(y.curr,mu,sig,l1obs, s1obs,n=length(y.curr))
# {
#   # propose new y vector
#   y.prop <- rnorm(n)
#   y.prop <- fn.compute.ystst(y.prop,l1obs, s1obs)  
#   # compute acceptance ratio
#   rat1 <- fn.eval.lik(y.prop,mu,sig)/fn.eval.lik(y.curr,mu,sig)
#   
#   prop.den.curr <-fn.attenuation(y.curr, Vt=Vt, l1obs,s1obs)/(fn.radius(y.curr)^(n-2))
#   prop.den.prop <-fn.attenuation(y.prop,Vt=Vt, l1obs,s1obs)/(fn.radius(y.prop)^(n-2)) 
#   rat2 <- prop.den.curr/prop.den.prop
#   rat<-rat1*rat2
#   #print(rat)
#   rat <- min(c(1,rat))
#   tmp <- runif(1) 
#   if (tmp < rat) {y.curr <- y.prop; a<-1} else {
#     a<-0
#   } 
#   #  update if needed
#   return(list(y.curr,a))
# }
# 
# 
# 
# 


# 
# y.prop <- rnorm(length(y.curr))
# y.prop <- fn.compute.ystst(y.prop,l1obs, s1obs)
# prod(dnorm(y.prop,mean=2,sd=4))
# fn.one.rep.y(x,mu=1,sig=2,l1obs, s1obs,n=n)


# date()
# n<-66
# total<-1e3
# center=0
# stDev=1
# l1obs<-0
# s1obs<-1 
# res.mat <- matrix(rep(0,n*total),ncol=n)
# y.curr <- rnorm(n)
# y.curr <- fn.compute.ystst(y.curr,l1obs, s1obs)
# y.curr<-list(y.curr,0)
# system.time(for (i in 1:total)
# {y.curr <- fn.one.rep.y(y.curr[[1]],mu=center,sig=stDev,l1obs, s1obs,n=n)
#  res.mat[i,] <- y.curr[[1]]}
# )
# 
# rlm(res.mat[550,]~1,scale.est='Huber')
# sum(round(diff(res.mat[,1]),10)!=0)/total
# 
# plot(round(diff(res.mat[,1]),10))

