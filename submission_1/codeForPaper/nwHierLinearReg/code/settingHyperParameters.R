#----------------------
#Setting priors for hierarchical model on NW data
#----------------------

library(MASS)
library(MCMCpack)
#-fit robust regressions to each state seperately
rm(list=ls())
path.to.NonHierworkspaces<-'../../nwLinearReg/workSpaces'
setwd(path.to.NonHierworkspaces)
load("nwdataPaper1PriorConstructionWorkSpace.RData")
load("nwdataPrepPaper1DataCleanWorkspace.RData") #not provided 
path.to.workspaces<-'../../nwHierLinearReg/workSpaces'
setwd(path.to.workspaces)
#View(priorDataSet)
X<-model.matrix(priorRlm)
sd(analysisSet$sqrt_Count2012)
dim(X)
# rm(c)
p<-4
states<-unique(priorDataSet$Primary_Agency_State)
nStates<-length(states)
stateFits<-list(); length(stateFits)<-nStates
n_is<-rep(NA, length(states))
#make sure there are enough cases for good estimate of sigma2??
for(i in 1:nStates){
  dataset<-subset(priorDataSet, priorDataSet$Primary_Agency_State==states[i])
  if(nrow(dataset)>40){
  stateFits[[i]]<-rlm(sqrt_Count2010~sqrt_Count2008+Associate_Count+ageMeasure,psi=psi.bisquare,scale.est="Huber",data=dataset, maxit=400)
  n_is[i]<-nrow(dataset)
  }
}


whichFit<-unlist(lapply(stateFits, FUN=function(x) class(x)[1]=='rlm'))
sum(whichFit)


which(is.na(n_is))
betaHats<-matrix(NA, nStates, 4)
sigma2Hats<-rep(NA,nStates)
for(i in c(1:nStates)[whichFit]){
  {
    betaHats[i,]<-coef(stateFits[[i]])
    sigma2Hats[i]<-stateFits[[i]]$s^2
  }
}
sum(is.na(betaHats[,1]))

cbind(sigma2Hats,n_is)
sum(is.na(n_is))
sum(is.na(betaHats[,1]))
sum(is.na(sigma2Hats))
all.equal(!is.na(betaHats[,1]),!is.na(n_is))

# #get rid of the state with huge sd?
whichFit[19]<-FALSE

betaHats<-betaHats[whichFit,]
n_is<-n_is[whichFit]
length(n_is)
dim(betaHats)

tst<-sigma2Hats[!is.na(sigma2Hats)]
curve(dinvgamma(x, a0Star,b0Star), from=0, to=.05)
points(tst/15^2, rep(20, length(tst)))
plot(sapply(tst/15^2, FUN=function(y) integrate(dinvgamma,lower=0, upper=y,a0Star,b0Star)$value))

wts<-n_is/sum(n_is)


sigma2Hats<-sigma2Hats[whichFit] 
str(sigma2Hats)
sort(sigma2Hats)
length(sigma2Hats)
nrow(betaHats)

# ------------------------
# prior for beta_i's
# ------------------------

nrow(priorDataSet)
vcov(priorRlm)
p<-nrow(vcov(priorRlm))
#SigmaBetaHat1<-priorRlm$s^2*solve(t(X)%*%X)
SigmaBetaHat1<-vcov(priorRlm)


delta_is<-apply(betaHats, 1, FUN=function(x) x-coef(priorRlm))

dList<-list()
for(i in 1:ncol(delta_is)){
  dList[[i]]<-delta_is[,i]%*%t(delta_is[,i])
}

SigmaDelta<-matrix(rep(0, p*p), p,p)
for(i in 1:length(dList)){
  SigmaDelta<-SigmaDelta+dList[[i]]*wts[i]
}
  
SigmaDelta  
K<-length(n_is)
#SigmaDelta<-K/(K-p)*SigmaDelta
#another way?
SigmaDelta2<-matrix(rep(0, p*p), p,p)
stateFits2<-stateFits[whichFit]
for(i in 1:length(stateFits2)){
  SigmaDelta2<-SigmaDelta2+vcov(stateFits2[[i]])*wts[i]
}
SigmaDelta2


z<-(det(SigmaDelta)/det(SigmaBetaHat1))^(1/p) #/nrow(X) #divide by nrow(X)??
swSq<-sum(wts^2)

nrow(X)
dim(X)
SigmaDelta
SigmaBetaHat1*z
round(SigmaDelta-SigmaBetaHat1*z, 4)

z/(1+z)

mnb<-z/(nrow(X))
swSq<-1
mnb*swSq #want this to be the mean for b*c
mnb*swSq<1/nrow(X)
fn.compute.ab<-function(mu, psi){
  a<-mu*psi
  b<-psi-a
  c(a,b)
}


ab<-fn.compute.ab(.6, 10)

a<-ab[1]
b<-ab[2]
a
b

all.equal(a/(a+b),.6) #mean
all.equal(a*b/((a+b)^2*(a+b+1)),.6*(1-.6)/((10+1))) #var
swSq<-1 #sum(wts^2)
mnb*swSq

#--------------------------------------------
#fix psi_bstr
#--------------------------------------------
psi_bstr<-20  
#--------------------------------------------
#
#fix mu_bstr
#---------------------------------
mu_bstr<-mnb*swSq
fn.compute.ab(mu_bstr,psi_bstr)

# -------------------------
# prior for sigma2_i's
# -------------------------

a0<-a0Star
sy<-sd(analysisSet[,"sqrt_Count2012"])
b0<-sy^2*b0Star
a0
b0
hist(sigma2Hats)
#given z~N(0,1) compute inverse gamma(ao, b0) random variable
invGam<-function(z, a0, b0){
  c<-pnorm(z)
  x<-qgamma(c, a0, scale=1/b0)
  1/x 
}


K<-length(sigma2Hats)
mean(sigma2Hats)
sd(sigma2Hats)
median(sigma2Hats)
mad(sigma2Hats)
summary(sigma2Hats)
#mean and sd of my inverse gamma
b0/(a0-1)
sqrt(b0^2/((a0-1)^2*(a0-2)))
priorRlm$s^2

curve(dinvgamma(x, a0,b0), from=0,to=30)
points( sigma2Hats,rep(.1,length(sigma2Hats)))
all.equal(sapply(tst/sy^2, FUN=function(y) integrate(dinvgamma,lower=0, upper=y,a0Star,b0Star)$value)
,sapply(tst, FUN=function(y) integrate(dinvgamma,lower=0, upper=y,a0,b0)$value))
sort(sapply(tst, FUN=function(y) integrate(dinvgamma,lower=0, upper=y,a0,b0)$value))


#prior parameters for beta on mu_rho
#idea: compute the z_i's; calculate the MLE of rho, make this the mean of mu_rho
fn.compute.Z<-function(sigma2,a0,b0){
  X<-integrate(function(x) {dinvgamma(x, shape=a0, scale = b0)},0, sigma2)$value #rate=b0
  qnorm(X)
  
}

Z<-sapply(sort(sigma2Hats)[1:17], FUN=fn.compute.Z, a0, b0)
K<-length(Z)
J<-matrix(1, K,K)
log.like.rho<-function(rho){
K*log(1-rho)+log(1+K*rho/(1-rho))+(1-rho)*t(Z)%*%Z+rho*t(Z)%*%J%*%Z  
}



curve(-log.like.rho(x),0.001,.999)
#maximize
op<-optim(.5, log.like.rho,control=list(fnscale=-1), lower=0, upper=.9999, method = c("L-BFGS-B"), hessian=TRUE)
muRho<-as.numeric(op$par)
varRho<--2*as.numeric(op$hessian^-1)
psiRho<-muRho*(1-muRho)/varRho-1
w1<-muRho*psiRho
w2<-psiRho-w1
w1;w2
curve(dbeta(x, w1,w2))
w1/(w1+w2)
w1*w2/((w1+w2)^2*(w1+w2+1))

#prior parameters for gamma on psi_rho
#how to choose these??
a_psir<-5
b_psir<-a_psir/psiRho
curve(dgamma(x, a_psir, b_psir), from=0, to=30)
a_psir/b_psir
a_psir/b_psir^2

# fn.compute.ab(.5,24)
# curve(dbeta(x, 12,12))
n<-nrow(X)
sigmahat<-priorRlm$s
meanPsiSq<-mean(((residuals(priorRlm)/sigmahat)*psi.bisquare(residuals(priorRlm)/sigmahat))^2)
meanPsiPrime<-mean(psi.bisquare(residuals(priorRlm)/sigmahat, deriv=1))
vhat<-(n/(n-p))*sigmahat^2*meanPsiSq/(meanPsiPrime^2)
Sigma0<-vhat*solve(t(X)%*%X)
vcov(priorRlm) #about the same
Sigma0<-n*Sigma0
Sigma0
mu0<-coef(priorRlm)

#---------------------------
#don't run if you want to generate data from the prior 
#---------------------------

rm(list=setdiff(ls(),c("mu_bstr","psi_bstr","w1","w2","a_psir","b_psir")))
save.image("wsnwdataHierPriorConstruction.RData")


# #remaining code generates data from the resulting model to see if the priors are reasonable. 
# getwd()
# 
# 
# stop()
# # --------------------------
# set.seed(1)
# J<-matrix(1, K,K)
# Id<-diag(K)
# samps<-2000
# 
# meanSimBetaHats<-matrix(NA, samps, p)
# sdSimBetaHats<-matrix(NA, samps, p)
# medSimBetaHats<-matrix(NA, samps, p)
# madSimBetaHats<-matrix(NA, samps, p)
# 
# meanSimSigma2Hats<-numeric(samps)
# sdSimSigma2Hats<-numeric(samps)
# 
# medSimSigma2Hats<-numeric(samps)
# madSimSigma2Hats<-numeric(samps) 
# 
# simBetaHat<-matrix(NA, K, p)
# simSigma2Hat<-numeric(K)
# rvec<-numeric(samps)
# for(i in 1:samps){
#   #generate sigma2_i's
#   mu_rho<-rbeta(1, w1, w2) 
#   psi_rho<-rgamma(1, a_psir, b_psir)  
#   ab<-fn.compute.ab(mu_rho, psi_rho)
#   a<-ab[1]
#   b<-ab[2]
#   rho<-rbeta(1, a, b)
#   rvec[i]<-rho
#   Sigma_rho<-rho*J+(1-rho)*Id
#   zSt<-mvrnorm(1, rep(0, K), Sigma_rho)
#   sigma2_is<-sapply(zSt, FUN=invGam, a0=a0, b0=b0)
#   
#   #generate beta_i's
#     #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
#     #psi_bstr<-rgamma(1, a_psib,b_psib)
#     v12<-fn.compute.ab( mu_bstr,  psi_bstr) #same every time
#     bc<-rbeta(1, v12[1],v12[2])
#     b<-bc/swSq
#     a<-1-bc
#     beta2<-mvrnorm(1, mu0, a*Sigma0)
#     beta_is<-mvrnorm(K, beta2, b*Sigma0)
#   
# #generate data 
#   for(j in 1:K){
#     Xd<-model.matrix(stateFits2[[j]])
#     y<-Xd%*%beta_is[j,]+rnorm(nrow(Xd),sd=sqrt(sigma2_is[j]))
#     fit<-rlm(Xd,y,psi=psi.bisquare, scale.est='Huber', maxit=1000)
#     simBetaHat[j,]<-coef(fit)
#     simSigma2Hat[j]<-fit$s^2
#   }
# #summarize the fits  
# meanSimBetaHats[i,]<-colMeans(simBetaHat)  
# sdSimBetaHats[i,]<-apply(simBetaHat,2,sd)  
# 
# medSimBetaHats[i,]<-apply(simBetaHat,2, median )  
# madSimBetaHats[i,]<-apply(simBetaHat,2,mad)  
#   
# meanSimSigma2Hats[i]<-mean(simSigma2Hat)  
# sdSimSigma2Hats[i]<-sd(simSigma2Hat)
#   
# medSimSigma2Hats[i]<-median(simSigma2Hat)  
# madSimSigma2Hats[i]<-mad(simSigma2Hat)  
# }  
# 
# #-----betas
# 
#   j<-1
# 
# #pdf('simBetas.pdf')
# for(j in 1:p){
#   plot(meanSimBetaHats[,j], sdSimBetaHats[,j], pch=19, main=paste('betaMeans',j), cex=.5, ylab='between group sd', xlab='mean')
#   points(mean(betaHats[,j]), sd(betaHats[,j]), col='red', pch=19, cex=1.5)
#  # points(median(betaHats[,j]), mad(betaHats[,j]), col='blue', pch=19, cex=1)
# }
#  #dev.off()
# j=1
# range(sdSimBetaHats[,j])
# sort(betaHats[,j])
# cbind(betaHats[,j], n_is)
# for(j in 1:p){
#   plot(medSimBetaHats[,j], madSimBetaHats[,j], pch=19, main=paste('betaMeds',j), cex=.5,ylab='between group mad', xlab='median')
#   points(median(betaHats[,j]), mad(betaHats[,j]), col='blue', pch=19,cex=1)
# }
# #---------
# 
# #----Sigma2s
# 
# #pdf("simSigma2s.pdf")
# plot(meanSimSigma2Hats,sdSimSigma2Hats, pch=19, main='meanSigma2s', cex=.5)
# points(mean(sigma2Hats), sd(sigma2Hats),  col='red', pch=19, cex=1.5)
# # points(median(sigma2Hats), mad(sigma2Hats),  col='blue', pch=19, cex=2)
# 
# 
# plot(medSimSigma2Hats,madSimSigma2Hats, pch=19, main='medianSigma2s', cex=.5)
# points(median(sigma2Hats), mad(sigma2Hats),  col='blue', pch=19,cex=1.5)
# #dev.off()
# 
