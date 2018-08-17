# Test simulation Source Code using simulated data from the model with known values

rm(list=ls())
setwd("~/Box Sync/research/snm/BlendedParadigmWriting/codeForPaper/simulation")
source('simulationSourceCode.R')

#simulate data directly from the model
set.seed(1)
nGroups<-90 # make an even number for ease
muTrue<-0
tau2True<-1
#sigma2True<-c(rep(c(1,4), (nGroups/2)))
mu<-1
v<-.00001
a<-mu^2/v-2
b<-mu*(a-1)

#first half the groups have ni=20, other half have ni=50. Half of the larger variance groups have the smaller sample size, half have larger sample size
# n<-c(rep(20, floor(nGroups/3)), rep(50, nGroups-floor(nGroups/3)),rep(100, nGroups-floor(nGroups/3)))
# 
# YList<-list()
# for(i in 1:nGroups){
#   YList[[i]]<-rnorm(n[i],thetaTrue[i], sqrt(sigma2True[i]))
# }
sigma2True<-4
YList<-fn.gen.data(sigma2True)
factorsList<-YList$factorsList
factorsMat<-matrix(unlist(factorsList), nrow=90, ncol=3, byrow = TRUE)
p<-factorsMat[,1]
n<-factorsMat[,2]
m<-factorsMat[,3]
thetaTrue<-YList$theta
hist(thetaTrue)
YList<-YList$yList
# Fit the complete data model----

tst<-fn.one.rep.complete.MCMC(YList, theta=thetaTrue, sigma2=sigma2True,mu=muTrue, tau2=tau2True,a=a,b=b)
#tst

#intitial values
#thetaInt<-sapply(YList, mean)
thetaInt<-rep(2, nGroups)
muInt<-5

st<-Sys.time()
completeFit<-fn.complete.MCMC(YList, nkeep=1000, nburn=1000, theta=thetaInt, sigma2=NULL, mu=muInt, tau2=NULL,a=a,b=b)
en<-Sys.time()-st
en
#
# Trace plots with true values indicated ----
#
par(mfrow=c(2,2))
for(i in 1:nGroups){
  plot(completeFit$theta[,i],type='l', main=paste('p=', p[i], 'n=',n[i], 'm=', m[i]), ylab=paste('theta', i), xlab='iteration')
  abline(h=thetaTrue[i], col='blue')
}


par(mfrow=c(2,2))
for(i in 1:nGroups){
  plot(completeFit$sigma2[,i],type='l',  main=paste('p=', p[i], 'n=',n[i], 'm=', m[i]), ylab=paste('sigma2', i), xlab='iteration')
abline(h=sigma2True, col='blue')
}

par(mfrow=c(1,1))
plot(completeFit$mu,type='l', main='', ylab='mu', xlab='iteration')
abline(h=muTrue, col='blue')

postMeansTheta<-apply(completeFit$theta, 2, mean)
mean(postMeansTheta)
var(postMeansTheta)
abline(h=mean(postMeansTheta))

plot(completeFit$tau2,type='l', main='', ylab='tau2', xlab='iteration')
abline(h=tau2True, col='blue')
abline(h=var(postMeansTheta), col='red')
abline(h=var(sapply(YList,mean)), col='green')



par(mfrow=c(1,1))
plot(postMeansTheta, thetaTrue, pch=16)
abline(0,1, col=4, lwd=2)

mean((postMeansTheta-thetaTrue)^2)
plot((postMeansTheta-thetaTrue)^2)

postMeansSigma2<-apply(completeFit$sigma2, 2, mean)
par(mfrow=c(1,1))
#plot(postMeansSigma2, sigma2True, ylim=c(0, 6), xlim=c(0,6), pch=16)
#abline(0,1)

#
# Incomplete model ------
#

#thetaInt<-sapply(YList, mean)
set.seed(1)
thetaInt<-rep(2, nGroups)
muInt<-5
source('simulationSourceCode.R')
st<-Sys.time()
incompleteFit1<-fn.Incomplete.MCMC(YList,
                             regEst='Tukey',
                             scaleEst='Huber',
                             nkeep=1000, 
                             nburn=1000, 
                             theta=thetaInt, 
                             sigma2=NULL, 
                             mu=muInt, 
                             tau2=NULL, 
                             printEach=100,
                             maxit=400,
                             a=a,
                             b=b)
en1<-Sys.time()-st
en1
#incompleteFit1$robustFits
acceptRates<-apply(incompleteFit1$yAccept, 2, mean)
plot(acceptRates, pch=16)
range(acceptRates)
factorsMat
sort(acceptRates, decreasing = FALSE)
#
# Trace plots with true values indicated ----
#
bhats<-sapply(incompleteFit1$robustFits,FUN=function(x) x$coefficients)
par(mfrow=c(2,2))
for(i in 1:nGroups){
  plot(incompleteFit1$theta[,i],type='l',main=paste('p=', p[i], 'n=',n[i], 'm=', m[i]), ylab=paste('theta', i), xlab='iteration')
  #lines(completeFit$theta[,i],type='l',col=2)
  abline(h=thetaTrue[i], col='blue')
  abline(h=bhats[i], col='red')
}


par(mfrow=c(2,2))
for(i in 1:nGroups){
  plot(incompleteFit1$sigma2[,i],type='l',  main=paste('n=',n[i]), ylab=paste('sigma2', i), xlab='iteration')
  abline(h=sigma2True, col='blue')
}

par(mfrow=c(1,1))
plot(incompleteFit1$mu,type='l', main='', ylab='mu', xlab='iteration')
lines(completeFit$mu,type='l',col='red')
abline(h=muTrue, col='blue')

par(mfrow=c(1,1))
plot(incompleteFit1$tau2,type='l', main='', ylab='mu', xlab='iteration')
lines(completeFit$tau2, col='red')
abline(h=tau2True, col='blue')



postMeansTheta<-apply(incompleteFit1$theta, 2, mean)
postMeansThetaComplete<-apply(completeFit$theta, 2, mean)
plot(postMeansTheta,postMeansThetaComplete, pch=16)
abline(0,1, col=4, lwd=2)

par(mfrow=c(1,1))
plot(postMeansTheta, thetaTrue, pch=16)
abline(0,1, col=4, lwd=2)
par(mfrow=c(1,1))
plot(postMeansThetaComplete, thetaTrue, pch=16)
abline(0,1, col=4, lwd=2)

robustFittedValues<-sapply(incompleteFit1$robustFits, FUN=function(x) x$coefficients)


# Evaluation metrics

# 1
mean((postMeansTheta-thetaTrue)^2)
mean((robustFittedValues-thetaTrue)^2)
mean((postMeansThetaComplete-thetaTrue)^2)


median((postMeansTheta-thetaTrue)^2)
median((robustFittedValues-thetaTrue)^2)
median((postMeansThetaComplete-thetaTrue)^2)

sum((postMeansTheta-thetaTrue)^2>(robustFittedValues-thetaTrue)^2)

plot((postMeansTheta-thetaTrue)^2,(robustFittedValues-thetaTrue)^2, ylim=c(0,1), xlim=c(0,.), pch=19)
abline(0,1, col=4)

boxplot((postMeansThetaComplete-thetaTrue)^2,(robustFittedValues-thetaTrue)^2,(postMeansTheta-thetaTrue)^2, ylim=c(0,.2))


normDiffs<-(postMeansThetaComplete-thetaTrue)^2
normDiffsFit<-lm(normDiffs~n+p+m)
summary(normDiffsFit)

robDiffs<-(postMeansTheta-thetaTrue)^2
robDiffsFit<-lm(robDiffs~n+p+m)
summary(robDiffsFit)

rlmDiffs<-(robustFittedValues-thetaTrue)^2
rlmDiffsFit<-lm(rlmDiffs~n+p+m)
summary(rlmDiffsFit)



mean((postMeansTheta-thetaTrue)[n==20]^2)
mean((robustFittedValues-thetaTrue)[n==20]^2)
mean((postMeansThetaComplete-thetaTrue)[n==20]^2)

mean((postMeansTheta-thetaTrue)[n==50]^2)
mean((robustFittedValues-thetaTrue)[n==50]^2)
mean((postMeansThetaComplete-thetaTrue)[n==50]^2)


mean((postMeansTheta-thetaTrue)[n==100]^2)
mean((robustFittedValues-thetaTrue)[n==100]^2)
mean((postMeansThetaComplete-thetaTrue)[n==100]^2)


#2 K-L divergance





var(postMeansTheta)
var(postMeansThetaComplete)



plot(incompleteFit1$tau2,type='l', main='', ylab='tau2', xlab='iteration')
lines(completeFit$tau2,type='l',col='red')

abline(h=tau2True, col='blue')
abline(h=var(postMeansTheta), col='red')
abline(h=var(sapply(YList,mean)), col='green')


mean((postMeansTheta-thetaTrue)^2)
plot((postMeansTheta-thetaTrue)^2)

postMeansSigma2<-apply(incompleteFit1$sigma2, 2, mean)
par(mfrow=c(1,1))
plot(postMeansSigma2, pch=16)
abline(h=sigma2True)





#Look at Trace Plots from the actual simulation

load('simulationWorkspace.RData')

sig2True=10
thetaTrue<-get(paste0('thetaTrueSig2_', sig2True))
huberFit<-get(paste0('restrictedFitHuberSig2_', sig2True))
tukeyFit<-get(paste0('restrictedFitTukeySig2_', sig2True))
nFit<-get(paste0('normalTheoryFitSig2_', sig2True))
par(mfrow=c(2,2))
for(i in 1:90){
  plot(huberFit$theta[,i], type='l')
  abline(h=thetaTrue[i], col='4')
}
par(mfrow=c(2,2))
for(i in 1:90){
  plot(tukeyFit$theta[,i], type='l')
  abline(h=thetaTrue[i], col='4')
}
  
par(mfrow=c(2,2))
for(i in 1:90){
  plot(nFit$theta[,i], type='l')
  abline(h=thetaTrue[i], col='4')
}




