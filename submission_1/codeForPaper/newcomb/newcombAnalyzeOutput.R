#######
#Running from start to finish generates plots for the newcomb application section in the paper
#paramSet controls which set of hyper-parameters to use
#plots saved in figures folder (user degined), table for predictive distributions saved in the paper_version2 folder

library(tables) #for easy importing to latex
rm(list=ls())
paramSet<-2 #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper; this number is appended on the end of each workspace
trim=2 #For LTS: what is k? trim=1: k=34, trim=2: k=62

#################################################
library(MASS)
data(newcomb)
#May need to change working directories

setwd("~/Box Sync/research/snm/BlendedParadigmWriting/codeForPaper/newcomb")


dirSampDirectory<-paste0(getwd(), "/directSampling/")


mcmcDirectory<-paste0(getwd(), "/mcmc/")

#tDistDirectory<-"~/Box Sync/research/snm/locationAndScale/dataAnalysis/newcomb/fittingTdist/"

tDistDirectory<-paste0(getwd(), "/fittingTdist/")


path.to.figures<-"~/Box Sync/research/snm/BlendedParadigmWriting/paper_version2/figures/"

nu<-5

######################################################
#read in the posterior marginal pdfs from Huber statistics from the direct sampling
#########################################################

load(paste(dirSampDirectory,"wsNewcombDirectSamplingHuber", paramSet, ".RData", sep=''))
eta
tau
alpha
beta

joint.post.Huber<-log.post.density.final
pdfThetaRestricted<-cbind(theta.grid,post.theta.pdf)
pdfSigma2Restricted<-cbind(sigma2.grid,post.sigma2.pdf)
theta.gridHuber<-theta.grid
sigma2.gridHuber<-sigma2.grid
delta.thetaHuber<-delta.theta
delta.sigma2Huber<-delta.sigma2
#estimated mean theta
delta.thetaHuber*sum(pdfThetaRestricted[,1]*pdfThetaRestricted[,2])

#compute likelihood function for sigma=s.obs
likelihood.estimate<-function(theta){
  sigma2=s.obs^2
  kDensityEst<-log(kernel_density((t.obs-theta)/sqrt(sigma2),log.s.obs-.5*log(sigma2))/sqrt(sigma2))
  return(exp(kDensityEst))
}
theta.gridHuber2<-seq(min(theta.gridHuber), max(theta.gridHuber), length.out=75)
# likelihood.huber<-sapply(theta.gridHuber2, FUN=likelihood.estimate)
# plot(theta.gridHuber2, likelihood.huber, type='l')
# lines(pdfThetaRestricted[,1], pdfThetaRestricted[,2], lty=2)
# lines(theta.gridHuber, dnorm(theta.gridHuber, eta, tau), col='purple')

######################################################
#read in the posterior draws from Full Likelihood and restricted likelihood with Hubers/Tukey estimators: done by MCMC
#########################################################

load(paste(mcmcDirectory,"WorkSpaceNewcombFullNormalModel", paramSet, ".RData", sep=''))
#load("WorkSpaceNewcombHuberAndProposal2Analysis.RData")
#load("WorkSpaceNewcombTukeyAndProposal2Analysis.RData")
load(paste(mcmcDirectory,"WorkSpaceNewcombHuberAndProposal2AnalysisWithDelAnalysis", paramSet, ".RData", sep=''))
load(paste(mcmcDirectory,"WorkSpaceNewcombTukeyAndProposal2AnalysisWithDelAnalysis", paramSet, ".RData", sep=''))


######################################################
#read in the posterior marginal pdfs using LMS statistics from the direct sampling
#########################################################
# 
load(paste(dirSampDirectory,"wsNewcombDirectSamplingLMS", paramSet, ".RData", sep=''))
eta
tau
alpha
beta

joint.post.LMS<-post.density.final
pdfThetaLMS<-cbind(theta.grid,post.theta.pdf)
pdfSigma2LMS<-cbind(sigma2.grid,post.sigma2.pdf)
delta.thetaLMS<-delta.theta
delta.sigma2LMS<-delta.sigma2
theta.gridLMS<-theta.grid
sigma2.gridLMS<-sigma2.grid




#estimated mean theta
delta.thetaLMS*sum(pdfThetaLMS[,1]*pdfThetaLMS[,2])
kernel_density
 theta.gridLMS2<-seq(min(theta.gridLMS),max(theta.gridLMS), length.out=75)
# likelihood.LMS<-sapply(theta.gridLMS2, FUN=likelihood.estimate)
# lines(theta.gridLMS2, likelihood.LMS, type='l', col=3)
# lines(pdfThetaLMS[,1], pdfThetaLMS[,2], col=3, lty=2)


######################################################
#read in the posterior marginal pdfs using LTS statistics from the direct sampling
#########################################################
# 
if(trim==1) load(paste(dirSampDirectory,"wsNewcombDirectSamplingLTS", paramSet, ".RData", sep=''))
if(trim==2) load(paste(dirSampDirectory,"wsNewcombDirectSamplingLTS", paramSet,'trim_', trim,  ".RData", sep=''))

eta
tau
alpha
beta

joint.post.LTS<-post.density.final
pdfThetaLTS<-cbind(theta.grid,post.theta.pdf)
pdfSigma2LTS<-cbind(sigma2.grid,post.sigma2.pdf)
delta.thetaLTS<-delta.theta
delta.sigma2LTS<-delta.sigma2
theta.gridLTS<-theta.grid
sigma2.gridLTS<-sigma2.grid



#estimated mean theta
delta.thetaLTS*sum(pdfThetaLTS[,1]*pdfThetaLTS[,2])
kernel_density
theta.gridLTS2<-seq(min(theta.gridLTS), max(theta.gridLTS), length.out=75)
# likelihood.LTS<-sapply(theta.gridLTS2, FUN=likelihood.estimate)
# lines(theta.gridLTS2, likelihood.LTS, type='l', col=5)
# lines(pdfThetaLTS[,1], pdfThetaLTS[,2], col=5, lty=2)






###############################################
#Read in the posterior samples from the t distribution with 3 (or 5) degrees of freedom
###############################################

# #load("fittingTdistributionWorkspace.RData")
# #for 3 df
# if(nu==3){
# load(paste(tDistDirectory,"fittingTdistributionWorkspaceWithDelAnalysis", paramSet, ".RData", sep=''))
# }

# for 5 df
if(nu==5){
 # load(paste(tDistDirectory,"fittingTdistributionWorkspaceWithDelAnalysis", paramSet, ".RData", sep=''))
  load(paste(tDistDirectory,'workSpaceFittingT5dfParamsetOnGrid', paramSet, ".RData", sep=''))
  joint.postT<-log.post.density.final
  pdfThetaT<-cbind(theta.grid,post.theta.pdf)
  pdfSigma2T<-cbind(sigma2.grid,post.sigma2.pdf)
  delta.thetaT<-delta.theta
  delta.sigma2T<-delta.sigma2
  theta.gridT<-theta.grid
  sigma2.gridT<-sigma2.grid
  

  #load(do.call(paste, list("workSpaceFittingT5dfParamSet", paramSet, ".RData", sep='')))
# muSamplesT#<-mcmcTSamps[,1]
# sigma2SamplesT#<-mcmcTSamps[,2]
# 
# muSamplesTDel#<-mcmcTDelSamps[,1]
# sigma2SamplesTDel#<-mcmcTDelSamps[,2]
}


mean(pdfSigma2T)
likelihood.estimate<-function(theta){
  sigma2=var(newcomb)*(nu/(nu-2))
  kDensityEst<-sum(dt((newcomb-theta)/sqrt(sigma2), df=nu, log=TRUE)) #-.5*log(sigma2))
  return(exp(kDensityEst))
}
# likelihood.T<-sapply(theta.gridT, FUN=likelihood.estimate)
# plot(theta.gridT, log(likelihood.T), type='l', col=8)
# lines(theta.gridLTS2, log(likelihood.LTS), type='l', col=3)
# lines(pdfThetaT[,1], pdfThetaT[,2], col=4, lty=4)

#----------
# Figures
#----------

par(mfrow=c(1,1))
pdf(paste(path.to.figures,'posteriorBetaNewcombDataParamSet', paramSet,'trim_', trim, '.pdf', sep=''))
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(density(muSamplesFull), type='l', lty=1, lwd=4,ylim=c(0,.8), xlim=c(22,30),main=expression(paste('Marginal Posteriors for ', beta)), 
     xlab=expression(beta), col=1)
#Deleted obs lines(density(muSamplesFullDel), col='red', lwd=4)
#lines(density(muSamplesT), type='l',lwd=4,col=2)
lines(pdfThetaT[,1],pdfThetaT[,2] ,type='l',col=2, lty=2, lwd=4)

#lines(pdfThetaLMS[,1],pdfThetaLMS[,2], type='l',lty=1,lwd=4,col=3)
sp<-spline(pdfThetaLMS[,1],pdfThetaLMS[,2], n=50)
lines(sp$x, sp$y, type='l',col=3, lty=2, lwd=4)

#lines(pdfThetaLTS[,1],pdfThetaLTS[,2], type='l',lty=2,lwd=4,col=3)

sp<-spline(pdfThetaLTS[,1],pdfThetaLTS[,2], n=50)
lines(sp$x, sp$y,col=3, lty=3, lwd=4)

#lines(theta.grid,post.theta.pdf)
lines(pdfThetaRestricted[,1],pdfThetaRestricted[,2],lty=1,lwd=4,col=4)
# lines(density(muSamplesHuber), type='l',lty=1,lwd=4,col=4)
lines(density(muSamplesTukey), type='l',lwd=4,col='darkblue', lty=4)
#,"Restricted-Tukey")
legend('topleft',legend=c("Normal", 't','LMS','LTS', "Huber", "Tukey"),title='Model', cex=1.5,       col=c(1,2,3,3,4,'darkblue'), lwd=4, lty=c(1,2,2,3,1,4))
par(mar=c(5.1, 4.1, 4.1, 2.1))
#axis(1, at=mean(newcomb), col=1, labels=expression(bar(X)), las=1, lwd=4,tck=.05)
#axis(1, at=mean(sort(newcomb)[-c(1,2)]), col=2, labels=expression(bar(X)-del), las=1, lwd=4,tck=.05)

# legend('topleft',legend=c("Full","K=2", "K=4","K=8", "Huber"), cex=.75, 
#        col=c(1,2,3,4,1),lty=c(1,1,1,1,2), lwd=c(1,1,1,1,1))
# axis(1, at=mean(newcomb), col=1, labels=expression(bar(X)), las=1, lwd=4,tck=.05)
#axis(1, at=33.02, col=1, labels='truth', las=2, lwd=4,tck=.05)

dev.off()


par(mfrow=c(1,1))
pdf(paste(path.to.figures,'posteriorSigma2NewcombDataParamSet', paramSet,'trim_', trim, '.pdf', sep=''))
plot(density(sigma2SamplesFull), type='l', lty=1, lwd=4,ylim=c(0,.15), xlim=c(0,170),main=expression(paste('Marginal Posteriors for ', sigma^2)), 
     xlab=expression(sigma^2), col=1)
#lines(density((nu/(nu-2))*sigma2SamplesT, adjust=2), type='l',lwd=4,col=2)
lines((nu/(nu-2))*pdfSigma2T[,1],pdfSigma2T[,2]/(nu/(nu-2)), type='l',lwd=4,col=2, lty=2)

#lines(pdfSigma2LMS[,1],pdfSigma2LMS[,2], type='l',lty=2,lwd=4,col=3)
sp<-spline(pdfSigma2LMS[,1],pdfSigma2LMS[,2], n=45)
lines(sp$x, sp$y, type='l',lty=2,lwd=4,col=3)
lines(pdfSigma2LTS[,1],pdfSigma2LTS[,2], type='l',lty=3,lwd=4,col=3)
#lines(density(sigma2SamplesT), type='l',lwd=4,col=2)
#lines(density((nu/(nu-2))*mcmcTSamps[,2], adjust=2), type='l',lwd=4,col=2)
#lines(pdfSigma2Restricted[,1],pdfSigma2Restricted[,2], type='l',lty=2,lwd=4,col=4)
lines(density(sigma2SamplesHuber), type='l',lty=1,lwd=4,col=4)
lines(density(sigma2SamplesTukey), type='l',lwd=4,col='darkblue',lty=4)
legend('topright',legend=c("Normal", 't','LMS','LTS', "Huber", "Tukey"),title='Model', cex=1.5,       col=c(1,2,3,3,4,'darkblue'), lwd=4, lty=c(1,2,2,3,1,4))

dev.off()
# plot(muSamplesHuber[1:10000],muSamplesT[1:10000], pch=19, cex=.2)
# abline(0,1)




#####################################
#Summarizing Predictive distributions
#####################################



#6 models
#full model (fit using mcmc)
#t-model    (fit using mcmc)
#restricted model-Huber's (fit using direct sampling mcmc and importance sampling re-weighting)
#restricted model-Tukey's (fit using mcmc)
#LMS fit using direct sampling
#LTS fit using direct sampling 


#function to compute the predicted distrubution on a grid of y values given posterior samples for the models with a normal as the base
fn.pred.dist.normal.base<-function(ygrid, muSamples, sigma2Samples){
  predDist<-sapply(ygrid, FUN=function(x) mean(dnorm(x, muSamples, sqrt(sigma2Samples))))
  return(cbind(ygrid,predDist))
}



#function to compute the predicted distrubution on a grid of yvalues given evaluation of posterior on a grid for the models with a normal as the base
fn.pred.dist.normal.base2<-function(ygrid,joint.post.density,delta.theta, delta.sigma2){
  #joint.post.density: n.grid.points by 3 matrix. Column one theta, column 2 sigma2, column 3 posterior value
  #it is assumed the grid points are eaually spaced according to delta.theta and delta.sigma2
  theta<- joint.post.density[,1]
  sigma2<- joint.post.density[,2]
  post<- joint.post.density[,3]
  predDist<-sapply(ygrid, FUN=function(x) sum(dnorm(x, theta, sqrt(sigma2))*post)*delta.theta*delta.sigma2)
  return(cbind(ygrid,predDist))
}

# ygrid<-seq(15,45, length.out=50)
# tst<-fn.pred.dist.normal.base2(ygrid,post.density.final,delta.theta, delta.sigma2)


#function to compute the predicted distrubution on a grid of y vales give posterior samples for the t model with nu df (nu is a global parameter)
tpdf<-function(x, mu, sigma2t){
  gamma(.5*(nu+1))/(gamma(.5*nu)*sqrt(pi*nu*sigma2t))*(1+(1/nu)*((x-mu)^2/sigma2t))^(-(.5*(nu+1)))
}
#estimate using samples
fn.pred.dist.t.base<-function(ygrid, muSamples, sigma2Samples){
  predDist<-sapply(ygrid, FUN=function(x) {mean(dt(((x-muSamples)/sqrt(sigma2Samples)),df=nu)/sqrt(sigma2Samples))})
  return(cbind(ygrid,predDist))
}

#estimate from a grid fit
fn.pred.dist.t.base2<-function(ygrid,joint.post.density,delta.theta, delta.sigma2){
  #joint.post.density: n.grid.points by 3 matrix. Column one theta, column 2 sigma2, column 3 posterior value
  #it is assumed the grid points are eaually spaced according to delta.theta and delta.sigma2
  theta<- joint.post.density[,1]
  sigma2<- joint.post.density[,2]
  post<- joint.post.density[,3]
  predDist<-sapply(ygrid, FUN=function(x) sum(tpdf(x, theta, sigma2)*post)*delta.theta*delta.sigma2)
 # predDist<-sapply(ygrid, FUN=function(x) sum(dt(((x-theta)/sqrt(sigma2)),df=nu)*post*delta.theta*delta.sigma2/sqrt(sigma2)))
  return(cbind(ygrid,predDist))
}


#fit predictive distribution on a grid of y
range(newcomb)
deltay<-.01
ygrid<-seq(-20, 80, by=deltay)
ygrid2<-seq(-20, 80, length.out=1000)
deltay2<-ygrid2[2]-ygrid2[1]

system.time(fullPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesFull[1:1000], sigma2SamplesFull[1:1000]))
restHuberPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesHuber[1:1000], sigma2SamplesHuber[1:1000])

restHuberPredictiveDist2<-fn.pred.dist.normal.base2(ygrid2,joint.post.Huber,delta.thetaHuber, delta.sigma2Huber)
#check for match
plot(restHuberPredictiveDist[,1], restHuberPredictiveDist[,2], type='l')
lines(restHuberPredictiveDist2[,1], restHuberPredictiveDist2[,2], col=2)

restTukeyPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesTukey[1:1000], sigma2SamplesTukey[1:1000])

#tPredictiveDist<-fn.pred.dist.t.base(ygrid, muSamplesT[1:1000],sigma2SamplesT[1:1000])
tPredictiveDist2<-fn.pred.dist.t.base2(ygrid2, joint.postT,delta.thetaT,delta.sigma2T)

#check for match
# plot(tPredictiveDist[,1],tPredictiveDist[,2], type='l')
# lines(tPredictiveDist2[,1],tPredictiveDist2[,2], col=2, type='l')

#all.equal(tPredictiveDist, tPredictiveDist22)


#LMS and LTS
lmsPredictiveDist<-fn.pred.dist.normal.base2(ygrid2,joint.post.LMS,delta.thetaLMS, delta.sigma2LMS)

ltsPredictiveDist<-fn.pred.dist.normal.base2(ygrid2,joint.post.LTS,delta.thetaLTS, delta.sigma2LTS)


# fullPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesFull, sigma2SamplesFull)
# restHuberPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesHuber, sigma2SamplesHuber)
# restTukeyPredictiveDist<-fn.pred.dist.normal.base(ygrid, muSamplesTukey, sigma2SamplesTukey)
# tPredictiveDist<-fn.pred.dist.t.base(ygrid, muSamplesT, sigma2SamplesT)




#should all be 1
sum(fullPredictiveDist[,2])*deltay
sum(restHuberPredictiveDist[,2])*deltay
sum(restHuberPredictiveDist2[,2])*deltay2
sum(restTukeyPredictiveDist[,2])*deltay
sum(tPredictiveDist2[,2])*deltay2
# sum(tPredictiveDist_2[,2])*deltay
sum(lmsPredictiveDist[,2])*deltay2
sum(ltsPredictiveDist[,2])*deltay2
######################
#Calculate the mean and variance of the predictive distributions
#####################



pdf(paste(path.to.figures, "predDistParamset", paramSet, "nu", nu,'trim_', trim, '.pdf', sep=''))
plot(fullPredictiveDist[,1], fullPredictiveDist[,2], type='l', xlim=c(-10,60), ylim=c(0,.13), ylab='Density', xlab=expression(y),lwd=4)
lines(tPredictiveDist2[,1],tPredictiveDist2[,2], col=2,lwd=4, lty=2)
lines(lmsPredictiveDist[,1],lmsPredictiveDist[,2], col=3,lty=2, lwd=4)
lines(ltsPredictiveDist[,1],ltsPredictiveDist[,2], col=3,lty=3, lwd=4)
lines(restHuberPredictiveDist[,1],restHuberPredictiveDist[,2], col=4, lwd=4)
#lines(restHuberPredictiveDist2[,1],restHuberPredictiveDist2[,2], col=1, lwd=4)
lines(restTukeyPredictiveDist[,1],restTukeyPredictiveDist[,2], col='darkblue', lty=4, lwd=4)
legend('topleft',legend=c("Normal", 't','LMS','LTS', "Huber", "Tukey"),title='Model', cex=1.3,       col=c(1,2,3,3,4,'darkblue'), lwd=4, lty=c(1,2,2,3,1,4))

dev.off()

pdf(paste(path.to.figures, "logpredDistParamset", paramSet, "nu", nu,'trim_',trim, '.pdf', sep=''))
plot(fullPredictiveDist[,1], log(fullPredictiveDist[,2]), type='l', xlim=c(-10,60),ylim=c(-12, 1), ylab='log density', xlab=expression(y),lwd=4)
lines(tPredictiveDist2[,1],log(tPredictiveDist2[,2]), col=2,lwd=4, lty=2)
lines(lmsPredictiveDist[,1],log(lmsPredictiveDist[,2]), col=3,lty=2, lwd=4)
lines(ltsPredictiveDist[,1],log(ltsPredictiveDist[,2]), col=3,lty=3, lwd=4)
lines(restHuberPredictiveDist[,1],log(restHuberPredictiveDist[,2]), col=4, lwd=4)
#lines(restHuberPredictiveDist2[,1],log(restHuberPredictiveDist2[,2]), col=1, lwd=4)
lines(restTukeyPredictiveDist[,1],log(restTukeyPredictiveDist[,2]), col='darkblue', lty=4, lwd=4)
legend('topleft',legend=c("Normal", 't','LMS','LTS', "Huber", "Tukey"),title='Model', cex=1.3,       col=c(1,2,3,3,4,'darkblue'), lwd=4, lty=c(1,2,2,3,1,4))

dev.off()


# #pdf('newcomb.pdf')
# plot(density(newcomb), xlab='',main='', lwd=4)
# #dev.off()


meanPredFull<-sum(fullPredictiveDist[,1]*fullPredictiveDist[,2])*deltay
meanPredFull
varPredFull<-sum((fullPredictiveDist[,1]-meanPredFull)^2*fullPredictiveDist[,2])*deltay
absPredFull<-sum(abs(fullPredictiveDist[,1]-meanPredFull)*fullPredictiveDist[,2])*deltay

meanPredHuber<-sum(restHuberPredictiveDist[,1]*restHuberPredictiveDist[,2])*deltay
meanPredHuber
varPredHuber<-sum((restHuberPredictiveDist[,1]-meanPredHuber)^2*restHuberPredictiveDist[,2])*deltay
sqrt(varPredHuber)
absPredHuber<-sum(abs(restHuberPredictiveDist[,1]-meanPredHuber)*restHuberPredictiveDist[,2])*deltay

meanPredTukey<-sum(restTukeyPredictiveDist[,1]*restTukeyPredictiveDist[,2])*deltay
meanPredTukey
varPredTukey<-sum((restTukeyPredictiveDist[,1]-meanPredTukey)^2*restTukeyPredictiveDist[,2])*deltay
sqrt(varPredTukey)
absPredTukey<-sum(abs(restTukeyPredictiveDist[,1]-meanPredTukey)*restTukeyPredictiveDist[,2])*deltay

meanPredT<-sum(tPredictiveDist[,1]*tPredictiveDist[,2])*deltay
meanPredT
varPredT<-sum((tPredictiveDist[,1]-meanPredT)^2*tPredictiveDist[,2])*deltay
absPredT<-sum(abs(tPredictiveDist[,1]-meanPredT)*tPredictiveDist[,2])*deltay
varPredT
sqrt(varPredT)

meanPredLMS<-sum(lmsPredictiveDist[,1]*lmsPredictiveDist[,2])*deltay
meanPredLMS
varPredLMS<-sum((lmsPredictiveDist[,1]-meanPredLMS)^2*lmsPredictiveDist[,2])*deltay
varPredLMS


meanPredLTS<-sum(ltsPredictiveDist[,1]*ltsPredictiveDist[,2])*deltay
meanPredLTS
varPredLTS<-sum((ltsPredictiveDist[,1]-meanPredLTS)^2*ltsPredictiveDist[,2])*deltay
varPredLTS


meanPredFull
meanPredHuber
meanPredTukey
meanPredT
meanPredLMS
meanPredLTS

varPredFull
varPredHuber
varPredTukey
varPredT
varPredLMS
varPredLTS


