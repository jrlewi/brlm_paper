#
# Post processing the simulation results -----
#

rm(list=ls())
st<-Sys.time()
#
# Load Simulation Results ----
#
concTau2Near1<-0 #1 then load simulation with prior on tau2 concentrated at 1, if 0, then load the simulation with noninformative prior

if(concTau2Near1){
  load('simulationWorkspace_concTau2Near1.RData')
} else {
load('simulationWorkspace.RData')
}





# sigma2True: vector of the values of sigma2 used to generate the data
# thetaTrueSig2_sig2True: vector of true theta_i's for each of the 90 groups under that data generated with true sigma2 equal to sig2True

# normalTheorySig2_sig2True, restrictedFitHuberSig2_sigTrue, restrictedFitTukeySig2_sigTrue: objects with MCMC samples from normal theory fit, restricted with Huber/Huber and restriced with Tukey/Huber fit to data generated using the true value of sigma2 of sig2True

ngridpts<-5000 #for K-L



#
# Compute Evaluation Metrics ----
#

#
# Evaluation Metric 1 ----
#
# (theta_i.hat-theta_i)^2



for(sig2True in sigma2True){
  
  #get true theta_i
  thetaTrue<-get(paste0('thetaTrueSig2_', sig2True))
  
  #full model Fit
  fit<-get(paste0('normalTheoryFitSig2_', sig2True))
  postMn<-apply(fit$theta, 2, mean)
  mse<-(postMn-thetaTrue)^2
  assign(paste0('mseFullSig2_',sig2True), mse)
  
  
  #Restricted Huber model Fit
  fit<-get(paste0('restrictedFitHuberSig2_', sig2True))
  postMn<-apply(fit$theta, 2, mean)
  mse<-(postMn-thetaTrue)^2
  assign(paste0('mseRestHuberSig2_',sig2True),  mse)
 
  
  #RLM Huber Fits
  ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
  mse<-(ests-thetaTrue)^2
  assign(paste0('mseRlmHuberSig2_',sig2True),   mse)

  
  #Restricted Tukey model Fit
  fit<-get(paste0('restrictedFitTukeySig2_', sig2True))
  postMn<-apply(fit$theta, 2, mean)
  mse<-(postMn-thetaTrue)^2
  assign(paste0('mseRestTukeySig2_',sig2True),  mse)

  #RLM Tukey Fits
  ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
  mse<-(ests-thetaTrue)^2
  assign(paste0('mseRlmTukeySig2_',sig2True),   mse)
}

objKeep1<-c(paste0('mseFullSig2_',sigma2True),
            paste0('mseRestHuberSig2_',sigma2True),
            paste0('mseRlmHuberSig2_',sigma2True),
            paste0('mseRestTukeySig2_',sigma2True),
            paste0('mseRlmTukeySig2_',sigma2True))



#
# Evaluation Metric 2 ----
#
# K-L Divergence E[log f(y_good_i|theta_i, sig2True)/ f(y_good_i) ], expected value taken with respect to  f(y_good_i|theta_i, sig2True)

#evaluate on a grid of y
# length.out<-500
# ymin<--30
# ymax<-abs(ymin)
# ygrid<-seq(ymin, ymax, length.out=length.out)


# function to compute the pridictive distribution on a grid within each group

fn.compute.pred.dist<-function(ygrid, thetaSamples, sigma2Samples){
  #ygrid: grid of y values to evaluate the pred distribution
  #thetaSamples, sigma2Samples MCMC samples from the given group
  #OR for the robust regressions, these are just the estimates of theta and sigma2 from the robust regressions
   if(length(ygrid)>1000){
  sapply(ygrid, FUN=function(x) {
    mean(dnorm(x, thetaSamples, sqrt(sigma2Samples)))
    } 
    ) } else {
       ygridN<-matrix(ygrid, nrow=length(thetaSamples), ncol=length(ygrid), byrow = TRUE)
    colMeans(dnorm(ygridN, thetaSamples, sqrt(sigma2Samples)))
     }
}


# lg<-1000
# ygrid1=seq(-10,10, length.out=lg)
# fit<-restrictedFitHuberSig2_0.5
# 
# system.time(tst1<-fn.compute.pred.dist(ygrid1, thetaSamples=fit$theta[,1], sigma2Samples=fit$sigma2[,1]))
# 
# lg<-1001
# ygrid2=seq(-10,10, length.out=lg)
# fit<-restrictedFitHuberSig2_0.5
# 
# system.time(tst2<-fn.compute.pred.dist(ygrid2, thetaSamples=fit$theta[,1], sigma2Samples=fit$sigma2[,1]))
# 
# plot(ygrid1,tst1, type='l')
# lines(ygrid2,tst2)


# all.equal(fn.compute.pred.dist(ygrid=c(1,2,3), thetaSamples=c(1), sigma2Samples=c(2)), dnorm(c(1,2,3), 1, sqrt(2)))
# tst<-fn.compute.pred.dist(ygrid, normalTheoryFitSig2_0.5$theta, normalTheoryFitSig2_0.5$sigma2)
# plot(ygrid, tst, type='l')
     
# function to compute the K-L metric, for one group

fn.compute.KL<-function(thetaTrue,thetaSamples, sigma2Samples, sig2True, ngridpts){
  #thetaTrue: scalar of the true value of theta for the given group
  #thetaSamples: MCMC sample of the theta_i from the given group
      #OR for robust regression it is just the estimate of theta for the given group
  #sigma2Sample: MCMC sample of the sigma2_i from the given group
    #OR for robust regression it is just the estimate of sigma2 for the given group
  #sig2True: scalar of the true value of sigma2
  #ygrid: grid of y values for the quadrature approximation of the expectation: not used anymore
  #ngridpts: number of grid points
  
  ygrid<-seq(thetaTrue-10*sqrt(sig2True),thetaTrue+10*sqrt(sig2True), length.out=ngridpts)
  predDist<-fn.compute.pred.dist(ygrid,  thetaSamples, sigma2Samples)
  trueDist<-dnorm(ygrid, thetaTrue, sqrt(sig2True))
  dif<-diff(ygrid)
  kl<-((log(trueDist)-log(predDist))*trueDist)[-1]
  KL<-sum((kl*dif)[kl!=Inf])
#   if(plot){
#     plot(ygrid, trueDist, type='l')
#     lines(ygrid, predDist)
#   }
  #out<-list()
  #out$ygrid<-ygrid
  #out$predDist<-predDist
  #out$KL<-KL
  #out
  KL
}

#plots

# ngridpts<-1000
# sig2True<-4
# thetaTrue<-get(paste0('thetaTrueSig2_', sig2True))
# mcmcFit<-get(paste0("restrictedFitHuberSig2_", sig2True))
# i=59
# factorsMat[i,]
# ygrid<-seq(thetaTrue[i]-10*sqrt(sig2True),thetaTrue[i]+10*sqrt(sig2True), length.out=ngridpts)
# pred1<-fn.compute.pred.dist(ygrid, thetaSamples=mcmcFit$theta[,i], sigma2Samples=mcmcFit$sigma2[,i])
# trueDist<-dnorm(ygrid, thetaTrue[i], sqrt(sig2True))
# 
# bhat<-sapply(mcmcFit$robustFits, FUN=function(x) x$coefficients)
# shat<-sapply(mcmcFit$robustFits, FUN=function(x) x$s)
# pred2<-dnorm(ygrid,bhat[i], shat[i] )
# 
# plot(ygrid, trueDist, col=4, type='l', ylim=c(0, max(trueDist, pred1)))
# lines(ygrid, pred1)
# lines(ygrid, pred2, col=2)
# sum(((log(trueDist)-log(pred1))*trueDist)[-1]*diff(ygrid))
# sum(((log(trueDist)-log(pred2))*trueDist)[-1]*diff(ygrid))



#exact: for use with the robust regression methods
fn.KL.twoNormals<-function(mu1,sigma2_1, mu2, sigma2_2){
  #KL(p,q) where p~N(mu_1,sigma2_1),q~N(mu_2, sigma2_2)
  sig_1<-sqrt(sigma2_1)
  sig_2<-sqrt(sigma2_2)
  log(sig_2/sig_1)+(sigma2_1+(mu1-mu2)^2)/(2*sigma2_2)-1/2  
}


# fn.compute.KL(0, normalTheoryFitSig2_0.5$theta[,1], normalTheoryFitSig2_0.5$sigma2[,1],1, ygrid)

# Function to compute the KL for each group

fn.compute.KL.each.group<-function(ngridpts, thetaTrueVect, sig2True,thetaSamplesMat, sigma2SamplesMat){
  # thetaTrueVect: the vector of true theta values for each group
  # thetaSampleMat: nTot by nGroups matrix of MCMC samples of the theta_i's
    #OR for robust regression a vector of estimates of the theta_is
  # sigma2SamplesMat: nTot by nGroups matrix of MCMC samples of the sigma2_i's
    #OR for robust regression a vector of estimates of the sigma2_is
  if(class(thetaSamplesMat)=='matrix'){
  if(class(sigma2SamplesMat)!='matrix'){ stop('error: thetaSamplesMat and sigma2SamplesMat must both be matrices of MCMC sampples or both be vectors of robust regression estimates')}  
  #create a list out of the columns of thetaSamplesMat and sigma2SamplesMat
  #this will be fed into mapply
  thetaSamplesMat<-split(thetaSamplesMat, col(thetaSamplesMat))
  sigma2SamplesMat<-split(sigma2SamplesMat, col(sigma2SamplesMat))
  } else {
    if(class(thetaSamplesMat)!='numeric' || class(sigma2SamplesMat)!='numeric'){
      stop('error: thetaSamplesMat and sigma2SamplesMat must both be matrices of MCMC sampples or both be vectors of robust regression estimates')
    }
    } 
mapply(fn.compute.KL,thetaTrue=thetaTrueVect,thetaSamples=thetaSamplesMat, sigma2Samples=sigma2SamplesMat, MoreArgs=list(sig2True=sig2True, ngridpts=ngridpts))
}


# 
# # # mini test
# ngridpts<-10
# sig2True<-.5
# thetaTrue<-get(paste0('thetaTrueSig2_', sig2True))
# mcmcFit<-get(paste0("restrictedFitHuberSig2_", sig2True))

# system.time(tst1<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True,thetaSamplesMat=mcmcFit$theta, sigma2SamplesMat=mcmcFit$sigma2))
# 
# st<-Sys.time()
# tst2<-numeric(90)
# for(i in 1:90){
#  tst2[i]<-fn.compute.KL(thetaTrue[i],thetaSamples=mcmcFit$theta[,i], sigma2Samples=mcmcFit$sigma2[,i], sig2True,ngridpts)
# }
# end<-Sys.time()-st
# end
# 
# 
# all.equal(tst1,tst2)
# #mini test for robust reg estimates
# bhat<-sapply(mcmcFit$robustFits, FUN=function(x) x$coefficients)
# shat2<-sapply(mcmcFit$robustFits, FUN=function(x) x$s^2)
# ngridpts<-1e2
# system.time(tst1Rlm<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True,thetaSamplesMat=bhat, sigma2SamplesMat=shat2))
# st<-Sys.time()
# tst2Rlm<-numeric(90)
# for(i in 1:90){
#   tst2Rlm[i]<-fn.compute.KL(thetaTrue[i],thetaSamples=bhat[i], sigma2Samples=shat2[i], sig2True,ngridpts)
# }
# end<-Sys.time()-st
# end
# 
# exactKLRlm<-numeric(90)
# for(i in 1:90){
#   exactKLRlm[i]<-fn.KL.twoNormals(mu1=thetaTrue[i],sigma2_1=sig2True, mu2=bhat[i], sigma2_2=shat2[i])
# }
# 
# 
# all.equal(tst1,tst2)
# 
# all.equal(tst1Rlm,tst2Rlm)
# 
# # plot(tst1, tst1Rlm, pch=16)
# # abline(0,1, col=4, lwd=2)
# 
# plot(exactKLRlm, tst1Rlm, pch=16)
# abline(0,1, col=4, lwd=2)
# 
# all.equal(exactKLRlm, tst1Rlm)
# all.equal(round(exactKLRlm,10), round(tst1Rlm,10))
# 
# mean(tst1)
# mean(tst1Rlm)
# plot(apply(mcmcFit$sigma2,2, mean)^.5, shat2^.5, pch=16)
# abline(0,1, lwd=2, col=4)

# i=55
# pred0<-fn.compute.pred.dist(ngridpts, bhat[i],shat2[i])
# 
# logpred<-dnorm(ngridpts, bhat[i],sqrt(shat2[i]), log=TRUE)
# pred<-dnorm(ngridpts, bhat[i],sqrt(shat2[i]))
# all.equal(pred0, pred)
# all.equal(pred0, exp(logpred))
# ftrue<-dnorm(ngridpts, thetaTrue[i], sqrt(sig2True))
# kl<-((log(ftrue)-logpred)*ftrue)[-1]
# kl2<-((log(ftrue)-log(pred))*ftrue)[-1]
# sum((kl*diff(ngridpts)), na.rm=TRUE)
# sum((kl2*diff(ngridpts)), na.rm=TRUE)
# sum((kl2*diff(ngridpts))[kl2!=Inf], na.rm=TRUE)
# tst1Rlm[i]
# exactKLRlm[i]


# Note: some numerical issues occur for the plug in predictive distribution for the robust regression estimates. The pred dist at the tail ygrid values is numerically zero. Could compute log pred for the robsust estimates. For unified code function fn.compute.pred to handle both MCMC samples and plug in, I handled this in compute.KL by adding KL<-sum((kl*dif)[kl!=Inf]). Both ways work here, but it is a note of caution. 


for(sig2True in sigma2True){
  print(sig2True)
  #get true theta_is
  thetaTrue<-get(paste0('thetaTrueSig2_', sig2True))
  
  #full model Fit KL
#   st<-Sys.time()
  fit<-get(paste0('normalTheoryFitSig2_', sig2True))
  kl<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True=sig2True,thetaSamplesMat=fit$theta, sigma2SamplesMat=fit$sigma2)
  assign(paste0('klFullSig2_',sig2True), kl)
#   Sys.time()-st
  
  #Restricted Huber model Fit KL
  fit<-get(paste0('restrictedFitHuberSig2_', sig2True))
  kl<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True=sig2True,thetaSamplesMat=fit$theta, sigma2SamplesMat=fit$sigma2)
  assign(paste0('klRestHuberSig2_',sig2True), kl)
  mean(kl)
  
  #RLM KL Huber 
  ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
  shat2<-sapply(fit$robustFits, FUN=function(x) x$s^2)
#   kl2<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True=sig2True,thetaSamplesMat=ests, sigma2SamplesMat=shat2)
  kl<-numeric(90)
  for(i in 1:90){
    kl[i]<-fn.KL.twoNormals(mu1=thetaTrue[i],sigma2_1=sig2True, mu2=ests[i], sigma2_2=shat2[i])
  } 
  
  assign(paste0('klRlmHuberSig2_',sig2True),   kl)
  mean(kl)

  #Restricted Tukey model Fit KL
  fit<-get(paste0('restrictedFitTukeySig2_', sig2True))
  kl<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True=sig2True,thetaSamplesMat=fit$theta, sigma2SamplesMat=fit$sigma2)
  assign(paste0('klRestTukeySig2_',sig2True), kl)


  #RLM KL Tukey
  ests<-sapply(fit$robustFits, FUN=function(x) x$coefficients)
  shat2<-sapply(fit$robustFits, FUN=function(x) x$s^2)
#   kl<-fn.compute.KL.each.group(ngridpts, thetaTrueVect=thetaTrue, sig2True=sig2True,thetaSamplesMat=ests, sigma2SamplesMat=shat2)
kl<-numeric(90)
for(i in 1:90){
  kl[i]<-fn.KL.twoNormals(mu1=thetaTrue[i],sigma2_1=sig2True, mu2=ests[i], sigma2_2=shat2[i])
} 
  assign(paste0('klRlmTukeySig2_',sig2True),   kl)
}


objKeep2<-c(paste0('klFullSig2_',sigma2True),
            paste0('klRestHuberSig2_',sigma2True),
            paste0('klRlmHuberSig2_',sigma2True),
            paste0('klRestTukeySig2_',sigma2True),
            paste0('klRlmTukeySig2_',sigma2True))


objKeep<-c(objKeep1, objKeep2, 'factorsMat',"sigma2True", 'concTau2Near1')
end<-Sys.time()-st
end

sdiff<-setdiff(ls(), objKeep)
rm(list=sdiff) #  everything else can be recovered by loading simulationWorkspace.RData.

if(concTau2Near1){
  save.image('simPostProcessed_concTau2Near1.RData')
} else {
  save.image('simPostProcessed.RData')
}


