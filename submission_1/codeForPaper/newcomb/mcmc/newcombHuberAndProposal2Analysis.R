########################################
#John Lewis
#Pre-lim analysis of newcomb's data to see if the calcs are correct
#########################################

strt<-Sys.time()

set.seed(1)

library(MCMCpack)

#number of samples
nkeep<-1e6
nburn<-1e5

paramSet<-2  #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper
source(do.call(paste, list('../',"newcombFixedParameters", paramSet, '.R', sep='')))


#Y_1,...,Y_n iid N(mu, sigma^2)
#let T(Y), S(Y) be Huber estimators of location and scale

#Prior Distributions for mu and sigma^2 are specified as follows:
#mu~N(eta, tau^2)
#sigma^2~IG(alpha, beta)
#Analysis is based on the posterior [mu, sigma^2|T(Y),S(Y)]



#sampling functions

sample.mu<-function(y,sigma.sq,eta, tau)
{
  condMean<-((n*mean(y)*tau^2+eta*sigma.sq)/(n*tau^2+sigma.sq))
  condVar<-(sigma.sq*tau^2)/(n*tau^2+sigma.sq)
  return(rnorm(1,mean=condMean, sd=sqrt(condVar)))
}


sample.sigma.sq<-function(y, mu, alpha, beta, n=length(y))
{
  alphaNew<-(0.5*n+alpha)
  betaNew<-(0.5*sum((y-mu)^2)+beta)
  return(rinvgamma(1,shape= alphaNew,scale=betaNew))
}

#setwd("~/My Box Files/research/snm/locationAndScale/mcmcCodes")
# setwd("../../../mcmcCodes")
source('locAndScaleHuberAndProposal2.R')



total<-nkeep+nburn
muSamples<-numeric(total)
sigma2Samples<-numeric(total)

fit1<-rlm(newcomb~1,scale.est='Huber')
l1obs<-coef(fit1)
s1obs<-fit1$s
sigma2Cur<-100
#note:: shouldn't start in the newcomb data because the radius is so much bigger because of the outliers
y.curr<-rnorm(n)
y.curr<-fn.compute.ystst(y.curr, l1obs=l1obs,s1obs=s1obs)

n<-length(newcomb)
accept<-0

###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(y.curr,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(y.curr, muCur, alpha=alpha, beta=beta, n=n)
  y.currAnda<-fn.one.rep.y(y.curr,muCur,sqrt(sigma2Cur),l1obs=l1obs, s1obs=s1obs,n=n)
  y.curr<-y.currAnda[[1]]
  muSamples[i]<-muCur
  sigma2Samples[i]<-sigma2Cur
  accept<-accept+y.currAnda[[2]]/total
  })
acceptHuber<-accept



muSamplesHuber<-muSamples[nburn:total]
sigma2SamplesHuber<-sigma2Samples[nburn:total]
#plot(muSamplesHuber,pch=19,cex=.2)
#plot(sigma2SamplesHuber,pch=19,cex=.2)


#####################################################################
#Using Data set with deleted outliers
####################################################################
#number of samples
#nkeep<-1e5
#nburn<-1000
total<-nkeep+nburn
muSamples<-numeric(total)
sigma2Samples<-numeric(total)
##################
#Delete outliers
#################
newcombDel<-sort(newcomb)[-c(1:2)]
n<-length(newcombDel)
fit1<-rlm(newcombDel~1,scale.est='Huber')
l1obs<-coef(fit1)
s1obs<-fit1$s
sigma2Cur<-100
#note:: shouldn't start in the newcomb data because the radius is so much bigger because of the outliers
y.curr<-rnorm(n)
y.curr<-fn.compute.ystst(y.curr, l1obs=l1obs,s1obs=s1obs)
accept<-0
#resource
source('locAndScaleHuberAndProposal2.R')
###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(y.curr,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(y.curr, muCur, alpha=alpha, beta=beta, n=n)
  y.currAnda<-fn.one.rep.y(y.curr,muCur,sqrt(sigma2Cur),l1obs=l1obs, s1obs=s1obs,n=n)
  y.curr<-y.currAnda[[1]]
  muSamples[i]<-muCur
  sigma2Samples[i]<-sigma2Cur
  accept<-accept+y.currAnda[[2]]/total
})
acceptHuberDel<-accept

muSamplesHuberDel<-muSamples[nburn:total]
sigma2SamplesHuberDel<-sigma2Samples[nburn:total]
#plot(muSamplesHuberDel,pch=19,cex=.2)
#plot(sigma2SamplesHuberDel,pch=19,cex=.2)




rm(list=setdiff(ls(), c("muSamplesHuber","sigma2SamplesHuber", 'acceptHuber',"muSamplesHuberDel","sigma2SamplesHuberDel", 'acceptHuberDel','paramSet')))
# #setwd("~/My Box Files/research/snm/locationAndScale/dataAnalysis/newcomb/mcmc")
setwd("../dataAnalysis/newcomb/mcmc")
save.image(do.call(paste, list('WorkSpaceNewcombHuberAndProposal2AnalysisWithDelAnalysis',paramSet, '.RData', sep='')))

# plot(muSamplesHuber)
# plot(muSamplesHuberDel)
# 
# plot(sigma2SamplesHuber)
# plot(sigma2SamplesHuberDel)

