########################################
#John Lewis
#Fitting the full normal to newcombs data
#Custom Code
#########################################

strt<-Sys.time()

set.seed(2)

library(MCMCpack)


#number of samples
nkeep<-1e6
nburn<-1e5

paramSet<-2  #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper
source(do.call(paste, list('../',"newcombFixedParameters", paramSet, '.R', sep='')))

alpha
beta
eta
tau

#Y_1,...,Y_n iid N(mu, sigma^2)
#let T(Y), S(Y) be Huber estimators of location and scale

#Prior Distributions for mu and sigma^2 are specified as follows:
#mu~N(eta, tau^2)
#sigma^2~IG(alpha, beta)
#Analysis is based on the posterior [mu, sigma^2|T(Y),S(Y)]


# #Specify the prior Parameters
# #prior parameters #prior parameters for the mean are chosen by looking a
# #estimates of the speed of light from before the experiment.
# #these were determined in speedOflightAnalysisClean2 in the orderstatistics folder


#sampling functions

sample.mu<-function(y,sigma.sq,eta, tau, n=length(y))
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



total<-nkeep+nburn
muSamplesFull<-numeric(total)
sigma2SamplesFull<-numeric(total)
sigma2Cur<-100
n<-length(newcomb)

###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(newcomb,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(newcomb, muCur, alpha=alpha, beta=beta)
  muSamplesFull[i]<-muCur
  sigma2SamplesFull[i]<-sigma2Cur
})


#plot(muSamplesFull, pch=19,cex=.2)
#plot(sigma2SamplesFull, pch=19,cex=.2)

#plot(density(muSamplesFull))
muSamplesFull<-muSamplesFull[nburn:total]
sigma2SamplesFull<-sigma2SamplesFull[nburn:total]


#############################
#Deleting the lowest two observations
#############################
newcombDel<-sort(newcomb)[3:66]

#number of samples
#nkeep<-1e5
#nburn<-1000
total<-nkeep+nburn
muSamplesFullDel<-numeric(total)
sigma2SamplesFullDel<-numeric(total)
sigma2Cur<-100
n<-length(newcombDel)

###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(newcombDel,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(newcombDel, muCur, alpha=alpha, beta=beta)
  muSamplesFullDel[i]<-muCur
  sigma2SamplesFullDel[i]<-sigma2Cur
})


#plot(muSamplesFullDel, pch=19,cex=.2)
#plot(sigma2SamplesFullDel, pch=19,cex=.2)

#plot(density(muSamplesFullDel))


muSamplesFullDel<-muSamplesFullDel[nburn:total]
sigma2SamplesFullDel<-sigma2SamplesFullDel[nburn:total]




rm(list=setdiff(ls(), c("muSamplesFull","sigma2SamplesFull","muSamplesFullDel","sigma2SamplesFullDel", "paramSet")))

save.image(do.call(paste, list('WorkSpaceNewcombFullNormalModel',paramSet, '.RData', sep='')))

