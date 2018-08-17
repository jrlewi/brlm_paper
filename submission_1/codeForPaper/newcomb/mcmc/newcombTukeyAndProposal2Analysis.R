########################################
#John Lewis
#Analysis newcomb's data using restricted likelihood of Tukey's location and Tukey's proposal 2
#########################################

strt<-Sys.time()

set.seed(2)

library(MCMCpack)


#number of samples
nkeep<-1e4
nburn<-1e4


paramSet<-2  #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper
source(do.call(paste, list('../',"newcombFixedParameters", paramSet, '.R', sep='')))

eta
tau
alpha
beta
data(newcomb)
n<-length(newcomb)
y<-newcomb


# data(newcomb)
# n<-length(newcomb)

#Y_1,...,Y_n iid N(mu, sigma^2)
#let T(Y), S(Y) be Tukey estimators of location and scale

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
source('locAndScaleTukeyAndProposal2.R')



total<-nkeep+nburn
muSamplesTukey<-numeric(total)
sigma2SamplesTukey<-numeric(total)

fit1<-rlm(newcomb~1,psi=psi.bisquare, scale.est='Huber')
l1obs<-coef(fit1)
s1obs<-fit1$s
sigma2Cur<-100
#note: shouldn't start in the newcomb data because the radius is so much bigger because of the outliers
y.curr<-rnorm(n)
y.curr<-fn.compute.ystst.tukey(y.curr, l1obs=l1obs,s1obs=s1obs)

n<-length(newcomb)
accept.tukey<-0
###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(y.curr,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(y.curr, muCur, alpha=alpha, beta=beta, n=n)
  y.currAnda<-fn.one.rep.y.tukey(y.curr,muCur,sqrt(sigma2Cur),l1obs=l1obs, s1obs=s1obs,n=n)
  y.curr<-y.currAnda[[1]]
  muSamplesTukey[i]<-muCur
  sigma2SamplesTukey[i]<-sigma2Cur
  accept.tukey<-accept.tukey+y.currAnda[[2]]/total
})
accept.tukey 
#muSamplesTukey<-c(muSamplesTukey,muSamplesTukey2)
#sigma2SamplesTukey<-c(sigma2SamplesTukey, sigma2SamplesTukey2)
# plot(muSamplesTukey,pch=19,cex=.2)
#plot(density(muSamplesTukey))
#plot(sigma2SamplesTukey,pch=19,cex=.2)


muSamplesTukey<-muSamplesTukey[nburn:total]
sigma2SamplesTukey<-sigma2SamplesTukey[nburn:total]


#####################################################################
#Using Data set with deleted outliers
####################################################################
#number of samples
#nkeep<-1e5
#nburn<-1000
total<-nkeep+nburn
muSamplesTukeyDel<-numeric(total)
sigma2SamplesTukeyDel<-numeric(total)

##################
#Delete outliers
#################
newcombDel<-sort(newcomb)[-c(1:2)]
n<-length(newcombDel)



fit1<-rlm(newcombDel~1,psi=psi.bisquare, scale.est='Huber')
l1obs<-coef(fit1)
s1obs<-fit1$s
sigma2Cur<-100
#note: shouldn't start in the newcomb data because the radius is so much bigger because of the outliers
y.curr<-rnorm(n)
y.curr<-fn.compute.ystst.tukey(y.curr, l1obs=l1obs,s1obs=s1obs)

accept.tukeyDel<-0
#resource
source('locAndScaleHuberAndProposal2.R')



###############################
system.time(for(i in 1:total){
  muCur<-sample.mu(y.curr,sigma2Cur,eta, tau)
  sigma2Cur<-sample.sigma.sq(y.curr, muCur, alpha=alpha, beta=beta, n=n)
  y.currAnda<-fn.one.rep.y.tukey(y.curr,muCur,sqrt(sigma2Cur),l1obs=l1obs, s1obs=s1obs,n=n)
  y.curr<-y.currAnda[[1]]
  muSamplesTukeyDel[i]<-muCur
  sigma2SamplesTukeyDel[i]<-sigma2Cur
  accept.tukeyDel<-accept.tukeyDel+y.currAnda[[2]]/total
})
accept.tukeyDel 
#muSamplesTukey<-c(muSamplesTukey,muSamplesTukey2)
#sigma2SamplesTukey<-c(sigma2SamplesTukey, sigma2SamplesTukey2)
#plot(muSamplesTukey,pch=19,cex=.2)
#plot(sigma2SamplesTukey,pch=19,cex=.2)


muSamplesTukeyDel<-muSamplesTukeyDel[nburn:total]
sigma2SamplesTukeyDel<-sigma2SamplesTukeyDel[nburn:total]

ed<-Sys.time()-strt
ed


#Elapsed time for 50000 samples: 196.54 
rm(list=setdiff(ls(), c("muSamplesTukey","sigma2SamplesTukey", 'accept.tukey',"muSamplesTukeyDel","sigma2SamplesTukeyDel", 'accept.tukeyDel', 'paramSet')))

setwd("../dataAnalysis/newcomb/mcmc")
save.image(do.call(paste, list('WorkSpaceNewcombTukeyAndProposal2AnalysisWithDelAnalysis',paramSet, '.RData', sep='')))
# plot(muSamplesTukey)
# plot(muSamplesTukeyDel)
# plot(sigma2SamplesTukey)
# plot(sigma2SamplesTukeyDel)

