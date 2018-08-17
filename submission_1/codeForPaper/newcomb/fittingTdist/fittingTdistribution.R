########################################
#John Lewis
#using a t-distribution in the analysis of the newcomb data
#using my own code
#########################################
set.seed(3)
library(MCMCpack)
library(msm)
library(MCMCpack)


#number of samples
nkeep<-1e6
nburn<-1e5


paramSet<-2  #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper
source(do.call(paste, list('../',"newcombFixedParameters", paramSet, '.R', sep='')))
#source('../newcombFixedParameters1.R')
eta
tau
alpha
beta
data(newcomb)
n<-length(newcomb)
y<-newcomb

#Y_1,...,Y_n iid t_nu(mu, sigma^2)

#Prior Distributions for mu and sigma^2 are specified as follows:
#mu~N(eta, tau^2)
#sigma^2~IG(alpha, beta)
#Analysis is based on the posterior [mu, sigma^2|Y]
#set up with auxillary parameters v_i: y_i~N(mu,V_i)
#V_i~inv-chi(nu,sigma^2)


#Specify the prior Parameters
#prior parameters #prior parameters for the mean are chosen by looking a
#estimates of the speed of light from before the experiment.
#these were determined in speedOflightAnalysisClean2 in the orderstatistics folder
##############################################################


#special fo t-distribution; keep
nu<-5 # degrees of freedom for the t
#should the prior on sigma2 change??
#yes
#want the prior on the var(Y)=(nu/(n-2))sigma2 ~IG(alpha, beta)
#so sigma2=((n-2)/nu)Var(Y)~IG(alpha,((n-2)/nu) beta)
beta<-beta*(nu-2)/nu 
#####################################################


#setwd("~/My Box Files/research/snm/locationAndScale/dataAnalysis/newcomb/fittingTdist")

# setwd("~/Box Sync/research/snm/BlendedParadigmWriting/codeForPaper/newcomb/fittingTdist")


#sampling functions
sample.v<-function(mu,sigma2){
  (nu*sigma2+(y-mu)^2)/rchisq(n,nu+1)
}
#v<-sample.v(30,30)
sample.mu<-function()
{
  a<-sum(y/v)
  b<-sum(1/v)
  condMean<-(eta/tau^2+a)/(1/tau^2+b)
  condVar<-(1/tau^2+b)^-1
  return(rnorm(1,mean=condMean, sd=sqrt(condVar)))
}



rwTune=5
logSigma2Target<-function(sigma2){
  log(sigma2)*(-alpha-1+0.5*n*nu)-beta/sigma2-.5*nu*sigma2*sum(1/v)
}

proposeSigma2Sample<-function(){
  rtnorm(1,mean=sigma2Cur,sd=rwTune, lower=0)
}

logSigma2Prop<-function(sigma2Prop,sigma2Cur){
  dtnorm(sigma2Prop,mean=sigma2Cur,sd=rwTune, lower=0, log=TRUE)
}

mhratioSigma2<-function(sigma2Cur,sigma2Prop){
  rat<-(logSigma2Target(sigma2Prop)-logSigma2Target(sigma2Cur))+(logSigma2Prop(sigma2Cur, sigma2Prop)-logSigma2Prop(sigma2Prop,sigma2Cur))
  return(min(1,exp(rat)))
}

mh.sigma2<-function(sigma2Cur)
{
 sigma2Prop<-proposeSigma2Sample()
 prop<-mhratioSigma2(sigma2Cur,sigma2Prop)
 if(runif(1)<prop){return(sigma2Prop)}
 else{return(sigma2Cur)}
 }
#mh.sigma2(sigma2Cur)



total<-nkeep+nburn
#sample vectors for the T distrbution model
muSamples<-numeric(total)
sigma2Samples<-numeric(total)
#vSamples<-matrix(total*length(y),total, length(y))

#initial values
muCur<-rnorm(1,eta,tau)
sigma2Cur<-rinvgamma(1,alpha,beta)
###############################
system.time(for(i in 1:total){
  v<-sample.v(muCur,sigma2Cur)
  muCur<-sample.mu()
  sigma2Cur<-mh.sigma2(sigma2Cur)
  muSamples[i]<-muCur
  sigma2Samples[i]<-sigma2Cur
  #vSamples[i,]<-v don't save the v's
  if(i%%1000==0) print(i)
})

#muSamplesT<-c(muSamplesT,muSamplesT2)
#sigma2SamplesT<-c(sigma2SamplesT, sigma2SamplesT2)

muSamplesT<-muSamples[nburn:total]
sigma2SamplesT<-sigma2Samples[nburn:total]
#vSamplesT<-vSamples[nburn:total,]


rm(list=setdiff(ls(), c("muSamplesT","sigma2SamplesT", 'nu',"paramSet")))
#don't keep "vSamplesT" because it is a masive matrix

#setwd("~/My Box Files/research/snm/locationAndScale/dataAnalysis/newcomb/fittingTdist")
save.image(do.call(paste, list('workSpaceFittingT',nu, 'dfParamset',paramSet, '.RData', sep='')))


# 
# 
# ###############################################
# #Using Data set with deleted outliers
# ###########################################
# 
# ##################
# #Delete outliers
# #################
# y<-sort(newcomb)[-c(1:2)]
# n<-length(y)
# 
# #sample vectors for the T distrbution model
# muSamples<-numeric(total)
# sigma2Samples<-numeric(total)
# vSamples<-matrix(total*length(y),total, length(y))
# 
# #initial values
# muCur<-rnorm(1,eta,tau)
# sigma2Cur<-rinvgamma(1,alpha,beta)
# ###############################
# system.time(for(i in 1:total){
#   v<-sample.v(muCur,sigma2Cur)
#   muCur<-sample.mu()
#   sigma2Cur<-mh.sigma2(sigma2Cur)
#   muSamples[i]<-muCur
#   sigma2Samples[i]<-sigma2Cur
#   vSamples[i,]<-v
# })
# 
# #muSamplesT<-c(muSamplesT,muSamplesT2)
# #sigma2SamplesT<-c(sigma2SamplesT, sigma2SamplesT2)
# 
# muSamplesTDel<-muSamples[nburn:total]
# sigma2SamplesTDel<-sigma2Samples[nburn:total]
# vSamplesTDel<-vSamples[nburn:total,]
# 
# 
# 
# rm(list=setdiff(ls(), c("muSamplesT","sigma2SamplesT","vSamplesT", 'nu',"muSamplesTDel","sigma2SamplesTDel","vSamplesTDel","paramSet")))
# 
# #setwd("~/My Box Files/research/snm/locationAndScale/dataAnalysis/newcomb/fittingTdist")
# save.image(do.call(paste, list('fittingTdistributionWorkspaceWithDelAnalysis',paramSet, '.RData', sep='')))
# # plot(muSamplesT)
# # plot(muSamplesTDel)
# # plot(sigma2SamplesT)
# # plot(sigma2SamplesTDel)
# # 
# # 
# # 
# # 
# 
# # muSamplesT1<-muSamplesT
# # sigma2SamplesT1<-sigma2SamplesT
# # 
# # plot(density(muSamplesT1))
# # lines(density(mcmcTSamps[,1]))
# # 
# # 
# # plot(density(sigma2SamplesT1))
# # lines(density(mcmcTSamps[,2]))