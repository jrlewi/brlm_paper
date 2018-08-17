########################################
#John Lewis
#using a t-distribution in the analysis of the newcomb data
#using Grid estimation
#########################################
set.seed(3)
library(MCMCpack)
library(msm)
library(MCMCpack)

#load this workspace along with fixed hyperparameters
load("~/Box Sync/research/snm/locationAndScale/dataAnalysis/newcomb/fittingTdist/workSpaceFittingT5dfParamset2.RData")
nu
eta
tau
alpha
beta #has already been adjusted for t-model
y

#Y_1,...,Y_n iid t_nu(mu, sigma^2)
#Prior Distributions for mu and sigma^2 are specified as follows:
#mu~N(eta, tau^2)
#sigma^2~IG(alpha, beta)
#Analysis is based on the posterior [mu, sigma^2|Y]
#set up with auxillary parameters v_i: y_i~N(mu,V_i)
#V_i~inv-chi(nu,sigma^2)

#define the grid
#delta.theta<-.02  #distance between points on theta grid 
ticks.theta<-200
theta.grid<-seq(20,35,length.out=ticks.theta)
delta.theta<-diff(theta.grid[1:2])
#ticks.theta<-length(theta.grid) #number of tick marks for theta

#delta.sigma2<-.02 #distance between points on sigma2 grid
ticks.sigma2<-200 #number of tick marks for sigma2
sigma2.grid<-seq(1,70,length.out=ticks.sigma2)
delta.sigma2<-diff(sigma2.grid[1:2])
#ticks.sigma2<-length(sigma2.grid) #number of tick marks for sigma2


#expand theta.grid and sigma.grid so that there are two vectors. Every possible pair of (theta,sigma2) is a single Row
#use expand.grid
expandGrid<-expand.grid(theta.grid,sigma2.grid)
#first column is the expanded theta vector, second column is the expanded sigma2 column



#unormalized posterior estimate function: estimates the posterior (unormalized) on the log scale
posterior.estimate<-function(theta,sigma2){
  #unormalized kernel density estimate of the postorior [theta, sigma2|t.obs,s.obs] on the log scale
  kDensityEst<-sum(log(dt((y-theta)/sqrt(sigma2), df=nu)/sqrt(sigma2)))+log(dnorm(theta,mean=eta,sd=tau))+log(dinvgamma(sigma2,alpha,beta))
  return(kDensityEst)
}


log.post.density<-mapply(FUN=posterior.estimate,expandGrid[,1],expandGrid[,2])
normalizing.constant<-sum(exp(log.post.density))*delta.theta*delta.sigma2
log.post.density.final<-cbind(expandGrid[,1],expandGrid[,2],exp(log.post.density)/normalizing.constant)


########################################
#Marginalizing
########################################


#Marginal for theta
post.theta.pdf<-sapply(theta.grid,FUN=function(x){delta.sigma2*sum(log.post.density.final[expandGrid[,1]==x,3])})
#Marginal for sigma2
post.sigma2.pdf<-sapply(sigma2.grid,FUN=function(x){delta.theta*sum(log.post.density.final[expandGrid[,2]==x,3])})

# dev.new()
#compare to t mcmc
plot(density(mcmcTSamps[,1]))
lines(theta.grid,post.theta.pdf, col=2)

plot(density(mcmcTSamps[,2]))
lines(sigma2.grid,post.sigma2.pdf,col=2)

save.image(do.call(paste, list('workSpaceFittingT5dfParamsetOnGrid', paramSet, ".RData", sep='')))



# 
# 
