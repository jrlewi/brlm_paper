#################
#John Lewis
#Robust inference via the Blended Paradigm
# Speed of Light Data Analysis
#
#Using Direct Sampling to evaluate the posterior distribution
# [theta,sigma^2|T(X),S(X)]
#  Focus Here is T(X) as Huber's location estimator and S(X) as Huber's proposal 2
#############################################


#################
#Libraries
#################
library(MASS)
library(ks)
library(MCMCpack)
################

data(newcomb)
n<-length(newcomb)
paramSet<-2  #set this to 1 or 2: 1 is the original set of hyper-parameters I came up with, 2 is from MacEachern's paper
source(do.call(paste, list('../',"newcombFixedParameters", paramSet, '.R', sep='')))
#source('../newcombFixedParameters1.R')

eta
tau
alpha
beta


#X_1,...,X_n iid N(theta, sigma^2)
#let T(X), S(X) be Huber estimators of location and scale

#Prior Distributions for theta and sigma^2 are specified as follows:
#theta~N(eta, tau^2)
#sigma^2~IG(alpha, beta)
#Analysis is based on the posterior [theta, sigma^2|T(X),S(X)]



#View the data 

# dev.new()
# pdf('newcomb.pdf')
# plot(density(newcomb), main='Observed Speed of Light Data', xlab='')
# dev.off()


#Specify the prior Parameters
#prior parameters #prior parameters for the mean are chosen by looking a
#estimates of the speed of light from before the experiment.
#these were determined in speedOflightAnalysisClean2 in the orderstatistics folder
# ##############################################################
# eta<--222.32 #value of the mean from previous estimates
# tau<-540.4248  #value of the standard deviations from previous estimates
# #the prior for sigma^2 is  an inverse gamma
# #inverse gamma parameters (Note R's parameterization)
# alpha<-5   #shape parameter
# beta<-50*(alpha-1)   #this gives a mean of 50
#####################################################
#beta^2/((alpha-1)^2*(alpha-2))

start<-Sys.time()



####################
#Summarize the Data with the observed value of the robust summary statistics
#####################
rlm.obs<-rlm(newcomb~1,scale.est='Huber')
t.obs<-rlm.obs$coefficients  #observed location estimate
s.obs<-rlm.obs$s             #observed scale estimate
log.s.obs<-log(s.obs)
##############################################################



#################
#Step 1: 
#Generation of N samples from the 'standard' distribution (in this case the standard normal)
#Summarize each by (T,S)_i i in 1,...,N. (These are the 'robust' summary statistics in this case)
##################
set.seed(1)
N<-1e5 #number of samples of size n to draw and summarize with T and S


########################
#define a function to output Huber's estimators 
#scale estimates on the log scale: this isn't necessary because the distribution of the estimates is fairly symmetric
#########################
huber.estimators<-function(X){
  rlm.out<-rlm(X~1, scale.est='Huber')
  return(c(as.numeric(rlm.out$coefficients), log(as.numeric(rlm.out$s))))
}
#############
#Use this function to calculate huber's statistics for each of the N samples
#############
X<-matrix(rnorm(N*n,0,1),nrow=N,ncol=n) #each row is a sample from N(0,1)
huber.matrix<-apply(X,MARGIN=1, FUN=huber.estimators)
###############################################################################

##############################################
#Step 2: Bandwidth Selection
##############################################

#For Now Just use the one dimensional plug in estimator computed using hpi function in the ks package (Wand and Jones 1994)


h1<-hpi(huber.matrix[1,], binned=TRUE)
h2<-hpi(huber.matrix[2,], binned=TRUE)
H<-diag(c(h1,h2))
H #represents the 'sd matrix' for tha Guassian Kernel


##############################################
#Step 3: Density Estimate for theta and sigma2
##############################################
#


#define the grid
#delta.theta<-.02  #distance between points on theta grid 
ticks.theta<-500
theta.grid<-seq(20,35,length.out=ticks.theta)
delta.theta<-diff(theta.grid[1:2])
#ticks.theta<-length(theta.grid) #number of tick marks for theta

#delta.sigma2<-.02 #distance between points on sigma2 grid
ticks.sigma2<-500 #number of tick marks for sigma2
sigma2.grid<-seq(5,70,length.out=ticks.sigma2)
delta.sigma2<-diff(sigma2.grid[1:2])
#ticks.sigma2<-length(sigma2.grid) #number of tick marks for sigma2


#expand theta.grid and sigma.grid so that there are two vectors. Every possible pair of (theta,sigma2) is a single Row
#use expand.grid
expandGrid<-expand.grid(theta.grid,sigma2.grid)
#first column is the expanded theta vector, second column is the expanded sigma2 column


#kernel density function using H as the bandwidth and huber.matrix as the data
kernel_density<-function(x1,x2){
  mean(dmvnorm(t(huber.matrix),mean=c(x1,x2),sigma=H^2))
}


#unormalized posterior estimate function: estimates the posterior (unormalized) on the log scale
posterior.estimate<-function(theta,sigma2){
  #unormalized kernel density estimate of the postorior [theta, sigma2|t.obs,s.obs] on the log scale
  kDensityEst<-log(kernel_density((t.obs-theta)/sqrt(sigma2),log.s.obs-.5*log(sigma2))/sqrt(sigma2))+log(dnorm(theta,mean=eta,sd=tau))+log(dinvgamma(sigma2,alpha,beta))
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
# pdf('marginaldistributionsRestricted.pdf')
# plot(theta.grid,post.theta.pdf, type='l')
# plot(sigma2.grid,post.sigma2.pdf, type='l')
# dev.off()

# pdfThetaRestricted<-cbind(theta.grid, post.theta.pdf)
# pdfSigma2Restricted<-cbind(sigma2.grid, post.sigma2.pdf)

#write the values to text files
# write.table(cbind(theta.grid,post.theta.pdf),file='postThetapdf.txt', row.names=FALSE,col.names=FALSE)
# write.table(cbind(sigma2.grid,post.sigma2.pdf),file='postSigma2pdf.txt', row.names=FALSE,col.names=FALSE)
# write.table(log.post.density.final,file='jointPosterior.txt', row.names=FALSE,col.names=FALSE)
save.image(do.call(paste, list('wsNewcombDirectSamplingHuber',paramSet, '.RData', sep='')))

end<-Sys.time()
end-start
##############
##END CODE####
##############



