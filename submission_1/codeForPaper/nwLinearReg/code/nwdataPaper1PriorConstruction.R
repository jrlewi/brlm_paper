#########################################
#John Lewis
#Prior Construction: Using the data from 2008-2010
#data was prepared on NW laptop: 2008-2010 for prior constuction and 2010-2012 for data
#########################################

load("../workSpaces/nwdataPrepPaper1DataCleanWorkspace.RData") #not included
library(plyr)
library(MASS)
library(MCMCpack)

#######################################
#For paper; remove states with fewer than 25 agencies
stCounts<-table(analysisSet$Primary_Agency_State)



pairs(~sqrt_Count2010+sqrt_Count2008+Associate_Count+ageMeasure, data=priorDataSet, pch=19, cex=.3)


# rlmfit<-rlm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure, psi=psi.bisquare, scale.est='Huber',data=analysisSet, maxit=1000)
# rlmfit$s^2
# curve(dinvgamma(x,(4+priorRlm$s^2)/priorRlm$s^2,4), from=0, to=10)
# 
#regression to for prior information

priorRlm<-rlm(sqrt_Count2010~sqrt_Count2008+Associate_Count+ageMeasure, psi=psi.bisquare, scale.est='Huber',data=priorDataSet, maxit=400)
priorRlm$s^2
priorOls<-lm(sqrt_Count2010~sqrt_Count2008+Associate_Count+ageMeasure, data=priorDataSet)
summary(priorOls)$sigma^2

head(priorDataSet)
nPrior<-nrow(priorDataSet)
#plot(priorRlm$residuals)
summary(priorRlm)
###############
#prior mean and variance for beta
###############
mu0<-as.numeric(coef(priorRlm))
#construction of covariance matrix estimate; refer to page 108 maronna
X<-model.matrix(priorRlm)
resids<-priorRlm$residuals
# resids2<-priorDataSet$sqrt_Count2010-fitted(priorRlm)
# plot(resids,resids)
n<-nrow(X)
p<-ncol(X)

sigmahat<-priorRlm$s
meanPsiSq<-mean(((resids/sigmahat)*psi.bisquare(resids/sigmahat))^2)
meanPsiPrime<-mean(psi.bisquare(resids/sigmahat, deriv=1))
vhat<-(n/(n-p))*sigmahat^2*meanPsiSq/(meanPsiPrime^2)
Sigma0<-vhat*solve(t(X)%*%X)
vcov(priorRlm) #about the same
Sigma0<-n*Sigma0
Sigma0

######################
#Prior parameters for sigma2
######################
#need an inverse gamma with sigmahat^2 as its mean
b0<-4
a0<-(b0+sigmahat^2)/sigmahat^2

b0/(a0-1)
sigmahat^2 #var(analysisSet$sqrt_Count2012)
b0^2/((a0-1)^2*(a0-2))

#plot the prior for sigma2
curve(dinvgamma(x,a0,b0), xlim=c(0,20))
# mean(rinvgamma(100000,a0,b0))

#########################
#prior value for df parameter for the t-model
#########################
nu<-3
###############################
#################################
#Centering and scaling the variables in the analysis set and adjusting the priors for beta and sigma2 accordingly
#################################

View(analysisSet)

#plot(analysisSet$Ending_HH_Count2010,analysisSet$Ending_HH_Count2012)
#center and scaling the analysis set
lapply(analysisSet,class)
centScale<-scale(analysisSet[,4:10])
analysisSetCs<-cbind(analysisSet[,1:3],centScale)
View(analysisSetCs)
plot(analysisSet$sqrt_Count2010,analysisSet$sqrt_Count2012)
plot(analysisSetCs$sqrt_Count2010,analysisSetCs$sqrt_Count2012)

#means of x in order  sqrt_Count2010 Associate_Count ageMeasure
xbar<-apply(analysisSet[,c("sqrt_Count2010", "Associate_Count", "ageMeasure")],2, mean)
#sd of x's
sx<-apply(analysisSet[,c("sqrt_Count2010", "Associate_Count", "ageMeasure")],2, sd)

ybar<-mean(analysisSet[,"sqrt_Count2012"])
sy<-sd(analysisSet[,"sqrt_Count2012"])

#prior adjustments for beta based on the center and scaled scale
#beta*=A%*%beta-c
#with beta=(beta_0, BETA1) where beta_0 is a scale for the intercept and BETA1 represent the reg. coefficients we have
#(I-(1/n)J)Y/sy=(I-(1/n)J)XS, with S=diag(sx2,....sxp): i.e. the sd's for each covariate
#then solving for Y=...we get the relationship #beta*=A%*%beta-c with A and c defined below:
#that is the we have the trasnformation to the centered and scaled scale

#beta~normal(mu_beta,Sigma_beta) implies the distribution for beta*=A%*%beta-c~normal(A%*%mu_beta-c, A%*%Sigma_beta%*%t(A))

A<-rbind(c(1, xbar), cbind(rep(0,p-1), diag(sx)))/sy
c<-c(ybar,rep(0,p-1))/sy
mu0Star<-A%*%mu0-c
Sigma0Star<-A%*%Sigma0%*%t(A)
round(Sigma0Star,4)
#prior adjustment to sigma2 #just change of variables sigma^2~IG(a0,b0) implies sigma^2/sy^2~IG(a0,b0/sy^2)
#sigma^2*=sigma^2/sy^2
a0Star<-a0
b0Star<-b0/sy^2

b0Star/(a0Star-1)
sigmahat^2/sy^2
b0Star^2/((a0Star-1)^2*(a0Star-2))
sqrt(b0Star^2/((a0Star-1)^2*(a0Star-2)))
curve(dinvgamma(x, a0Star,b0Star))

################################
#mask analysisSetCs
###############################

str(analysisSetCs)
class(analysisSetCs$Agency_Type)
levels(analysisSetCs$Agency_Type)
#Agency_Type to Type #
nTypes<-length(levels(analysisSetCs$Agency_Type))
orderTypes<-c(1:nTypes)[c(7, 2:6,1,8:nTypes)]
levels(analysisSetCs$Agency_Type)<-do.call(paste, list('Type',orderTypes, sep=' '))
levels(analysisSetCs$Agency_Type)
View(analysisSetCs)

#Agency_Contract_type to Con #
nCon<-length(levels(analysisSetCs$Agency_Contract_Type))
set.seed(123)
orderCon<-c(1:nCon)[sample(1:nCon)]
levels(analysisSetCs$Agency_Contract_Type)<-do.call(paste, list('Con',orderCon, sep=' '))
levels(analysisSetCs$Agency_Contract_Type)
View(analysisSetCs)



#Primary_Agency_State to State #
nState<-length(levels(analysisSetCs$Primary_Agency_State))
set.seed(123)
orderState<-c(1:nState)[sample(1:nState)]
levels(analysisSetCs$Primary_Agency_State)<-do.call(paste, list('State',orderState, sep=' '))
levels(analysisSetCs$Primary_Agency_State)
View(analysisSetCs)

levels(analysisSet$Primary_Agency_State)


rm(list=setdiff(ls(),c('mu0Star','Sigma0Star', 'a0Star','b0Star', 'nu', 'analysisSetCs' )))
save.image("../workSpaces/nwdataPaper1PriorConstructionWorkSpace.RData")
