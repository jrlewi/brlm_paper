#
# Simulation for paper ----
#

# 'Complete' represents full data posteriors
# 'Incomplete' represents our restricted versions


rm(list=ls())
source('simulationSourceCode.R')
st<-Sys.time()
set.seed(123) # for reproducible results

concTau2Near1<-0 #1: then concentrate tau^2 prior near 1. 0: then don't and a=b=0

if(concTau2Near1){
  mu<-1
  v<-.00001
  a<-mu^2/v-2 #most prob between .99 and 1.01
  b<-mu*(a-1)
} else {
  a<-0; b<-0
}


sigma2True<-c(.5,1,4,10) #true values of sigma2 to try for the generation of data

#size of MCMC chain
nburn<-1e4
nkeep<-1e4
  
#print progess after printEach interations
printEach<-floor(nkeep/10)
#max number of iterations for the rlm fit
maxit<-400

for (sig2True in sigma2True){

#
# Generate the data ----
#

YList<-fn.gen.data(sig2True)
factorsList<-YList$factorsList
factorsMat<-matrix(unlist(factorsList), nrow=90, ncol=3, byrow = TRUE)
p<-factorsMat[,1]
n<-factorsMat[,2]
m<-factorsMat[,3]
thetaTrue<-YList$theta

YList<-YList$yList

assign(paste0('YList', 'Sig2_', sig2True), YList)
assign(paste0('thetaTrue', 'Sig2_', sig2True), thetaTrue)

#
# Fit the complete normal theory model---
#


#intitial values

thetaInt<-rnorm(90, sapply(YList, mean))

muInt<-rnorm(1, mean(thetaInt))

print(paste('begin normal theory fit- sigma2', sig2True))

completeFit<-fn.complete.MCMC(YList, 
                              nkeep=nkeep, 
                              nburn=nburn, 
                              theta=thetaInt, 
                              sigma2=NULL, 
                              mu=muInt, 
                              tau2=NULL,
                              printEach=printEach, 
                              a=a, 
                              b=b)

print(paste('end normal theory fit- sigma2', sig2True))


assign(paste0('normalTheoryFit', 'Sig2_', sig2True), completeFit)
rm(completeFit)

#
# Fit Incomplete method with Huber/Huber----
# Also fits individual robust regs to each group
#

thetaInt<-rnorm(90, sapply(YList, mean))
muInt<-rnorm(1, mean(thetaInt))

print(paste('begin restricted Huber fit- sigma2', sig2True))

incompleteFitHuber<-fn.Incomplete.MCMC(YList,
                                   regEst='Huber',
                                   scaleEst='Huber',
                                   nkeep=nkeep, 
                                   nburn=nburn, 
                                   theta=thetaInt, 
                                   sigma2=NULL, 
                                   mu=muInt, 
                                   tau2=NULL, 
                                   printEach=printEach,
                                   maxit=maxit,
                                   a=a, 
                                   b=b)

print(paste('end restricted Huber fit- sigma2', sig2True))

assign(paste0('restrictedFitHuber', 'Sig2_', sig2True), incompleteFitHuber)
rm(incompleteFitHuber)


#
# Fit Incomplete method with Tukey/Huber----
# Also fits individual robust regs to each group
#


thetaInt<-rnorm(90, sapply(YList, mean))
muInt<-rnorm(1, mean(thetaInt))

print(paste('begin restricted Tukey fit- sigma2', sig2True))

incompleteFitTukey<-fn.Incomplete.MCMC(YList,
                                       regEst='Tukey',
                                       scaleEst='Huber',
                                       nkeep=nkeep, 
                                       nburn=nburn, 
                                       theta=thetaInt, 
                                       sigma2=NULL, 
                                       mu=muInt, 
                                       tau2=NULL, 
                                       printEach=printEach,
                                       maxit=maxit,
                                       a=a,
                                       b=b)

print(paste('end restricted Tukey fit- sigma2', sig2True))

assign(paste0('restrictedFitTukey', 'Sig2_', sig2True), incompleteFitTukey)
rm(incompleteFitTukey)
rm(YList)

} #end loop



# Save all needed results as a workspace ----
if(concTau2Near1){
  save.image('simulationWorkspace_concTau2Near1.RData')
} else {
save.image('simulationWorkspace.RData')
}

Sys.time()-st

