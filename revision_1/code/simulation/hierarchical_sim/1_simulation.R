#
# Simulation for paper ----
#

# 'Complete' represents full data posteriors
# 'Incomplete' represents our restricted versions
library(parallel)



date() # SteveMod

# rm(list=ls())
# source('//trad/profs$/maceachern.1/Desktop/LEWIS.LEE/codeForPaper/codeForPaper/simulation/StevesimulationSourceCode.R')
source('SimulationSourceCode.R')
st <- Sys.time()

a<-0; b<-0



sigma2True<- 4 #true values of sigma2 to try for #the generation of data
#sigma2True <- 4  # SteveMod
#ass <- 1.25 # SteveMod
ass.vec <- c(1.25, 5 , 10) #values of shape parameter (alpha) on sigma_i^2 to use. #c(1.25, 2.5, 5 , 10, 20)
scale_vec <- c(0.5,1, 2) #c(0.5, 1/sqrt(2),1, sqrt(2), 2) #values of scale to use to define the scale parameter on sigma_i^2 where  beta_ss = alpha_ss * sig2True * scale. A value of 1 shrinks to the correct value of sigma2.  

#size of MCMC chain
nburn <- 1500
nkeep <- 1500

sims <- 21:30
set.seed(min(sims)) # for reproducible results

#print progess after printEach interations
printEach <- floor(nkeep/2)
#max number of iterations for the rlm fit
maxit <- 400

results_dir <- file.path(getwd(), 'results')
dir.create(results_dir, showWarnings = FALSE)

for(sig2True in sigma2True){
#sigma2True <- 4  # SteveMod
#sig2True <- 4    # SteveMod

dir_to_save_data <- file.path(getwd(), paste0('data_', "sig2_", sig2True))
dir.create(dir_to_save_data, showWarnings = FALSE)

for(data_sim in sims){
# Generate the data ----  
dat_file <- file.path(dir_to_save_data,  paste0('data_',data_sim,'.rds')) 
if(file.exists(dat_file)){ #checks to see if data file already exists - so that can run future sims with different priors without overwriting existing data and messing up MSE/KL calculations for already run sims. 
YList <- readRDS(dat_file)
} else {
YList <- fn.gen.data(sig2True)
saveRDS(YList, dat_file)
}
  factorsList<-YList$factorsList
  factorsMat<-matrix(unlist(factorsList), nrow=90, ncol=3, byrow = TRUE)
  p<-factorsMat[,1]
  n<-factorsMat[,2]
  m<-factorsMat[,3]
  thetaTrue<-YList$theta
  YList<-YList$yList

for(scale_a_ss in scale_vec){  
# cl <- makeCluster(detectCores(), type="FORK")
# registerDoParallel(cl) 
for(ass in ass.vec){    # SteveMod
# SteveMod 

# dir_to_save_in <- file.path(getwd(), paste0("a_s_", ass,'__', "scale_a_s_", round(scale_a_ss,2), '__', "sig2_", sig2True))
# dir.create(dir_to_save_in, showWarnings = FALSE)
  
bss <- ass * sig2True * scale_a_ss
#
# End SteveMod

#
# Fit the complete normal theory model---
#

#intitial values

thetaInt <- rnorm(90, sapply(YList, mean))

muInt <- rnorm(1, mean(thetaInt))

print(paste('Begin normal theory fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))

completeFit <- fn.complete.MCMC(YList, 
                              nkeep=nkeep, 
                              nburn=nburn, 
                              theta=thetaInt, 
                              sigma2=NULL, 
                              mu=muInt, 
                              tau2=NULL,
                              printEach=printEach, 
                              a=a, 
                              b=b, ass=ass, bss=bss) # SteveMod


# print(paste('End normal theory fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))

saveRDS(completeFit, file = file.path(results_dir, paste0('normal_', data_sim, '__as_', ass,'__', "scale_as_", round(scale_a_ss,2), '__', "sig2_", sig2True, '.rds')))

#
# Fit Incomplete method with Huber/Huber----
# Also fits individual robust regs to each group
#

thetaInt<-rnorm(90, sapply(YList, mean))
muInt<-rnorm(1, mean(thetaInt))

print(paste('Begin restricted Huber fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))

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
                                   b=b, ass=ass, bss=bss) # SteveMod

# print(paste('End restricted Huber fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))


saveRDS(incompleteFitHuber, file = file.path(results_dir, paste0('huber_', data_sim, '__as_', ass,'__', "scale_as_", round(scale_a_ss,2), '__', "sig2_", sig2True, '.rds')))


#
# Fit Incomplete method with Tukey/Huber----
# Also fits individual robust regs to each group
#


thetaInt<-rnorm(90, sapply(YList, mean))
muInt<-rnorm(1, mean(thetaInt))

print(paste('Begin restricted Tukey fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))



incompleteFitTukey <- fn.Incomplete.MCMC(YList,
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
                                       b=b, ass=ass, bss=bss) # SteveMod

# print(paste('End restricted Tukey fit ', data_sim, ': a_s set to', ass, ' scale_a_ss set to', scale_a_ss, ' and sigma2 set to', sig2True))

saveRDS(incompleteFitTukey, file = file.path(results_dir, paste0('tukey_', data_sim, '__as_', ass,'__', "scale_as_", round(scale_a_ss,2), '__', "sig2_", sig2True, '.rds')))

}
# stopCluster(cl)
  }
}
}



# # Save all needed results as a workspace ----
# #if(concTau2Near1){
# #  save.image('simulationWorkspace_concTau2Near1.RData')
# #} else {
# #save.image('simulationWorkspace.RData')
# #}
# # save.image('StevesimWork.RData')
# save.image('//trad/profs$/maceachern.1/Desktop/LEWIS.LEE/codeForPaper/codeForPaper/StevesimWork.RData') # SteveMod

Sys.time()-st

date() # SteveMod

