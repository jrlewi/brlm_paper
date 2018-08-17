############################
#Source Code for the Hierarchical models
#########################################

n<-1000 #sample size for the training data
run<-4

reps<-25 #number of repetitions

nburn<-2e4 #set length of mcmc chains
nkeep<-5e4

nu<-3 #df for t-model
#setwd to code folder

path.to.nonHierWorkSpaces<-"../../nwLinearReg/workSpaces/"

#source results from nonhierarchical fits to get the training/holdout sets

#if statement for n=2000 
#load the y's from the non-Hierarchical Fits
if(n>=1000){
load(paste0(path.to.nonHierWorkSpaces, 'wsNwAnalysis_n',n,'run',run, '.RData'))
} else {
  load(paste0(path.to.nonHierWorkSpaces, 'wsNwAnalysis_n',n,'.RData'))
}

#only keep info on y's
rm(list=setdiff(ls(), list(paste0('holdIndicesMatrix', n),paste0('yholdoutMat',n),paste0('yOpenMat',n),paste0('yOpenType1Mat',n),"path.to.nonHierWorkSpaces",'n',"run", "reps", "nburn", "nkeep", "nkeept",'nu')))

maxit<-5000 #for max iterations of the rlm fitting function

#source the fitting functions 

source('nwHierNormalModelFittingFunction.R')
source('nwHierRestrictedModelFittingFunction.R')
source('nwTdistHierModelFittingFunction.R')

#load prior info; this also loads the analysis set called analysisSetCs
load(paste0(path.to.nonHierWorkSpaces,"nwdataPaper1PriorConstructionWorkSpace.RData"))
load("../workSpaces/wsnwdataHierPriorConstruction.RData")

nGroups<-length(unique(analysisSetCs$Primary_Agency_State))
N<-nrow(analysisSetCs)
rownames(analysisSetCs)<-1:N #this is so holdindices/trainIndices match with rownames

#function to compute predicted values and/or fits
fits<-function(betahats, X){X%*%betahats}#; names(fits)<-rownames(X); fits}

#function to compute marginals for restricted and full models
fn.compute.marginals.hierModelNormal<-function(betalsamples, sigma2lsamples, yhold,Xhold){
  #betalsamples the array of betals: in the specific format: the 3 dimension is the groups. columns represent samples, row represent slopes
  #sigma2lsampls: #columns represnt groups, rows represent samples
  
  betalSampsList<-lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList)<-names(yhold)
  #holdout means
  nMuList<-mapply(fits, betalSampsList, Xhold)
  
  #sd's across samples for each houldout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)
  
  mapply(FUN=function(yh,mean,sd){
    rowMeans(dnorm(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE)))
  }
         ,yhold,nMuList,nSigmaList)
  
}


#------------
#function to compute marginals for The T model; 
#-----------

tdensity<-function(y, mean, sigma){
  (gamma(.5*(nu+1))/(gamma(.5*nu)*sigma*sqrt(nu*pi)))*(1+((y-mean)/sigma)^2/nu)^(-.5*(nu+1))
}



fn.compute.marginals.hierModelTmodel<-function(betalsamples, sigma2lsamples, yhold,Xhold){
  #betalsamples the array of betals: in the specific format: the 3rd dimension is the groups. columns represent samples, row represent slopes
  #sigma2lsampls: #columns represnt groups, rows represent samples
  
  betalSampsList<-lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList)<-names(yhold)
  #holdout means
  nMuList<-mapply(fits, betalSampsList, Xhold)
  
  #sd's across samples for each houldout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)
  
  mapply(FUN=function(yh,mean,sd){
    rowMeans(tdensity(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE)))
  }
         ,yhold,nMuList,nSigmaList)
  
}


