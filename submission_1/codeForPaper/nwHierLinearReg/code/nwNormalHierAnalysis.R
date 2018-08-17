###########################
#Normal Hierarchical Model Analysis
###########################
rm(list=ls())
#setwd(code folder)
source("sourceCodeForAnalysis.R")
#n's from here (25,100,1000,2000)
n
run
reps
nburn
nkeep


Sys.info()


#run multiple runs
seeds<-c(3,4,5,6,7)
#set seed after this load it brings in a seed already
set.seed(seeds[run])


p<-length(mu0Star)
#for run in run....


#-------------------
#outputs to save
#--------------------
nTheoryGroupBetaMeans<-array(NA, c(p,nGroups, reps))
nTheoryGroupBetaSDs<-array(NA, c(p,nGroups, reps))
nTheorybetalConverge<-array(NA, c(p,nGroups, reps))

nTheoryBetaMeans<-array(NA, c(p, reps))
nTheoryBetaSDs<-array(NA, c(p, reps))
nTheoryBETAconverge<-array(NA, c(p, reps))

nTheorySigma2Means<-array(NA, c(nGroups, reps))
nTheorySigma2SDs<-array(NA, c(nGroups, reps))
nTheorySigma2Converge<-array(NA, c(nGroups, reps))

nTheorybstarConverge<-numeric(reps)
nTheoryMuRhoConverge<-numeric(reps)
nTheoryPsiRhoConverge<-numeric(reps)
nTheoryRhoConverge<-numeric(reps)

#predictions: This is a list of lists: For each rep there is a list of the predictions for each state (list length ==22) split by state
nTheoryPredListofLists<-replicate(reps, list()) 

#marginals
nTheoryMarginalsListofLists<-replicate(reps, list())



olsPredListofLists<-replicate(reps, list()) 

#marginals
olsMarginalsListofLists<-replicate(reps, list())


stateCountsList<-list();length(stateCountsList)<-reps

#------------------------------------




system.time(  
for(i in 1:reps){


hold<-get(paste0("holdIndicesMatrix",n))[i,]
holdoutSet<-analysisSetCs[hold,]

trainingSet<-analysisSetCs[-hold,]
stateCounts<-table(trainingSet$Primary_Agency_State)
#countsBigEnough<-length(stateCounts) #all states remainin have large enough counts
stateCountsList[[i]]<-stateCounts

#prepare data for fitting function; get ols preds and marginals on holdoutset
#y is list of responses from each group
#X is list of design matrices for each group
y<-list(); length(y)<-sum(length(stateCounts)); names(y)<-names(stateCounts)
X<-list(); length(X)<-sum(length(stateCounts)); names(X)<-names(stateCounts)
betaHats<-matrix(NA, p,length(y))
sigHats<-numeric(length(y))
step_Z<-numeric(length(y))


#prepare the holdout data for predictions
#yhold is list of holdout responses from each group
#Xhold is list of holfout design matrices for each group
yhold<-list(); length(yhold)<-sum(length(stateCounts)); names(yhold)<-names(stateCounts)
Xhold<-list(); length(Xhold)<-sum(length(stateCounts)); names(Xhold)<-names(stateCounts)



for(group in 1:length(y)){
  subdata<-trainingSet[trainingSet$Primary_Agency_State==names(y)[group],]
  fit<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdata)
  y[[group]]<-subdata$sqrt_Count2012
  X[[group]]<-model.matrix(fit)
  sigHats[group]<-summary(fit)$s
  #step_Z[group]<-fn.compute.Z(summary(fit)$s^2, a0Star,b0Star) #abs(fn.compute.Z(summary(fit)$s^2))/5
  #if(is.na(step_Z[group])){step_Z[group]<-1}
  betaHats[,group]<-coef(fit)
  subdataHold<-holdoutSet[holdoutSet$Primary_Agency_State==names(yhold)[group],]
  if(nrow(subdataHold)>0){
    #only to get the holdout design matrix
    fitHold<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdataHold)
    yhold[[group]]<-subdataHold$sqrt_Count2012
    Xhold[[group]]<-model.matrix(fitHold)
  
    olsPreds<-Xhold[[group]]%*%coef(fit)
    olsPredListofLists[[i]][[group]]<-olsPreds
    olsMarginalsListofLists[[i]][[group]]<-dnorm(yhold[[group]],olsPreds,summary(fit)$sigma)
  
  } else{
    print(paste(group,'has no agencies in holdout set for rep', i))
  }
  
}


#tunning parameters
#step_logbstar<-abs(log(mu_bstr)/4)
step_logbstar<-abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
mu_rho_step<-.3 #(w1/(w1+w1))/sqrt(w1*w2/((w1+w2)^2*(w1+w1+1)))
psi_rho_step<-a_psir^.5 #mean/sd
rho_step<-.1
#step_Z<-abs(step_Z)/5
nis<-unlist(lapply(y, length), use.names=FALSE)
step_Z<-fn.compute.Z(mean(sigHats^2), a0Star, b0Star)/(sqrt(nis))

#full fit
nTheory<-heirNormTheoryLm(y,
                                   X,
                                   nkeep,nburn,
                                   mu0Star,
                                   Sigma0Star,
                                   a0Star, 
                                   b0Star,
                                   #alpha_mustr,
                                   #beta_mustr,
                                   # a_psib,
                                   # b_psib,
                                   mu_bstr,
                                   psi_bstr,
                                   swSq=1,
                                   w1,
                                   w2, 
                                   a_psir,
                                   b_psir)


#post means beta l
betalSamps<-aperm(nTheory$betal, c(1,3,2))
#betal Converge?
for(grp in 1:nGroups){
  nTheorybetalConverge[,grp,i]<-abs(geweke.diag(mcmc(t(betalSamps[,,grp])))$z)
}

postMeansBetal<-apply(betalSamps,c(1,3) , mean)
nTheoryGroupBetaMeans[,,i]<-postMeansBetal
postSDsBetal<-apply(betalSamps,c(1,3) , sd)
nTheoryGroupBetaSDs[,,i]<-postSDsBetal

#post means Beta
postMeansBETA<-colMeans(nTheory$Beta)
nTheoryBetaMeans[,i]<-postMeansBETA
postSDsBeta<-apply(nTheory$Beta,2 , sd)
nTheoryBetaSDs[,i]<-postSDsBeta

#Beta converge?
nTheoryBETAconverge[,i]<-abs(geweke.diag(mcmc(nTheory$Beta))$z)



#post means sigma2s
postMeansSigma2s<-colMeans(nTheory$sigma2s)
nTheorySigma2Means[,i]<-postMeansSigma2s

#post sds sigma2s
postSDsSigma2s<-apply(nTheory$sigma2s,2,sd)
nTheorySigma2SDs[,i]<-postSDsSigma2s
#sigma2s converge?
nTheorySigma2Converge[,i]<-abs(geweke.diag(mcmc(nTheory$sigma2s))$z)

#bstar converge
nTheorybstarConverge[i]<-abs(geweke.diag(mcmc(nTheory$bstar))$z)
#mu_rho converge
nTheoryMuRhoConverge[i]<-abs(geweke.diag(mcmc(nTheory$mu_rho))$z)
#psi_rho_converge 
nTheoryPsiRhoConverge[i]<-abs(geweke.diag(mcmc(nTheory$psi_rho))$z)
#rho_converge
nTheoryRhoConverge[i]<-abs(geweke.diag(mcmc(nTheory$rho))$z)

#-------------------------
#preds on holdout set
#-------------------------
postMeansBetalList<-split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
nTheoryPreds<-mapply(fits, postMeansBetalList,Xhold)
nTheoryPredListofLists[[i]]<-nTheoryPreds


#------------------------
#computing marginal likelihoods for each element in holdout sample
#------------------------
nTheoryMarginalsListofLists[[i]]<-fn.compute.marginals.hierModelNormal(betalSamps, nTheory$sigma2s, yhold,Xhold)
print(i)
}
)



#------




assign(paste0("nTheoryGroupBetaMeans",n, "_run", run),nTheoryGroupBetaMeans)
assign(paste0("nTheoryGroupBetaSDs",n, "_run", run),    nTheoryGroupBetaSDs)
assign(paste0("nTheorybetalConverge",n, "_run", run),    nTheorybetalConverge)
assign(paste0("nTheoryBetaMeans",n, "_run", run),    nTheoryBetaMeans)
assign(paste0("nTheoryBetaSDs",n, "_run", run),    nTheoryBetaSDs)
assign(paste0("nTheoryBETAconverge",n, "_run", run),nTheoryBETAconverge)
assign(paste0("nTheorySigma2Means",n, "_run", run),nTheorySigma2Means)
assign(paste0("nTheorySigma2SDs",n, "_run", run),nTheorySigma2SDs)
assign(paste0("nTheorySigma2Converge",n, "_run", run),nTheorySigma2Converge)
assign(paste0("nTheorybstarConverge",n, "_run", run),nTheorybstarConverge)
assign(paste0("nTheoryMuRhoConverge",n, "_run", run),nTheoryMuRhoConverge)
assign(paste0("nTheoryPsiRhoConverge",n, "_run", run),nTheoryPsiRhoConverge)
assign(paste0("nTheoryRhoConverge",n, "_run", run),nTheoryRhoConverge)
assign(paste0("nTheoryPredListofLists",n, "_run", run),nTheoryPredListofLists)
assign(paste0("nTheoryMarginalsListofLists",n, "_run", run),nTheoryMarginalsListofLists)
assign(paste0("olsPredListofLists",n, "_run", run),olsPredListofLists)
assign(paste0("olsMarginalsListofLists",n, "_run", run),olsMarginalsListofLists)
assign(paste0("nTheoryStateCountsList",n,  "_run", run),stateCountsList)

assign(paste0("yholdoutMat",n,  "_run", run),get(paste0("yholdoutMat",n)))
assign(paste0("yOpenMat",n,  "_run", run),get(paste0("yOpenMat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))
assign(paste0("holdIndicesMatrix",n,  "_run", run),get(paste0("holdIndicesMatrix",n)))


rm(list=setdiff(ls(), c(ls()[grep('run',ls())], 'path.to.code','path.to.nonHierWorkSpaces','n')))



#--------------------------
#save the workspace
#--------------------------
save.image(paste0('../workSpaces/wsNwNormalHierarchAnalysis_n',n,'run', run, '.RData'))


