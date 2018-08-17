###########################
#Normal Hierarchical Model Analysis
###########################
rm(list=ls())
#setwd to code folder
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
tModelGroupBetaMeans<-array(NA, c(p,nGroups, reps))
tModelGroupBetaSDs<-array(NA, c(p,nGroups, reps))
tModelbetalConverge<-array(NA, c(p,nGroups, reps))

tModelBetaMeans<-array(NA, c(p, reps))
tModelBetaSDs<-array(NA, c(p, reps))
tModelBETAconverge<-array(NA, c(p, reps))

tModelSigma2Means<-array(NA, c(nGroups, reps))
tModelSigma2SDs<-array(NA, c(nGroups, reps))
tModelSigma2Converge<-array(NA, c(nGroups, reps))

tModelbstarConverge<-numeric(reps)
tModelMuRhoConverge<-numeric(reps)
tModelPsiRhoConverge<-numeric(reps)
tModelRhoConverge<-numeric(reps)

#predictions: This is a list of lists: For each rep there is a list of the predictions for each state (list length ==22) split by state
tModelPredListofLists<-replicate(reps, list()) 

#marginals
tModelMarginalsListofLists<-replicate(reps, list())


# 
# olsPredListofLists<-replicate(reps, list()) 
# 
# #marginals
# olsMarginalsListofLists<-replicate(reps, list())


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
      #fit<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdata)
      fit<-rlm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,scale.est='Huber',data=subdata, maxit=1000)
      y[[group]]<-subdata$sqrt_Count2012
      X[[group]]<-model.matrix(fit)
      #sigHats[group]<-summary(fit)$s
      sigHats[group]<-fit$s
#       step_Z[group]<-fn.compute.Z(summary(fit)$s^2, a0Star,b0Star) #abs(fn.compute.Z(summary(fit)$s^2))/5
#       if(is.na(step_Z[group])){step_Z[group]<-1}
      betaHats[,group]<-coef(fit)
      subdataHold<-holdoutSet[holdoutSet$Primary_Agency_State==names(yhold)[group],]
      if(nrow(subdataHold)>0){
        #only to get the holdout design matrix
        #fitHold<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdataHold)
        fitHold<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdataHold)
        yhold[[group]]<-subdataHold$sqrt_Count2012
        Xhold[[group]]<-model.matrix(fitHold)
        
        #olsPreds<-Xhold[[group]]%*%coef(fit)
        #olsPredListofLists[[i]][[group]]<-olsPreds
        #olsMarginalsListofLists[[i]][[group]]<-dnorm(yhold[[group]],olsPreds,summary(fit)$sigma)
        
      }
      
    }
    
#     
#     #tunning parameters
#     step_logbstar<-.1
#     mu_rho_step<-.1
#     psi_rho_step<-.3*(a_psir/b_psir)
#     rho_step<-.05
#     step_Z<-abs(step_Z)/5
#     #step_Z[7]<-.4
#     # step_Z[12]<-.4
#     
    
    #tunning parameters
    #step_logbstar<-abs(log(mu_bstr)/4)
    step_logbstar<-abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
    mu_rho_step<-.3 #(w1/(w1+w1))/sqrt(w1*w2/((w1+w2)^2*(w1+w1+1)))
    psi_rho_step<-a_psir^.5 #mean/sd
    rho_step<-.1
    #step_Z<-abs(step_Z)/5
    nis<-unlist(lapply(y, length), use.names=FALSE)
    step_Z<-fn.compute.Z(mean(sigHats^2), a0Star, b0Star)/(sqrt(nis))
    
    
   #T model Fit
    tModel<-heirTTheoryLm(y,
                              X,
                              nkeep,nburn,
                              mu0Star,
                              Sigma0Star,
                              a0Star, 
                          ((nu-2)/nu)*b0Star,
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
                              b_psir,
                          nu)
    
    
    #post means beta l
    betalSamps<-aperm(tModel$betal, c(1,3,2))
    #betal Converge?
    for(grp in 1:nGroups){
      tModelbetalConverge[,grp,i]<-abs(geweke.diag(mcmc(t(betalSamps[,,grp])))$z)
    }
    
    postMeansBetal<-apply(betalSamps,c(1,3) , mean)
    tModelGroupBetaMeans[,,i]<-postMeansBetal
    postSDsBetal<-apply(betalSamps,c(1,3) , sd)
    
    tModelGroupBetaSDs[,,i]<-postSDsBetal
    
    #post means Beta
    postMeansBETA<-colMeans(tModel$Beta)
    tModelBetaMeans[,i]<-postMeansBETA
    postSDsBeta<-apply(tModel$Beta,2 , sd)
    tModelBetaSDs[,i]<-postSDsBeta
    
    #Beta converge?
    tModelBETAconverge[,i]<-abs(geweke.diag(mcmc(tModel$Beta))$z)
    
    
    
    #post means sigma2s
    postMeansSigma2s<-colMeans(tModel$sigma2s)
    tModelSigma2Means[,i]<-postMeansSigma2s
    
    #post sds sigma2s
    postSDsSigma2s<-apply(tModel$sigma2s,2,sd)
    tModelSigma2SDs[,i]<-postSDsSigma2s
    #sigma2s converge?
    tModelSigma2Converge[,i]<-abs(geweke.diag(mcmc(tModel$sigma2s))$z)
    
    #bstar converge
    tModelbstarConverge[i]<-abs(geweke.diag(mcmc(tModel$bstar))$z)
    #mu_rho converge
    tModelMuRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$mu_rho))$z)
    #psi_rho_converge 
    tModelPsiRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$psi_rho))$z)
    #rho_converge
    tModelRhoConverge[i]<-abs(geweke.diag(mcmc(tModel$rho))$z)
    
    #-------------------------
    #preds on holdout set
    #-------------------------
    postMeansBetalList<-split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
    tModelPreds<-mapply(fits, postMeansBetalList,Xhold)
    tModelPredListofLists[[i]]<-tModelPreds
    
    
    #------------------------
    #computing marginal likelihoods for each element in holdout sample
    #------------------------
    tModelMarginalsListofLists[[i]]<-fn.compute.marginals.hierModelTmodel(betalSamps, tModel$sigma2s, yhold,Xhold)
    print(i)
  }
)



#------




assign(paste0("tModelGroupBetaMeans",n, "_run", run),tModelGroupBetaMeans)
assign(paste0("tModelGroupBetaSDs",n, "_run", run),    tModelGroupBetaSDs)
assign(paste0("tModelbetalConverge",n, "_run", run),    tModelbetalConverge)
assign(paste0("tModelBetaMeans",n, "_run", run),    tModelBetaMeans)
assign(paste0("tModelBetaSDs",n, "_run", run),    tModelBetaSDs)
assign(paste0("tModelBETAconverge",n, "_run", run),tModelBETAconverge)
assign(paste0("tModelSigma2Means",n, "_run", run),tModelSigma2Means)
assign(paste0("tModelSigma2SDs",n, "_run", run),tModelSigma2SDs)
assign(paste0("tModelSigma2Converge",n, "_run", run),tModelSigma2Converge)
assign(paste0("tModelbstarConverge",n, "_run", run),tModelbstarConverge)
assign(paste0("tModelMuRhoConverge",n, "_run", run),tModelMuRhoConverge)
assign(paste0("tModelPsiRhoConverge",n, "_run", run),tModelPsiRhoConverge)
assign(paste0("tModelRhoConverge",n, "_run", run),tModelRhoConverge)
assign(paste0("tModelPredListofLists",n, "_run", run),tModelPredListofLists)
assign(paste0("tModelMarginalsListofLists",n, "_run", run),tModelMarginalsListofLists)
assign(paste0("tModelStateCountsList",n,  "_run", run),stateCountsList)

assign(paste0("yholdoutMat",n,  "_run", run),get(paste0("yholdoutMat",n)))
assign(paste0("yOpenMat",n,  "_run", run),get(paste0("yOpenMat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))
assign(paste0("holdIndicesMatrix",n,  "_run", run),get(paste0("holdIndicesMatrix",n)))


rm(list=setdiff(ls(), c(ls()[grep('run',ls())], 'path.to.code','path.to.nonHierWorkSpaces','n','nu')))



#--------------------------
#save the workspace
#--------------------------
save.image(paste0('../workSpaces/wsNwTModelHierarchAnalysis_n',n,'run', run, '.RData'))

