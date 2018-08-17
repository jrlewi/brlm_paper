###########################
#Restricted Hierarchical Model Analysis
###########################
rm(list=ls())
#setwd to code folder
source("sourceCodeForAnalysis.R")

#regEst<-'Tukey'
regEst<-'Huber'
regEst



############################
#define the psi function
############################
if(regEst=='Huber') {
  psi<-get('psi.huber') #internal
  fn.psi<-get('fn.psi.huber')
  
} else { 
  if(regEst=='Tukey'){
    psi<-get('psi.bisquare') #internal
    fn.psi<-get('fn.psi.bisquare')
  } else {stop("only set up for Huber or Tukey regression estimates")}}


#n's from here (25,100,1000,2000)
n
run
reps
nburn
nkeep


Sys.info()


#run multiple runs
seeds<-c(3,4,5,6)
#set seed after this load it brings in a seed already
set.seed(seeds[run])


p<-length(mu0Star)
#for run in run....


#-------------------
#outputs to save
#--------------------
restrictedGroupBetaMeans<-array(NA, c(p,nGroups, reps))
restrictedGroupBetaSDs<-array(NA, c(p,nGroups, reps))
restrictedbetalConverge<-array(NA, c(p,nGroups, reps))

restrictedBetaMeans<-array(NA, c(p, reps))
restrictedBetaSDs<-array(NA, c(p, reps))
restrictedBETAconverge<-array(NA, c(p, reps))

restrictedSigma2Means<-array(NA, c(nGroups, reps))
restrictedSigma2SDs<-array(NA, c(nGroups, reps))
restrictedSigma2Converge<-array(NA, c(nGroups, reps))

restrictedbstarConverge<-numeric(reps)
restrictedMuRhoConverge<-numeric(reps)
restrictedPsiRhoConverge<-numeric(reps)
restrictedRhoConverge<-numeric(reps)

#predictions: This is a list of lists: For each rep there is a list of the predictions for each state (list length ==22) split by state
restrictedPredListofLists<-replicate(reps, list()) 

#marginals
restrictedMarginalsListofLists<-replicate(reps, list())



rlmPredListofLists<-replicate(reps, list()) 

#marginals
rlmMarginalsListofLists<-replicate(reps, list())

#acceptance rates
#yAccept<-array(NA, c(nburn+nkeep, nGroups, reps)) #saving all 0/1s
yAccept<-array(NA, c(reps, nGroups)) #just the column means
#------------------------------------
stateCountsList<-list();length(stateCountsList)<-reps



system.time(  
  for(i in 1:reps){
    
    
    hold<-get(paste0("holdIndicesMatrix",n))[i,]
    holdoutSet<-analysisSetCs[hold,]
    
    trainingSet<-analysisSetCs[-hold,]
    stateCounts<-table(trainingSet$Primary_Agency_State)
    #print(stateCounts)
    #countsBigEnough<-length(stateCounts) #all states remainin have large enough counts
    stateCountsList[[i]]<-stateCounts
    
    #prepare data for fitting function; get rlm preds and marginals on holdoutset
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
      fit<-rlm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,psi=psi,scale.est='Huber',data=subdata, maxit=1000)
      y[[group]]<-subdata$sqrt_Count2012
      X[[group]]<-model.matrix(fit)
      sigHats[group]<-fit$s
#       step_Z[group]<-fn.compute.Z(fit$s^2, a0Star,b0Star) 
#       if(is.na(step_Z[group])){step_Z[group]<-1}
      betaHats[,group]<-coef(fit)
      subdataHold<-holdoutSet[holdoutSet$Primary_Agency_State==names(yhold)[group],]
      if(nrow(subdataHold)>0){
        #only to get the holdout design matrix
        fitHold<-lm(sqrt_Count2012~sqrt_Count2010+Associate_Count+ageMeasure,data=subdataHold)
        yhold[[group]]<-subdataHold$sqrt_Count2012
        Xhold[[group]]<-model.matrix(fitHold)
        
        rlmPreds<-Xhold[[group]]%*%coef(fit)
        rlmPredListofLists[[i]][[group]]<-rlmPreds
        rlmMarginalsListofLists[[i]][[group]]<-dnorm(yhold[[group]],rlmPreds,fit$s)
      
      } else {
        print(paste(group,'has no agencies in holdout set for rep', i))
      }
      
    }
    
    

    
    
#     #tunning parameters
#     step_logbstar<-.1
#     mu_rho_step<-.1
#     psi_rho_step<-.3*(a_psir/b_psir)
#     rho_step<-.05
#     step_Z<-abs(step_Z)/3
#     #step_Z[7]<-.4
#     # step_Z[12]<-.4
    
    
    #tunning parameters
    #step_logbstar<-abs(log(mu_bstr)/4)
    step_logbstar<-abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
    mu_rho_step<-.3 #(w1/(w1+w1))/sqrt(w1*w2/((w1+w2)^2*(w1+w1+1)))
    psi_rho_step<-a_psir^.5 #mean/sd
    rho_step<-.1
    #step_Z<-abs(step_Z)/5
    nis<-unlist(lapply(y, length), use.names=FALSE)
    step_Z<-fn.compute.Z(mean(sigHats^2), a0Star, b0Star)/(sqrt(nis))
    
    
    
    
    
    #restricted fit
    restricted<-hierNormTheoryRestLm2(y,
                                                X,
                                                regEst,
                                                scaleEst='Huber',
                                                nkeep, nburn,
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
                                                b_psir,
                                                maxit=maxit)
    
    
    #post means beta l
    betalSamps<-aperm(restricted$betal, c(1,3,2))
    #betal Converge?
    for(grp in 1:nGroups){
      restrictedbetalConverge[,grp,i]<-abs(geweke.diag(mcmc(t(betalSamps[,,grp])))$z)
    }
    
    postMeansBetal<-apply(betalSamps,c(1,3) , mean)
    restrictedGroupBetaMeans[,,i]<-postMeansBetal
    postSDsBetal<-apply(betalSamps,c(1,3) , sd)
    restrictedGroupBetaSDs[,,i]<-postSDsBetal
    
    #post means Beta
    postMeansBETA<-colMeans(restricted$Beta)
    restrictedBetaMeans[,i]<-postMeansBETA
    postSDsBeta<-apply(restricted$Beta,2 , sd)
    restrictedBetaSDs[,i]<-postSDsBeta
    
    #Beta converge?
    restrictedBETAconverge[,i]<-abs(geweke.diag(mcmc(restricted$Beta))$z)
    
    
    
    #post means sigma2s
    postMeansSigma2s<-colMeans(restricted$sigma2s)
    restrictedSigma2Means[,i]<-postMeansSigma2s
    
    #post sds sigma2s
    postSDsSigma2s<-apply(restricted$sigma2s,2,sd)
    restrictedSigma2SDs[,i]<-postSDsSigma2s
    #sigma2s converge?
    restrictedSigma2Converge[,i]<-abs(geweke.diag(mcmc(restricted$sigma2s))$z)
    
    #bstar converge
    restrictedbstarConverge[i]<-abs(geweke.diag(mcmc(restricted$bstar))$z)
    #mu_rho converge
    restrictedMuRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$mu_rho))$z)
    #psi_rho_converge 
    restrictedPsiRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$psi_rho))$z)
    #rho_converge
    restrictedRhoConverge[i]<-abs(geweke.diag(mcmc(restricted$rho))$z)
    
    #-------------
    #Acceptance rates
    #------------
    #yAccept[,,i]<-restricted$yAccept #saving all 0/1s
    yAccept[i,]<-restricted$yAccept #just the column means
   # print(restricted$yAccept)
    #-------------------------
    #preds on holdout set
    #-------------------------
    postMeansBetalList<-split(postMeansBetal, rep(1:ncol(postMeansBetal), each = nrow(postMeansBetal)))
    restrictedPreds<-mapply(fits, postMeansBetalList,Xhold)
    restrictedPredListofLists[[i]]<-restrictedPreds
    
    
    #------------------------
    #computing marginal likelihoods for each element in holdout sample
    #------------------------
    restrictedMarginalsListofLists[[i]]<-fn.compute.marginals.hierModelNormal(betalSamps, restricted$sigma2s, yhold,Xhold)
    print(i)
  }
)

#yststRep20

#------




assign(paste0(regEst,"restrictedGroupBetaMeans",n, "_run", run),restrictedGroupBetaMeans)
assign(paste0(regEst,"restrictedGroupBetaSDs",n, "_run", run),    restrictedGroupBetaSDs)
assign(paste0(regEst,"restrictedbetalConverge",n, "_run", run),    restrictedbetalConverge)
assign(paste0(regEst,"restrictedBetaMeans",n, "_run", run),    restrictedBetaMeans)
assign(paste0(regEst,"restrictedBetaSDs",n, "_run", run),    restrictedBetaSDs)
assign(paste0(regEst,"restrictedBETAconverge",n, "_run", run),restrictedBETAconverge)
assign(paste0(regEst,"restrictedSigma2Means",n, "_run", run),restrictedSigma2Means)
assign(paste0(regEst,"restrictedSigma2SDs",n, "_run", run),restrictedSigma2SDs)
assign(paste0(regEst,"restrictedSigma2Converge",n, "_run", run),restrictedSigma2Converge)
assign(paste0(regEst,"restrictedbstarConverge",n, "_run", run),restrictedbstarConverge)
assign(paste0(regEst,"restrictedMuRhoConverge",n, "_run", run),restrictedMuRhoConverge)
assign(paste0(regEst,"restrictedPsiRhoConverge",n, "_run", run),restrictedPsiRhoConverge)
assign(paste0(regEst,"restrictedRhoConverge",n, "_run", run),restrictedRhoConverge)
assign(paste0(regEst,"restrictedPredListofLists",n, "_run", run),restrictedPredListofLists)
assign(paste0(regEst,"restrictedMarginalsListofLists",n, "_run", run),restrictedMarginalsListofLists)
assign(paste0(regEst,"rlmPredListofLists",n, "_run", run),rlmPredListofLists)
assign(paste0(regEst,"rlmMarginalsListofLists",n, "_run", run),rlmMarginalsListofLists)
assign(paste0(regEst,"restrictedstateCountsList",n,  "_run", run),stateCountsList)
assign(paste0(regEst,"yAccept",n,  "_run", run),yAccept)

assign(paste0("yholdoutMat",n,  "_run", run),get(paste0("yholdoutMat",n)))
assign(paste0("yOpenMat",n,  "_run", run),get(paste0("yOpenMat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))
assign(paste0("yOpenType1Mat",n,  "_run", run),get(paste0("yOpenType1Mat",n)))

assign(paste0("holdIndicesMatrix",n,  "_run", run),get(paste0("holdIndicesMatrix",n)))



rm(list=setdiff(ls(), c(ls()[grep('run',ls())], 'path.to.code','path.to.nonHierWorkSpaces','n','regEst')))



#--------------------------
#save the workspace
#--------------------------
save.image(paste0('../workSpaces/wsNwRestrictedHierarchAnalysis', regEst, '_n',n,'run', run, '.RData', sep=''))

