#Consolidate output from Hierarchical models; combining multiple workspaces from each model into one workspace

#analyses for n=1000,2000

#For n=1000, and 2000 I have  multiple (4) runs for each model, saved in aptley named workspaces
#--
#Load all of the results
#--
path.to.workspaces<-"../workSpaces/"
path.to.nonHierWorkSpaces<-"../../nwLinearReg/workSpaces/"
setwd(path.to.workspaces)

#load prior info; this also loads the analysis set called analysisSetCs
load(paste0(path.to.nonHierWorkSpaces,"nwdataPaper1PriorConstructionWorkSpace.RData"))
load("wsnwdataHierPriorConstruction.RData")

N<-nrow(analysisSetCs)

# n<-1000
# rest<-'Huber'

for(rest in c('Huber', 'Tukey')){
  for(run in 1:4){
    for(n in c(1000, 2000)){
  load(paste0('wsNwRestrictedHierarchAnalysis', rest, '_n', n, 'run', run, '.RData'))
  
}
}
}

for(run in 1:4){
   for(n in c(1000, 2000)){
  load(paste0('wsNwNormalHierarchAnalysis', '_n', n, 'run', run, '.RData'))
  load(paste0('wsNwTModelHierarchAnalysis', '_n', n, 'run', run, '.RData'))
}
}



for(n in c(1000,2000)){

 #huber restricted
  hubMarg<-c(get(paste0("HuberrestrictedMarginalsListofLists",n,"_run1")),
             get(paste0("HuberrestrictedMarginalsListofLists",n,"_run2")),
             get(paste0("HuberrestrictedMarginalsListofLists",n,"_run3")),
             get(paste0("HuberrestrictedMarginalsListofLists",n,"_run4"))
                 )
  
  #unlisting
  hubMarMat<-matrix(unlist(hubMarg), length(hubMarg),byrow=TRUE)
  
assign(paste0('All_HuberrestrictedMarginalsMat',n), hubMarMat)

  hubPred<-c(get(paste0("HuberrestrictedPredListofLists",n,"_run1")),
             get(paste0("HuberrestrictedPredListofLists",n,"_run2")),
             get(paste0("HuberrestrictedPredListofLists",n,"_run3")),
             get(paste0("HuberrestrictedPredListofLists",n,"_run4"))
             )

#unlisting
hubPredMat<-matrix(unlist(hubPred), length(hubPred),byrow=TRUE)
assign(paste0('All_HuberrestrictedPredMat',n), hubPredMat)

  #--------------------
#huber RLM
  #-------------------
  
  #huber rlm
  hubRlmMarg<-c(get(paste0("HuberrlmMarginalsListofLists",n,"_run1")),
             get(paste0("HuberrlmMarginalsListofLists",n,"_run2")),
             get(paste0("HuberrlmMarginalsListofLists",n,"_run3")),
             get(paste0("HuberrlmMarginalsListofLists",n,"_run4"))
  )
  
  #unlisting
  hubRlmMarMat<-matrix(unlist(hubRlmMarg), length(hubRlmMarg),byrow=TRUE)
  
  assign(paste0('All_HuberrlmMarginalsMat',n),  hubRlmMarMat)
  
  hubRlmPred<-c(get(paste0("HuberrlmPredListofLists",n,"_run1")),
             get(paste0("HuberrlmPredListofLists",n,"_run2")),
             get(paste0("HuberrlmPredListofLists",n,"_run3")),
             get(paste0("HuberrlmPredListofLists",n,"_run4"))
  )
  
  #unlisting
  hubRlmPredMat<-matrix(unlist(hubRlmPred), length(hubRlmPred),byrow=TRUE)
  assign(paste0('All_HuberrlmPredMat',n),hubRlmPredMat)  
  #------------------------------------------------------
  

#--------------------------------------------------------------------
  tukMarg<-c(get(paste0("TukeyrestrictedMarginalsListofLists",n,"_run1")),
           get(paste0("TukeyrestrictedMarginalsListofLists",n,"_run2")),
           get(paste0("TukeyrestrictedMarginalsListofLists",n,"_run3")),
           get(paste0("TukeyrestrictedMarginalsListofLists",n,"_run4"))
             )
#unlisting
tukMarMat<-matrix(unlist(tukMarg), length(tukMarg),byrow=TRUE)
  assign(paste0('All_TukeyrestrictedMarginalsMat',n), tukMarMat)

  tukPred<-c(get(paste0("TukeyrestrictedPredListofLists",n,"_run1")),
           get(paste0("TukeyrestrictedPredListofLists",n,"_run2")),
           get(paste0("TukeyrestrictedPredListofLists",n,"_run3")),
           get(paste0("TukeyrestrictedPredListofLists",n,"_run4"))
             )
#unlisting
tukPredMat<-matrix(unlist(tukPred), length(tukPred),byrow=TRUE)
  assign(paste0('All_TukeyrestrictedPredMat',n), tukPredMat)

#Tukey RLM
  
  tukRlmMarg<-c(get(paste0("TukeyrlmMarginalsListofLists",n,"_run1")),
                get(paste0("TukeyrlmMarginalsListofLists",n,"_run2")),
                get(paste0("TukeyrlmMarginalsListofLists",n,"_run3")),
                get(paste0("TukeyrlmMarginalsListofLists",n,"_run4"))
  )
  
  #unlisting
  tukRlmMarMat<-matrix(unlist(tukRlmMarg), length(tukRlmMarg),byrow=TRUE)
  
  assign(paste0('All_TukeyrlmMarginalsMat',n),  tukRlmMarMat)
  
  tukRlmPred<-c(get(paste0("TukeyrlmPredListofLists",n,"_run1")),
                get(paste0("TukeyrlmPredListofLists",n,"_run2")),
                get(paste0("TukeyrlmPredListofLists",n,"_run3")),
                get(paste0("TukeyrlmPredListofLists",n,"_run4"))
  )
  
  #unlisting
  tukRlmPredMat<-matrix(unlist(tukRlmPred), length(tukRlmPred),byrow=TRUE)
  assign(paste0('All_TukeyrlmPredMat',n),tukRlmPredMat)  

#-----------------------------------------------------------------  


#Normal Theory
nTheoryMarg<-c(get(paste0("nTheoryMarginalsListofLists",n,"_run1")),
           get(paste0("nTheoryMarginalsListofLists",n,"_run2")),
           get(paste0("nTheoryMarginalsListofLists",n,"_run3")),
           get(paste0("nTheoryMarginalsListofLists",n,"_run4"))
)

#unlisting
nTheoryMarMat<-matrix(unlist(nTheoryMarg), length(nTheoryMarg),byrow=TRUE)

assign(paste0('All_nTheoryMarginalsMat',n), nTheoryMarMat)


nTheoryPred<-c(get(paste0("nTheoryPredListofLists",n,"_run1")),
           get(paste0("nTheoryPredListofLists",n,"_run2")),
           get(paste0("nTheoryPredListofLists",n,"_run3")),
           get(paste0("nTheoryPredListofLists",n,"_run4"))
)
nTheoryPredMat<-matrix(unlist(nTheoryPred), length(nTheoryPred),byrow=TRUE)

assign(paste0('All_nTheoryPredMat',n), nTheoryPredMat)

#OLS
  #Normal Theory
  olsMarg<-c(get(paste0("olsMarginalsListofLists",n,"_run1")),
                 get(paste0("olsMarginalsListofLists",n,"_run2")),
                 get(paste0("olsMarginalsListofLists",n,"_run3")),
                 get(paste0("olsMarginalsListofLists",n,"_run4"))
  )
  
  #unlisting
  olsMarMat<-matrix(unlist(olsMarg), length(olsMarg),byrow=TRUE)
  
  assign(paste0('All_olsMarginalsMat',n), olsMarMat)
  
  
  olsPred<-c(get(paste0("olsPredListofLists",n,"_run1")),
                 get(paste0("olsPredListofLists",n,"_run2")),
                 get(paste0("olsPredListofLists",n,"_run3")),
                 get(paste0("olsPredListofLists",n,"_run4"))
  )
  olsPredMat<-matrix(unlist(olsPred), length(olsPred),byrow=TRUE)
  
  assign(paste0('All_olsPredMat',n), olsPredMat)
  
  #----------------------
  #same for t model
  #---------------------  
  
  tModelMarg<-c(get(paste0("tModelMarginalsListofLists",n,"_run1")),
                 get(paste0("tModelMarginalsListofLists",n,"_run2")),
                 get(paste0("tModelMarginalsListofLists",n,"_run3")),
                 get(paste0("tModelMarginalsListofLists",n,"_run4"))
  )
  
  #unlisting
  tModelMarMat<-matrix(unlist(tModelMarg), length(tModelMarg),byrow=TRUE)
  
  assign(paste0('All_tModelMarginalsMat',n), tModelMarMat)
  
  
  tModelPred<-c(get(paste0("tModelPredListofLists",n,"_run1")),
                 get(paste0("tModelPredListofLists",n,"_run2")),
                 get(paste0("tModelPredListofLists",n,"_run3")),
                 get(paste0("tModelPredListofLists",n,"_run4"))
  )
  tModelPredMat<-matrix(unlist(tModelPred), length(tModelPred),byrow=TRUE)
  
  assign(paste0('All_tModelPredMat',n), tModelPredMat)

  
  
#-------------------------------------------------------- 
  
  
  
  

#combine HoldIndecesMartix
#Note: don't use original hold Indices matrix from non-hier fit: this will have the same indices just not in the same order
# hold<-rbind(get(paste0("holdIndicesMatrix",n, "_run1")),
#                get(paste0("holdIndicesMatrix",n, "_run2")),
#                get(paste0("holdIndicesMatrix",n, "_run3")),
#                get(paste0("holdIndicesMatrix",n, "_run4"))
# )
# assign(paste0("All_holdIndicesMatrix",n), hold)

holdHierMatrix1<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix2<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix3<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix4<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix5<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix6<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix7<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix8<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix9<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix10<-matrix(NA,length(hubMarg),N-n)
holdHierMatrix11<-matrix(NA,length(hubMarg),N-n)

for(run in 1:nrow(holdHierMatrix1)){
holdHierMatrix1[run,]<-as.numeric(unlist(lapply(hubMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix2[run,]<-as.numeric(unlist(lapply(hubPred[[run]], rownames), use.names=FALSE)) 
#holdHierMatrix3[run,]<-as.numeric(unlist(lapply(hubRlmMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix3[run,]<-as.numeric(unlist(lapply(hubRlmPred[[run]], rownames), use.names=FALSE)) 

holdHierMatrix4[run,]<-as.numeric(unlist(lapply(tukMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix5[run,]<-as.numeric(unlist(lapply(tukPred[[run]], rownames), use.names=FALSE)) 
#holdHierMatrix7[run,]<-as.numeric(unlist(lapply(tukRlmMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix6[run,]<-as.numeric(unlist(lapply(tukRlmPred[[run]], rownames), use.names=FALSE)) 

holdHierMatrix7[run,]<-as.numeric(unlist(lapply(nTheoryMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix8[run,]<-as.numeric(unlist(lapply(nTheoryPred[[run]], rownames), use.names=FALSE)) 

#holdHierMatrix9[run,]<-as.numeric(unlist(lapply(olsMarg[[run]], rownames), use.names=FALSE)) 
holdHierMatrix9[run,]<-as.numeric(unlist(lapply(olsPred[[run]], rownames), use.names=FALSE)) 

holdHierMatrix10[run,]<-as.numeric(unlist(lapply(tModelMarg[[run]], names), use.names=FALSE)) 
holdHierMatrix11[run,]<-as.numeric(unlist(lapply(tModelPred[[run]], rownames), use.names=FALSE)) 

}
if(sum(is.na(holdHierMatrix1)>0)) stop('hold matrix NA')
   for(i in 1:11){
if(!all.equal(holdHierMatrix1,get(paste0("holdHierMatrix",i)))) { 
  stop('hold matrix not equal')
} else {print(i)}
}

assign(paste0("All_holdHierIndexMat", n), holdHierMatrix1)

}
   
   
save.image('consolidateOutputFromHierAnalysesFIX.RData')
