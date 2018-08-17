#######################
#John Lewis
#Summarize the output from nwAnalysisSimulation.R
#Task: pull in results from simulations and caluclate the marginal likelihoods 
#######################

#For n=1000 and 2000 I have  multiple (4) runs, saved in aptley named workspaces

#library(tables) #for easy importing of tabels to latex



#------------------------------------------
#Load all of the results
#-------------------------------------------



setwd("../workSpaces")

load("nwdataPaper1PriorConstructionWorkSpace.RData")


for(n in c(25,50,100,200,500)){
  load(paste0('wsNwAnalysis_n',n,'.RData'))
}
#--------------------------------------------
#load each run, save results and combine, move to next n
#------------------------------------------

#for n=1000, 2000, I have four seperate runs

#objects that will be combined via rbind
saveNames1<-c('yholdoutMat','yOpenMat','yOpenType1Mat','margNTheory','margRest',"margRestHuber",'margT','margRlm',"margRlmHuber",'margOls','holdIndicesMatrix',"rlmPredMat", 'rlmPredHuberMat','olsPredMat','nTheoryPredMat','restPredMat','restPredHuberMat','tPredMat')
#objects that will be combined via cbind
saveNames2<-c('rlmEstimates',"rlmEstimatesHuber",'nTheoryEstimates','restrictedEstimates',"restrictedEstimatesHuber",'tmodelEstimates','olsEstimates')
#objects that will be combined via append
saveNames3<-c("acceptT", "acceptY", "acceptYHuber")


for(n in c(1000, 2000)){ #maybe add 2500
  
  #initialize
  for(sv in 1:length(saveNames1)){
    assign(saveNames1[sv], vector('list', 1))
  }
  #initialize
  for(sv in 1:length(saveNames2)){
    assign(saveNames2[sv], vector('list', 1))
  }
  #initialize
  for(sv in 1:length(saveNames2)){
    assign(saveNames2[sv], vector('list', 1))
  }
  
  #for run=1
  run<-1
  load(paste0('wsNwAnalysis_n', n, 'run', run, '.RData'))
  for(sv in saveNames1){
    assign(sv, get(paste0(sv,n)))
  }
  for(sv in saveNames2){
    assign(sv, get(paste0(sv,n)))
  }
  for(sv in saveNames3){
    assign(sv, get(paste0(sv,n)))
  }
  #combine with run=2:4
  for(run in 2:4){  
    load(paste0('wsNwAnalysis_n', n, 'run', run, '.RData'))
    for(sv in saveNames1){
      assign(sv,rbind(get(sv), get(paste0(sv,n))))
    }
    for(sv in saveNames2){
      assign(sv,cbind(get(sv), get(paste0(sv,n))))
    }
    for(sv in saveNames3){
      assign(sv,append(get(sv), get(paste0(sv,n))))
    }
  }
  #append n onto all of saveNames
  for(sv in c(saveNames1, saveNames2, saveNames3)){
    assign(paste0(sv,n), get(sv))
  }
}

reps<-100
N<-3303

save.image('../workSpaces/consolidatedOutputForPlots.RData')




#Not run
# #----------------------------------------
# #summarizing the results
# #----------------------------------------
# 
# 
# 
# rounding<-5 #number of decimal places to round
# multiply<-100 #multiply mse's by this number
# model.names<-c('rlm-T','rlm-H','ols', 'normal', 'rest.-T','rest.-H', 't') #names in the original order
# ord<-c(1,5,2,6,7,3,4) #the order to display the results in latex
# #put in the order I want
# model.names<-model.names[ord]
# sampsizes<-c('25','50','100','200','500','1000','1500','2000')
# 
# predNames<-c("rlmPredMat", 'rlmPredHuberMat','olsPredMat','nTheoryPredMat','restPredMat','restPredHuberMat','tPredMat')
# #put in the order I want
# predNames<-predNames[ord]
# 
# margNames<-c('margRlm',"margRlmHuber",'margOls', 'margNTheory','margRest',"margRestHuber",'margT')
# #put in the order I want
# margNames<-margNames[ord]
# cbind(margNames, predNames)
# 
# acceptNames<-c('acceptY','acceptYHuber','acceptT')
# 
# 

