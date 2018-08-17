#nw Analyze: compare the non-hierarchical and hierarchical fits:

#install.packages('tables')
library(tables) #for easy importing of tabels to latex
rm(list=ls())

#non-hierarchical

path.to.nonHierworkspaces<-"../../nwLinearReg/workSpaces/"
setwd(path.to.nonHierworkspaces)
load("nwdataPaper1PriorConstructionWorkSpace.RData")
load('consolidatedOutputForPlots.RData')


#hierarchical
path.to.workspaces<-"../../nwHierLinearReg/workSpaces/"
setwd(path.to.workspaces)
load("wsnwdataHierPriorConstruction.RData")
#load consolidated WS results
load('consolidateOutputFromHierAnalysesFIX.RData') #fixed the hyper params


rounding<-5 #number of decimal places to round
multiply<-100 #multiply mse's by this number
model.names<-c('rlm-T','rlm-H','ols', 'normal', 'rest.-T','rest.-H', 't') #names in the original order
ord<-c(1,5,2,6,7,3,4) #the order to display the results in latex
#put in the order I want
model.names<-model.names[ord]
sampsizes<-c('25','50','100','200','500','1000','1500','2000')

predNames<-c("rlmPredMat", 'rlmPredHuberMat','olsPredMat','nTheoryPredMat','restPredMat','restPredHuberMat','tPredMat')
#put in the order I want
predNames<-predNames[ord]

margNames<-c('margRlm',"margRlmHuber",'margOls', 'margNTheory','margRest',"margRestHuber",'margT')
#put in the order I want
margNames<-margNames[ord]
cbind(margNames, predNames)

acceptNames<-c('acceptY','acceptYHuber','acceptT')


predNamesHier<-c("All_TukeyrlmPredMat", 'All_HuberrlmPredMat','All_olsPredMat','All_nTheoryPredMat','All_TukeyrestrictedPredMat','All_HuberrestrictedPredMat','All_tModelPredMat')
#put in the order I want
predNamesHier<-predNamesHier[ord]

margNamesHier<-c('All_TukeyrlmMarginalsMat',"All_HuberrlmMarginalsMat",'All_olsMarginalsMat', 'All_nTheoryMarginalsMat','All_TukeyrestrictedMarginalsMat',"All_HuberrestrictedMarginalsMat",'All_tModelMarginalsMat')
#put in the order I want
margNamesHier<-margNamesHier[ord]
cbind(margNames,margNamesHier, predNames,predNamesHier,model.names)


#--
#which indices are open, and open type1
zeroOnCsScale<-min(analysisSetCs$sqrt_Count2012) 
openAgencies<-which(analysisSetCs$sqrt_Count2012>zeroOnCsScale)
openType1<-which(analysisSetCs$sqrt_Count2012>zeroOnCsScale & analysisSetCs$Agency_Type=='Type1')
type1Agencies<-which(analysisSetCs$Agency_Type=='Type1')
length(type1Agencies)/nrow(analysisSetCs)
head(openType1)
head(openAgencies)
#---

#trimming proportions
trim.vals<-c(0,.05,.10,.15,.20,.25,.30)


#returns the vector of 100 mean TLM for a single n and one marg.name
#note need a different (slighty) function for Hier and non-hier: the holdIndices Matrices are in different orders for Hierarchical versus non hierarchical(they are the same: just in different orders); I checked this with the previous fn.get.mean.sd.marg function defined in the respective analyzeOutputWithin.R files
fn.get.mean.sd.margHier<-function(agenciesToUse,n, trim.val,marg.name, base.model){
  
  hold<-get(paste0("All_holdHierIndexMat",n))
  allMeans<-array(NA, dim=c(nrow(hold), length(trim.val)))
  inScopeAgencies<-apply(hold, 1,function(x) which(x %in% agenciesToUse))
  #force to be a list
  if(class(inScopeAgencies)=='matrix'){
    inScopeAgencies<-lapply(seq_len(ncol(inScopeAgencies)), function(i) inScopeAgencies[,i])
  }
  for(rep in 1:nrow(hold)){
    marg<-log(get(paste0(base.model,n))[rep,inScopeAgencies[[rep]]]) #on log scale
    for(tm in trim.val){
      include<-which(marg>=quantile(marg, tm))
      for(mn in marg.name){
        tmp<-log(get(paste0(mn, n))[rep,inScopeAgencies[[rep]]][include]) #put on log scale
        allMeans[rep, which(trim.val==tm)]<-mean(tmp[tmp!=-Inf])
      }
    }
  }
  out<-list()
  colnames(allMeans)<-trim.val
  mns<-apply(allMeans, c(2), mean)
  sds<-apply(allMeans, c(2), sd)
  names(mns)<-trim.val
  colnames(allMeans)<-trim.val
  #rownames(mns)<-marg.name
  names(sds)<-trim.val
  #rownames(sds)<-marg.name
  out$allMeans<- allMeans
  out$means<-mns
  out$sds<-sds
  out
}


fn.get.mean.sd.margNon<-function(agenciesToUse, n, trim.val,marg.name, base.model){
  hold<-get(paste0("holdIndicesMatrix",n))
  allMeans<-array(NA, dim=c(nrow(hold), length(trim.val)))
  inScopeAgencies<-apply(hold, 1,function(x) which(x %in% agenciesToUse))
  #force to be a list
  if(class(inScopeAgencies)=='matrix'){
    inScopeAgencies<-lapply(seq_len(ncol(inScopeAgencies)), function(i) inScopeAgencies[,i])
  }
  for(rep in 1:nrow(hold)){
    marg<-log(get(paste0(base.model,n))[rep,inScopeAgencies[[rep]]]) #on log scale
    for(tm in trim.val){
      include<-which(marg>=quantile(marg, tm))
      for(mn in marg.name){
        tmp<-log(get(paste0(mn, n))[rep,inScopeAgencies[[rep]]][include]) #put on log scale
        allMeans[rep, which(trim.val==tm)]<-mean(tmp[tmp!=-Inf])
      }
    }
  }
  out<-list()
  colnames(allMeans)<-trim.val
  mns<-apply(allMeans, c(2), mean)
  sds<-apply(allMeans, c(2), sd)
  names(mns)<-trim.val
  #rownames(mns)<-marg.names
  names(sds)<-trim.val
  #rownames(sds)<-marg.names
  out$allMeans<-allMeans
  out$means<-mns
  out$sds<-sds
  out
}


fn.pairedTtest<-function(allMeans1,allMeans2){
  if(ncol(allMeans1)!=ncol(allMeans2)){stop('not same dimension')}
ttestList<-list(); length(ttestList)<-ncol(allMeans1) 
for(i in 1:ncol(allMeans1)){
  ttestList[[i]]<-t.test(allMeans1[,i],allMeans2[,i], paired=TRUE)
}
ttestList  
}

n<-1000
trim.vals<-c(0,.05,.10,.15,.20,.25,.30)
trim.vals<-c(.15,.20,.25,.30)
#concentrate on type1 agencies and Tukey's restricted model model

#tukey
nonHierN1000<-fn.get.mean.sd.margNon(type1Agencies, n=1000, trim.vals,margNames[2], margNames[5])

nonHierN2000<-fn.get.mean.sd.margNon(type1Agencies, n=2000, trim.vals,margNames[2], margNames[5])

HierN1000<-fn.get.mean.sd.margHier(type1Agencies, n=1000, trim.vals,margNamesHier[2], margNamesHier[5])

HierN2000<-fn.get.mean.sd.margHier(type1Agencies, n=2000, trim.vals,margNamesHier[2], margNamesHier[5])


#huber
nonHierHubN1000<-fn.get.mean.sd.margNon(type1Agencies, n=1000, trim.vals,margNames[4], margNames[5])

nonHierHubN2000<-fn.get.mean.sd.margNon(type1Agencies, n=2000, trim.vals,margNames[4], margNames[5])

HierHubN1000<-fn.get.mean.sd.margHier(type1Agencies, n=1000, trim.vals,margNamesHier[4], margNamesHier[5])

HierHubN2000<-fn.get.mean.sd.margHier(type1Agencies, n=2000, trim.vals,margNamesHier[4], margNamesHier[5])

all.equal(rbind(colMeans(HierN1000$allMeans),colMeans(nonHierN1000$allMeans))
,rbind(HierN1000$means,nonHierN1000$means))
rbind(colMeans(HierN2000$allMeans),colMeans(nonHierN2000$allMeans))

pairedTn1000<-fn.pairedTtest(nonHierN1000$allMeans,HierN1000$allMeans)
pairedTn2000<-fn.pairedTtest(nonHierN2000$allMeans,HierN2000$allMeans)

lapply(pairedTn1000, FUN=function(x) c(x$estimate, x$p.value))
lapply(pairedTn2000, FUN=function(x) c(x$estimate, x$p.value))


#-Normal model
#tukey
nonHierNormalN1000<-fn.get.mean.sd.margNon(type1Agencies, n=1000, trim.vals,margNames[7], margNames[5])

nonHierNormalN2000<-fn.get.mean.sd.margNon(type1Agencies, n=2000, trim.vals,margNames[7], margNames[5])

HierNormalN1000<-fn.get.mean.sd.margHier(type1Agencies, n=1000, trim.vals,margNamesHier[7], margNamesHier[5])

HierNormalN2000<-fn.get.mean.sd.margHier(type1Agencies, n=2000, trim.vals,margNamesHier[7], margNamesHier[5])



###################################
#Combining to latex tables
###################################
# path.to.tables<-"C:/Users/lewisjr2/Box Sync/research/snm/BlendedParadigmWriting/paper_version2"
# setwd(path.to.tables)

#tukey
msdNon1 <- paste(format(round(nonHierN1000$means,3), dropTrailing=FALSE)," (",format(round(nonHierN1000$sds,3), dropTrailing=FALSE),")",sep="")
msdHier1 <- paste(format(round(HierN1000$means,3), dropTrailing=FALSE)," (",format(round(HierN1000$sds,3), dropTrailing=FALSE),")",sep="")

msdNon2 <- paste(format(round(nonHierN2000$means,3), dropTrailing=FALSE)," (",format(round(nonHierN2000$sds,3), dropTrailing=FALSE),")",sep="")
msdHier2 <- paste(format(round(HierN2000$means,3), dropTrailing=FALSE)," (",format(round(HierN2000$sds,3), dropTrailing=FALSE),")",sep="")

#huber

msdHubNon1 <- paste(format(round(nonHierHubN1000$means,3), dropTrailing=FALSE)," (",format(round(nonHierHubN1000$sds,3), dropTrailing=FALSE),")",sep="")
msdHubHier1 <- paste(format(round(HierHubN1000$means,3), dropTrailing=FALSE)," (",format(round(HierHubN1000$sds,3), dropTrailing=FALSE),")",sep="")

msdHubNon2 <- paste(format(round(nonHierHubN2000$means,3), dropTrailing=FALSE)," (",format(round(nonHierHubN2000$sds,3), dropTrailing=FALSE),")",sep="")
msdHubHier2 <- paste(format(round(HierHubN2000$means,3), dropTrailing=FALSE)," (",format(round(HierHubN2000$sds,3), dropTrailing=FALSE),")",sep="")



tab<-rbind(msdNon1,msdHier1,msdNon2,msdHier2,msdHubNon1,msdHubHier1,msdHubNon2,msdHubHier2)
# tab<-round(tab, 3)
# tab<-rbind(nonHierN1000$means,HierN1000$means,nonHierN2000$means,HierN2000$means)
# tab<-round(tab, 3)
rownames(tab)<-rep(c('Non-Hier.','Hier.'),4)
tab
colnames(tab)<-trim.vals
getwd()
table_options(justification='c')
tab1<-latex(tab,title='$f$',
            label='tlmTable',
            rgroup=c("Tukey ($n=1000$)","Tukey ($n=2000$) ","Huber  ($n=1000$)","Huber ($n=2000$)"),
            n.rgroup=c(2,2,2,2),
            caption.loc='bottom',
            math.col.names=TRUE,
            rgroupTexCmd='mdseries',
            size="small",
            file='tlmTable2.tex',
            caption='Mean (standard deviation) of TLM for `Type 1\' agencies for the restricted non-hierarchical and hierarchical models for $n=1000$ and $2000$. The table is used for comparision of the two methods and a summary appears in the text.' 
            )
#tlmTable2, fized a slight bug in the restricted model fitting function...didn't change results much.



##Adding the Guassian Values
#Normal n=1000
msdNormalNon1 <- paste(format(round(nonHierNormalN1000$means,3), dropTrailing=FALSE)," (",format(round(nonHierNormalN1000$sds,3), dropTrailing=FALSE),")",sep="")
msdNormalHier1 <- paste(format(round(HierNormalN1000$means,3), dropTrailing=FALSE)," (",format(round(HierNormalN1000$sds,3), dropTrailing=FALSE),")",sep="")

#normal n=2000
msdNormalNon2 <- paste(format(round(nonHierNormalN2000$means,3), dropTrailing=FALSE)," (",format(round(nonHierNormalN2000$sds,3), dropTrailing=FALSE),")",sep="")
msdNormalHier2 <- paste(format(round(HierNormalN2000$means,3), dropTrailing=FALSE)," (",format(round(HierNormalN2000$sds,3), dropTrailing=FALSE),")",sep="")




tab2<-rbind(msdNon1,msdHier1,msdNon2,msdHier2,msdHubNon1,msdHubHier1,msdHubNon2,msdHubHier2,msdNormalNon1,msdNormalHier1,msdNormalNon2,msdNormalHier2)

rownames(tab2)<-rep(c('Non-Hier.','Hier.'),6)
tab2
colnames(tab2)<-trim.vals
getwd()
table_options(justification='c')
tab22<-latex(tab2,title='$f$',
            label='tlmTableWithGaus',
            rgroup=c("Tukey ($n=1000$)","Tukey ($n=2000$) ","Huber  ($n=1000$)","Huber ($n=2000$)","Gaussian  ($n=1000$)","Gaussian ($n=2000$)"),
            n.rgroup=c(2,2,2,2,2,2),
            caption.loc='bottom',
            math.col.names=TRUE,
            rgroupTexCmd='mdseries',
            size="small",
            file='tlmTable2WithGaus.tex',
            caption='Mean (standard deviation) of TLM for `Type 1\' agencies for the restricted non-hierarchical and hierarchical models for $n=1000$ and $2000$. The table is used for comparision of the two methods and a summary appears in the text. Also supplied are the values for the Gaussian model fits' 
)





