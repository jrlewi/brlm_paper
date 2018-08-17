
#This code produces the plots in the paper.
#Summarize the output from nwAnalysisSimulation.R

#analyses for n=25,50,100,200,500,1000,2000

#For n=1000, and 2000 I have  multiple (4) runs saved in aptley named workspaces: note: objects in each workspace of the same name: must load each ws by itself, transfer results, then load next ws.
rm(list=ls())
#add your own path to save figures. To just make figures without saving find and replace 'pdf' with '#pdf', 'dev.off' with '#dev.off' and "setwd(path.to.figures)" with "#setwd(path.to.figures)"


#--
#Load all of the results
#--
path.to.workspaces<-"../workSpaces"
setwd(path.to.workspaces)

load("nwdataPaper1PriorConstructionWorkSpace.RData")

#load summary results
load('consolidatedOutputForPlots.RData')

path.to.figures<-"/Users/johnlewis/Box Sync/research/snm/BlendedParadigmWriting/paper_version2/figures"
path.to.figures<-"~/Box Sync/research/snm/BlendedParadigmWriting/paper_version2/figures"
#path.to.figures<-"your path here"


rounding<-5 #number of decimal places to round
multiply<-100 #multiply mse's by this number
model.names<-c('rlm-T','rlm-H','ols', 'normal', 'incomp.-T','incomp.-H', 't') #names in the original order
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


#--
#which indices are open, and open type1
zeroOnCsScale<-min(analysisSetCs$sqrt_Count2012) #closed agencies count on masked data
openAgencies<-which(analysisSetCs$sqrt_Count2012>zeroOnCsScale)
openType1<-which(analysisSetCs$sqrt_Count2012>zeroOnCsScale & analysisSetCs$Agency_Type=='Type1')
type1Agencies<-which(analysisSetCs$Agency_Type=='Type1')
length(type1Agencies)/nrow(analysisSetCs)
head(openType1)
head(openAgencies)
#---

#trimming proportions
trim.vals<-c(0,.05,.10,.15,.20,.25,.30)



#check that these match
cbind(predNames,margNames)

margNames[c(2,4,5,7)]

#------------------------------
#for each fitting and holout set trim the marginal likelihoods based on a base model: i.e. trim the same indices as the base model. calculate the mean and sd of the mean marginal likelihoods across the different fitting/holdout sets
#---------

fn.get.mean.sd.marg<-function(agenciesToUse, n, trim.val,marg.names, base.model){
  hold<-get(paste0("holdIndicesMatrix",n))
  allMeans<-array(NA, dim=c(nrow(hold), length(marg.names), length(trim.val)))
  inScopeAgencies<-apply(hold, 1,function(x) which(x %in% agenciesToUse))
  #force to be a list
  if(class(inScopeAgencies)=='matrix'){
    inScopeAgencies<-lapply(seq_len(ncol(inScopeAgencies)), function(i) inScopeAgencies[,i])
  }
  for(rep in 1:nrow(hold)){
    marg<-log(get(paste0(base.model,n))[rep,inScopeAgencies[[rep]]]) #on log scale
    for(tm in trim.val){
    include<-which(marg>=quantile(marg, tm))
      for(mn in marg.names){
        tmp<-log(get(paste0(mn, n))[rep,inScopeAgencies[[rep]]][include]) #put on log scale
        allMeans[rep, which(marg.names==mn), which(trim.val==tm)]<-mean(tmp[tmp!=-Inf])
      }
    }
  }
  out<-list()
  mns<-apply(allMeans, c(2,3), mean)
  sds<-apply(allMeans, c(2,3), sd)
  colnames(mns)<-trim.val
  rownames(mns)<-marg.names
  colnames(sds)<-trim.val
  rownames(sds)<-marg.names
  out$means<-mns
  out$sds<-sds
  out
}



tst3<-fn.get.mean.sd.marg(1:N, n=2000, trim.val=trim.vals,marg.names=margNames, base.model=margNames[7])
tst3

margNames
model.names
fn.plot.mean.st.marginals<-function(agenciesToUse,agenciesLab, n, trims,whch,base, ylim=NULL, xlim=NULL, legend=TRUE,cex=1,pch=19:25, cex.lab=1, cex.main=1){
  par(mar=c(5.1, 4.4, 4.1, 2.1))
  mn.sd<-fn.get.mean.sd.marg(agenciesToUse, n, trims,margNames[whch], margNames[base])
  lg<-length(trims)
  jitters<-seq(-1,1, length.out=length(margNames[whch]))
  trimMat<-matrix(NA, length(trims), length(margNames[whch]))
  for(i in 1:length(margNames[whch])){
    trimMat[,i]<-seq(from=1,to=4*length(trims), by=4)+jitters[i]
  }
  x<-trimMat
  y<-t(mn.sd$means)
  sd<-t(mn.sd$sds)
  if(is.null(ylim)) {ylim=c(min(y-sd),max(y+sd))}
  matplot(trimMat,t(mn.sd$means), type='p',pch=pch,col=color[whch], bg=color[whch], xaxt='n', xlab=paste0('trimming fraction'), ylab=paste0('mean log marginals'), lty=1:lg, font.main = 3, main=paste('n=', n, sep=""),ylim=ylim,xlim=xlim, cex.lab=cex.lab,cex.main=cex.main)
  abline(v=trimMat[-nrow(trimMat),ncol(trimMat)]+.5*(trimMat[2,1]-trimMat[1,ncol(trimMat)]), col='gray',lty=2)
  axis(1, at=seq(from=1,to=4*length(trims), by=4), labels=trims)
 if(legend) {legend('topleft', legend=model.names[whch],pch=pch, col=color[whch],pt.bg=color[whch], cex=cex,xpd=TRUE)} #lty=1: inset=c(-0.2,0),
  segments(x, y-sd,x, y+sd)
#  epsilon = 0.02
#  segments(x-epsilon,y-sd,x+epsilon,y-sd)
#  segments(x-epsilon,y+sd,x+epsilon,y+sd)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}
par()$mar

#-------------------------------------------------------
#plots
#-------------------------------------------------------
setwd(path.to.figures)

color<-c('magenta', 'red', 'green','blue', 'dark green', 'gray','black')

lg<-1:length(color)



sizes<-c(25,100,1000,2000) #training set sizes to plot
baseModelNumber<-5
margModelNumbers<-1:5#c(2,4,5,6,7)
trim.vals<-c(0,0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
margNames
ylims<-c(-1,2)

pdf(paste0('logMargAllAgencies','BaseModel', model.names[baseModelNumber],'.pdf'), width=8)
n<-sizes[1]
mn<-margNames[baseModelNumber]
fn.plot.mean.st.marginals(1:N,"All Agencies", n, trims=trim.vals, whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=TRUE,cex=1)

for(n in sizes[-1]){
  mn<-margNames[baseModelNumber]
    fn.plot.mean.st.marginals(1:N,"All Agencies", n, trims=trim.vals, whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=FALSE, cex=1)
  }
dev.off()

pdf(paste0('logMargOpenAgencies','BaseModel', model.names[baseModelNumber],'.pdf'), width=8)
n<-sizes[1]
mn<-margNames[baseModelNumber]
fn.plot.mean.st.marginals(openAgencies,"Open Agencies", n, trims=trim.vals, whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=TRUE,cex=1)
for(n in sizes[-1]){
  mn<-margNames[baseModelNumber]
    fn.plot.mean.st.marginals(openAgencies,"Open Agencies", n, trims=trim.vals,whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=FALSE,cex=1)
}
dev.off()

pdf(paste0('logMargType1Agencies','BaseModel', model.names[baseModelNumber],'.pdf'), width=8)
labSize<-2.15
mainSize<-2
legSize<-1.7
n<-sizes[1]
mn<-margNames[baseModelNumber]
fn.plot.mean.st.marginals(type1Agencies,"Type 1 Agencies",  n, trims=trim.vals, whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=TRUE,cex=legSize, cex.lab=labSize, cex.main=mainSize)

for(n in sizes[-1]){
  mn<-margNames[baseModelNumber]
  fn.plot.mean.st.marginals(type1Agencies,"Type 1 Agencies", n, trims=trim.vals,whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=FALSE,cex=legSize, cex.lab=labSize, cex.main=mainSize)
}  
dev.off()



pdf(paste0('logMargOpenType1Agencies','BaseModel', model.names[baseModelNumber],'.pdf'), width=8)

n<-sizes[1]
mn<-margNames[baseModelNumber]
fn.plot.mean.st.marginals(openType1,"Open, Type 1 Agencies",  n, trims=trim.vals, whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=TRUE,cex=1)
for(n in sizes[-1]){
  mn<-margNames[baseModelNumber]
    fn.plot.mean.st.marginals(openType1,"Open, Type 1 Agencies", n, trims=trim.vals,whch=margModelNumbers, base=baseModelNumber, ylim=ylims, legend=TRUE,cex=1)
}  
dev.off()




