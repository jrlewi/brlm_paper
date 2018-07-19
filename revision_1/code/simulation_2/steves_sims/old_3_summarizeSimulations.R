#
# Summaries of Evaluation Metrics ---
#
library('plyr')
library('reshape2')
rm(list=ls())
concTau2Near1<-0

#if(concTau2Near1){
#  load('simPostProcessed_concTau2Near1.RData')
#  #load('simPostProcessed2_concTau2Near1.RData')
#} else {
#  load('simPostProcessed.RData')
#  #load('simPostProcessed2.RData')
#}
load('//trad/profs$/maceachern.1/Desktop/LEWIS.LEE/codeForPaper/codeForPaper/StevesimPostProc.RData') # SteveMod


#load('simPostProcessed2_sampleVersion.RData')
sumNames<-c('Full', 'RestHuber', 'RestTukey', 'RlmHuber','RlmTukey')
#get(paste0('mse', sumNames[1],'Sig2_', sigma2True[1]))


#combine the results into data.frame
mseDf<-matrix(NA, ncol=9)
klDf<-matrix(NA, ncol=9)
#trimMspeDf<-matrix(NA, ncol=9)
#tlmDf<-matrix(NA, ncol=9)

for(i in 1:length(sigma2True)){
mseDfTemp<-cbind(factorsMat, sigma2True[i])
klDfTemp<-cbind(factorsMat, sigma2True[i])
#trimMspeDfTemp<-cbind(factorsMat, sigma2True[i])
#tlmDfTemp<-cbind(factorsMat, sigma2True[i])

#meanTlmDfTemp<-cbind(factorsMat, sigma2True[i])
 for(j in 1:5){
   mseDfTemp<-cbind(mseDfTemp, get(paste0('mse', sumNames[j],'Sig2_', sigma2True[i])))
   klDfTemp<-cbind(klDfTemp, get(paste0('kl', sumNames[j],'Sig2_', sigma2True[i])))
 #  trimMspeDfTemp<-cbind(trimMspeDfTemp, get(paste0('trimMspe', sumNames[j],'Sig2_', sigma2True[i]))/sigma2True[i]) 
  # tlmDfTemp<-cbind(tlmDfTemp, get(paste0('tlm', sumNames[j],'Sig2_', sigma2True[i]))) 
#   meanTlmDfTemp<-cbind(meanTlmDfTemp, get(paste0('meanTlm', sumNames[j],'Sig2_', sigma2True[i])))
 }
mseDf<-rbind(mseDf, mseDfTemp)
klDf<-rbind(klDf, klDfTemp)
#trimMspeDf<-rbind(trimMspeDf, trimMspeDfTemp)
#tlmDf<-rbind(tlmDf, tlmDfTemp)

#meanTlmDf<-rbind(meanTlmDf,  meanTlmDfTemp)
}
row.names(mseDf)<-1:nrow(mseDf)
row.names(klDf)<-1:nrow(klDf)
#row.names(trimMspeDf)<-1:nrow(trimMspeDf)
#row.names(tlmDf)<-1:nrow(tlmDf)

#row.names(meanTlmDf)<-1:nrow(meanTlmDf)
mseDf<-as.data.frame(mseDf[-1,])
klDf<-as.data.frame(klDf[-1,])
#trimMspeDf<-as.data.frame(trimMspeDf[-1,])
#tlmDf<-as.data.frame(tlmDf[-1,])
#meanTlmDf<-as.data.frame(meanTlmDf[-1,])
names(mseDf)<-c('p','n','m', 'sig2',sumNames)
names(klDf)<-c('p','n','m', 'sig2',sumNames)
#names(trimMspeDf)<-c('p','n','m', 'sig2',sumNames)
#names(tlmDf)<-c('p','n','m', 'sig2',sumNames)
#names(meanTlmDf)<-c('p','n','m', 'sig2',sumNames)
row.names(mseDf)<-1:nrow(mseDf)
row.names(klDf)<-1:nrow(klDf)
#row.names(trimMspeDf)<-1:nrow(trimMspeDf)
#row.names(tlmDf)<-1:nrow(tlmDf)
#row.names(meanTlmDf)<-1:nrow(meanTlmDf)
head(mseDf)
head(klDf)
#head(trimMspeDf)
#head(tlmDf)
#head(meanTlmDf)
# fit1<-lm(Full~n*p*m, data=mseDf)
# summary(fit1)

round(apply(mseDf,2,mean),5) # SteveMod
round(apply(klDf,2,mean),5) # SteveMod

#
# Plot Main Effects ----
# (by method: i.e. 2-way interactions of p,n,m, and sigma2, with method)

# MSE--squared error (theta_i-theta_ihat)^2 metric

mse_p<-daply(mseDf, .(p), colwise(mean))
#msesd_p<-daply(mseDf, .(p), colwise(sd))
# colMeans(mseDf[mseDf[,'p']==.1,])

mse_n<-daply(mseDf, .(n), colwise(mean))
# colMeans(mseDf[mseDf[,'n']==25,])
mse_m<-daply(mseDf, .(m), colwise(mean))
#
mse_sig2<-daply(mseDf, .(sig2), colwise(mean))



# mse main effect p
#png('mseMainEffects.png', width=800)

if(concTau2Near1){
  pdf('mseMainEffects_concTau2Near1.pdf')
} else {
pdf('mseMainEffects.pdf')
}
matplot(mse_p[,4:8], type='b', pch=19,xaxt='n', xlab='p', ylab='SE', col=1:5)
title('Main effect of p')
axis(side=1, at=1:nrow(mse_p), labels=row.names(mse_p))


# mse main effect n
matplot(mse_n[,4:8], type='b', pch=19,xaxt='n', xlab='n', ylab='SE')
title('Main effect of n')
axis(side=1, at=1:nrow(mse_n), labels=row.names(mse_n))


# mse main effect m
matplot(mse_m[,4:8], type='b', pch=19,xaxt='n', xlab='m', ylab='SE')
title('Main effect of m')
axis(side=1, at=1:nrow(mse_m), labels=row.names(mse_m))


# mse main effect sigma2
matplot(mse_sig2[,4:8], type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='SE')
title(expression(paste('Main effect of ', sigma^2)))
axis(side=1, at=1:nrow(mse_sig2), labels=row.names(mse_sig2))
dev.off()

plot(unlist(mse_sig2[,4]), ylim=c(0,.5), type='b')
lines(unlist(mse_sig2[,5]), col=2)
lines(unlist(mse_sig2[,6]), col=3)
lines(unlist(mse_sig2[,7]), col=4)
lines(unlist(mse_sig2[,8]), col=5)

# #
# # trimed MSPE--squared error E[(y-theta_hat)^2] metric----
# #
# 
# trimMspe_p<-daply(trimMspeDf, .(p), colwise(mean))
# trimMspe_n<-daply(trimMspeDf, .(n), colwise(mean))
# trimMspe_m<-daply(trimMspeDf, .(m), colwise(mean))
# trimMspe_sig2<-daply(trimMspeDf, .(sig2), colwise(mean))
# 
# 
# 
# # mspe main effect p
# #png('trimMspeMainEffects.png', width=800)
# 
# 
# if(concTau2Near1){
#   pdf('trimMspeMainEffects_concTau2Near1.pdf')
# } else {
#   pdf('trimMspeMainEffects.pdf')
# }
# 
# # tmp_p<-matrix(unlist(trimMspe_p[,4:8]), nrow=3, ncol=5)-matrix(colMeans(matrix(unlist(trimMspe_p[,4:8]), nrow=3, ncol=5)), nrow=3,ncol=5, byrow=TRUE)
# 
# matplot(trimMspe_p[,4:8], type='b', pch=19,xaxt='n', xlab='p', ylab='Trimmed MSPE/sigma2', col=1:5)
# # matplot(tmp_p, type='b', pch=19,xaxt='n', xlab='p', ylab='Trimmed MSPE/sigma2', col=1:5)
# title('Main effect of p')
# axis(side=1, at=1:nrow(trimMspe_p), labels=row.names(trimMspe_p))
# 
# 
# # mse main effect n
# matplot(trimMspe_n[,4:8], type='b', pch=19,xaxt='n', xlab='n', ylab='Trimmed MSPE/sigma2')
# title('Main effect of n')
# axis(side=1, at=1:nrow(trimMspe_n), labels=row.names(trimMspe_n))
# 
# 
# # mse main effect m
# matplot(trimMspe_m[,4:8], type='b', pch=19,xaxt='n', xlab='m', ylab='Trimmed MSPE/sigma2')
# title('Main effect of m')
# axis(side=1, at=1:nrow(trimMspe_m), labels=row.names(trimMspe_m))
# 
# 
# # mse main effect sigma2
# matplot(trimMspe_sig2[,4:8], type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='Trimmed MSPE/sigma2')
# title(expression(paste('Main effect of ', sigma^2)))
# axis(side=1, at=1:nrow(trimMspe_sig2), labels=row.names(trimMspe_sig2))
# dev.off()
# 
# 


#
#K-L ----
#

kl_p<-daply(klDf, .(p), colwise(mean))
# colmeans(klDf[klDf[,'p']==.1,])
names(klDf)
kl_n<-daply(klDf[klDf[,'sig2']==1,], .(n), colwise(mean))


kl_n1<-daply(klDf[klDf[,'sig2']==1,], .(n), colwise(mean))
kl_n10<-daply(klDf[klDf[,'sig2']==10,], .(n), colwise(mean))


# colmeans(klDf[klDf[,'n']==25,])
kl_m<-daply(klDf, .(m), colwise(mean))
#
kl_sig2<-daply(klDf, .(sig2), colwise(mean))



# kl main effect p

if(concTau2Near1){
  pdf('KLMainEffects_concTau2Near1.pdf')
  
} else {
  pdf('KLMainEffects.pdf')
}


matplot(kl_p[,4:8], type='b', pch=19,xaxt='n', xlab='p', ylab='K-L')
title('Main effect of p')
axis(side=1, at=1:nrow(kl_p), labels=row.names(kl_p))

#leave full off
matplot(kl_p[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='p', ylab='K-L')
title('Main effect of p')
axis(side=1, at=1:nrow(kl_p), labels=row.names(kl_p))


# 
# # kl main effect n
matplot(kl_n[,4:8], type='b', pch=19,xaxt='n', xlab='n', ylab='K-L')
title('Main effect of n')
axis(side=1, at=1:nrow(kl_n), labels=row.names(kl_n))


# # kl main effect n

pdf('klMainEffectBySig2.pdf')
matplot(kl_n1[,5:8], type='b', col=2:5, pch=19,xaxt='n', xlab='n', ylab='K-L')
title('Main effect of n, sigma2=1')
axis(side=1, at=1:nrow(kl_n1), labels=row.names(kl_n1))

matplot(kl_n10[,5:8], type='b', col=2:5, pch=19,xaxt='n', xlab='n', ylab='K-L')
title('Main effect of n- sigma2=10')
axis(side=1, at=1:nrow(kl_n10), labels=row.names(kl_n10))

dev.off()




#leave full off
matplot(kl_n[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='n', ylab='K-L')
title('Main effect of n')
axis(side=1, at=1:nrow(kl_n), labels=row.names(kl_n))


# kl main effect m
matplot(kl_m[,4:8], type='b', pch=19,xaxt='n', xlab='m', ylab='K-L')
title('Main effect of m')
axis(side=1, at=1:nrow(kl_m), labels=row.names(kl_m))

#leave full off
matplot(kl_m[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='m', ylab='K-L')
title('Main effect of m')
axis(side=1, at=1:nrow(kl_m), labels=row.names(kl_m))




# kl main effect sigma2
matplot(kl_sig2[,4:8], type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='K-L')
title(expression(paste('Main effect of ', sigma^2)))
axis(side=1, at=1:nrow(kl_sig2), labels=row.names(kl_sig2))


# leave off full
matplot(kl_sig2[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='K-L')
title(expression(paste('Main effect of ', sigma^2)))
axis(side=1, at=1:nrow(kl_sig2), labels=row.names(kl_sig2))

dev.off()
# 
# plot(colMeans(klDf)[6:9])


# plot(unlist(kl_sig2[,5]),col=2, ylim=c(0.07,.12), type='b')
# lines(unlist(kl_sig2[,6]), col=3)
# lines(unlist(kl_sig2[,7]), col=4)
# lines(unlist(kl_sig2[,8]), col=5)

# 
# 
# 
# #Investigate interactions with KL
# 
# head(klDf)
# 
# klDf2<-melt(klDf, id.vars = c('p','n','m','sig2'))
# dim(klDf)
# klDf2[720:730,]
# 
# names(klDf2)[c(5,6)]<-c('model', 'kl')
# names(klDf2)
# apply(klDf2, 2,class)
# colNames<-c('p','n','m','sig2')
# klDf2[,colNames] <- lapply(klDf2[,colNames] , factor)
# 
# 
# # Effect of n
# fit1<-lm(kl~.^2, data=klDf2[,c(1,2,5,6)])
# summary(fit1)
# formula(fit1)
# 
# names(klDf2)
# newdf<-expand.grid(p=factor(.3),n=levels(klDf2[,'n']), model=levels(klDf2[,'model']))
# cbind(newdf,predict(fit1,newdata = newdf))
# 
# #
# #
# #
# #
# #
# #
# 
# 
# # plot(unlist(trimMspe_sig2[,4]), ylim=c(0,.5))
# # lines(unlist(trimMspe_sig2[,5]), col=2)
# # lines(unlist(trimMspe_sig2[,6]), col=3)
# # lines(unlist(trimMspe_sig2[,7]), col=4)
# # lines(unlist(trimMspe_sig2[,8]), col=5)






#
#
#
#Tlm----
#
#
# 
# tlm_p<-daply(tlmDf, .(p), colwise(mean))
# # colMeans(tlmDf[tlmDf[,'p']==.1,])
# 
# tlm_n<-daply(tlmDf, .(n), colwise(mean))
# # colMeans(tlmDf[tlmDf[,'n']==25,])
# tlm_m<-daply(tlmDf, .(m), colwise(mean))
# #
# tlm_sig2<-daply(tlmDf, .(sig2), colwise(mean))
# #colMeans(tlmDf[tlmDf[,'sig2']==10,])
# 
# 
# #tlm_sig2_n<-daply(tlmDf, .(sig2,n), colwise(mean))
# 
# 
# 
# if(concTau2Near1){
#   pdf('tlmMainEffect_concTau2Near1.pdf')
# } else {
#   pdf('tlmMainEffect.pdf')
# }
# 
# 
# # tlm main effect p
# matplot(tlm_p[,4:8], type='b', pch=19,xaxt='n', xlab='p', ylab='TLM')
# title('Main effect of p')
# axis(side=1, at=1:nrow(tlm_p), labels=row.names(tlm_p))
# 
# #leave full off
# matplot(tlm_p[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='p', ylab='TLM')
# title('Main effect of p')
# axis(side=1, at=1:nrow(tlm_p), labels=row.names(tlm_p))
# 
# # tlm main effect n
# matplot(tlm_n[,4:8], type='b', pch=19,xaxt='n', xlab='n', ylab='TLM')
# title('Main effect of n')
# axis(side=1, at=1:nrow(tlm_n), labels=row.names(tlm_n))
# 
# #leave full off
# matplot(tlm_n[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='n', ylab='TLM')
# title('Main effect of n')
# axis(side=1, at=1:nrow(tlm_n), labels=row.names(tlm_n))
# 
# 
# # tlm main effect m
# matplot(tlm_m[,4:8], type='b', pch=19,xaxt='n', xlab='m', ylab='mean TLM')
# title('Main effect of m')
# axis(side=1, at=1:nrow(tlm_m), labels=row.names(tlm_m))
# 
# #leave full off
# matplot(tlm_m[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab='m', ylab='mean TLM')
# title('Main effect of m')
# axis(side=1, at=1:nrow(tlm_m), labels=row.names(tlm_m))
# 
# 
# 
# 
# # tlm main effect sigma2
# matplot(tlm_sig2[,4:8], type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='mean TLM')
# title(expression(paste('Main effect of ', sigma^2)))
# axis(side=1, at=1:nrow(tlm_sig2), labels=row.names(tlm_sig2))
# 
# 
# # leave off full
# matplot(tlm_sig2[,5:8],col=2:5, lty=2:5, type='b', pch=19,xaxt='n', xlab=expression(sigma^2), ylab='K-L')
# title(expression(paste('Main effect of ', sigma^2)))
# axis(side=1, at=1:nrow(tlm_sig2), labels=row.names(tlm_sig2))
# 
# dev.off()
# 
