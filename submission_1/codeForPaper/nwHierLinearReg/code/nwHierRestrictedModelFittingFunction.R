 ##################
#Restricted posterior Code
# condition on individual statistics from each group
#base model is in HierNormTheoryLm.R


#T() is a robust summary of the data; Here T(X) represents the regression and scale estimates from each group individually



source("nwHierNormalModelSamplingFunctions.R")
source("nwHierNormalModelFittingFunction.R")
source("../../linearRegressionHuberAndProposal2.R")



####################################################
#new: a bit faster; uses fn.one.rep.y2 which uses fn.attenuation2-a slightly faster version of fn.attenuation
hierNormTheoryRestLm2<-function(y,
                                X,
                                regEst='Huber',
                                scaleEst='Huber',
                                nkeep=1e4, nburn=1e3,
                                mu0,
                                Sigma0,
                                a0, 
                                b0,
                                mu_bstr,
                                psi_bstr,
                                swSq=1,
                                w1,
                                w2, 
                                a_psir,
                                b_psir,
                                maxit=400)
{
  #y is a list of the responses for each group
  #X is the design Matrix-a list of the design matrices Xi for each group
  #X[[i]] is the design matrix for each group
  
  mu0<<-mu0
  Sigma0<<-Sigma0
  a0<<-a0
  b0<<-b0
  #alpha_mustr<<-alpha_mustr
  #beta_mustr<<-beta_mustr
  #a_psib<<-a_psib
  #b_psib<<-b_psib
  swSq<<-swSq
  w1<<-w1; w2<<-w2
  a_psir<<-a_psir
  b_psir<<-b_psir
  v12<<-fn.compute.ab(mu_bstr,psi_bstr)
  v1<<-v12[1]
  v2<<-v12[2]
  ############
  projList=NULL #leave this as null; projection onto deviation space for each group
  XtX<-lapply(X, FUN=function(X) t(X)%*%X)
  p<<-length(mu0) #the number of reg. coefficients per group
  pTot<<-length(X)*p
  ni<<-sapply(y, length)
  nGroups<<-pTot/p
  Sigma0Inv<<-solve(Sigma0)
  total<-nkeep+nburn
  #initialize outputs
  #list of betaSamples: each element of list is betaSamples from corresponding group
  betaGroupSamples<-array(NA, c(p,nGroups, total))
  BetaSamples<-matrix(NA,total,p)
  sigma2Samples<-matrix(NA,total,nGroups)
  #mu_bstrSamples<-numeric(total)
  #psi_bstrSamples<-numeric(total)
  bstarSamples<-numeric(total)
  mu_rhoSamples<-numeric(total)
  psi_rhoSamples<-numeric(total)
  rhoSamples<-numeric(total)
  yAccept<-matrix(NA,total,nGroups)
  logprop.curr<-matrix(NA,total,nGroups)
  yMhRatios<-matrix(NA,total,nGroups)
  
  #initial values
  #starting values of parameters
  #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
  #psi_bstr<-rgamma(1, a_psib,b_psib)  
  #v<-fn.compute.ab(mu_bstr,psi_bstr)
  bstar<-mu_bstr#rbeta(1, v[1],v[2])
  Beta<-mvrnorm(1,mu=mu0, Sigma=(1-bstar)*Sigma0)
  betalMat<-matrix(Beta, p,nGroups)
  mu_rho<-rbeta(1,w1,w2)
  psi_rho<-rgamma(1, a_psir,b_psir)
  ab<-fn.compute.ab(mu_rho,psi_rho)
  #rho<-rbeta(1,ab[1],ab[2])
  rho<-runif(1, .1, .9)
  Sigma_rho<-(1-rho)*diag(nGroups)+matrix(1, nGroups,nGroups)*rho
  Z<-mvrnorm(1, rep(0,nGroups), Sigma_rho)
  
  
  
  
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
  
  if(scaleEst!='Huber'){
    stop('scale estimate only set up for Hubers Prop2')
  }
  fn.chi<-fn.chi.prop2
  
  robustList<-list();length(robustList)<-nGroups
  bHatObsList<-list();length(bHatObsList)<-nGroups
  sigHatObsList<-list();length(sigHatObsList)<-nGroups
  
  #fit the robust linear model to each group separately
  #modified for small groups?
  for(groups in 1:nGroups){
    robust<-rlm(X[[groups]],y[[groups]],psi=psi, scale.est=scaleEst, maxit=maxit)
    robustList[[groups]]<-robust
    bHatObsList[[groups]]<-robust$coef
    sigHatObsList[[groups]]<-robust$s
  }
  
  #####################
  #choose starting value for yi
  #####################
  log.prop.den.curr<-numeric(nGroups)
  yCurr<-list(); length(yCurr)<-nGroups
  for(i in 1:nGroups){
    yCurr[[i]]<-numeric(ni[i])
  }
  if(is.null(projList)) {projList<-list(); length(projList)<-nGroups
                         QtList<-list(); length(QtList)<-nGroups
                         for(i in 1:nGroups){
                           Q<-qr.Q(qr(X[[i]]))
                           projList[[i]]<-diag(ni[i])-tcrossprod(Q,Q)
                            QtList[[i]]<-t(Q)
                         }
  }
  #finding starting yi's for each group 
  for(i in 1:nGroups){
    y.prop <- rnorm(ni[i])
    Xi<-X[[i]]
    bHatObsi<-bHatObsList[[i]] #observed stats in each group
    sigHatObsi<-sigHatObsList[[i]]
    proji<-projList[[i]]
    Qti<-QtList[[i]]
    yCurri <- fn.comp.ystst(y.prop,Xi,l1obs=bHatObsi,s1obs=sigHatObsi,psi,scale.est=scaleEst,maxit)
    yCurri<-yCurri
    yCurr[[i]]<-yCurri
    log.prop.den.curr[i]<-log.prop.den2(yCurri,Xi, proji,Qti,bHatObsi, sigHatObsi,fn.psi, fn.chi,ni[i],p)
  }
  
  
  ################## 
  #MCMC 
  #################
  for(iter in 1:total){
    
    #[theta|everything] step
    samp<-fn.hier.one.rep(yCurr,
                          X,
                          XtX,
                          v1,#mu_bstr,
                          v2,#psi_bstr,
                          bstar,
                          Beta,
                          betalMat, # first to update
                          Z,
                          mu_rho,
                          psi_rho,
                          rho)
    #update temp values
    Beta<-samp$Beta
    betalMat<-samp$betalMat
    Z<-samp$Z
    sigma2<-samp$sigma2
    #mu_bstr<-samp$mu_bstr
    #psi_bstr<-samp$psi_bstr
    bstar<-samp$bstar
    mu_rho<-samp$mu_rho
    psi_rho<-samp$psi_rho
    rho<-samp$rho
    #update outputs
    BetaSamples[iter,]<-Beta
    betaGroupSamples[,,iter]<-betalMat
    sigma2Samples[iter,]<-sigma2
   # mu_bstrSamples[iter]<-mu_bstr
    #psi_bstrSamples[iter]<-psi_bstr
    bstarSamples[iter]<-bstar
    mu_rhoSamples[iter]<-mu_rho
    psi_rhoSamples[iter]<-psi_rho
    rhoSamples[iter]<-rho
    
    #[y|everything + robust statistics] step; each group updated separately
    for(gp in 1:nGroups){
      yicurr<-yCurr[[gp]]
      Xi<-X[[gp]]
      betaCuri<-betalMat[,gp]
      sigma2Curi<-sigma2[gp]
      bHatObsi<-bHatObsList[[gp]]
      sigHatObsi<-sigHatObsList[[gp]]
      log.prop.den.curri<-log.prop.den.curr[gp]
      proj<-projList[[gp]]
      Qt<-QtList[[gp]]
      ySample<-fn.one.rep.y2(yicurr,betaCuri,sqrt(sigma2Curi),bHatObsi, sigHatObsi,Xi, log.prop.den.curri, proj,Qt,fn.psi,fn.chi, psi,scaleEst,maxit)
      yCurr[[gp]]<-ySample[[1]]
      yAccept[iter,gp]<-ySample[[2]]
      log.prop.den.curr[gp]<-ySample[[3]]
      logprop.curr[iter,gp]<-ySample[[3]]
      yMhRatios[iter,gp]<-ySample[[4]]
    }
  }
  
  out<-list()
  out$Beta<-BetaSamples[-c(1:nburn),]
  out$betal<-betaGroupSamples[,,-c(1:nburn)]
  out$sigma2s<-sigma2Samples[-c(1:nburn),]
 # out$mu_bstr<-mu_bstrSamples[-c(1:nburn)]
#  out$psi_bstr<-psi_bstrSamples[-c(1:nburn)]
  out$bstar<-bstarSamples[-c(1:nburn)]
  out$mu_rho<-mu_rhoSamples[-c(1:nburn)]
  out$psi_rho<-psi_rhoSamples[-c(1:nburn)]
  out$rho<-rhoSamples[-c(1:nburn)]
#   out$yAccept<-yAccept #colMeans(yAccept)
  out$yAccept<-colMeans(yAccept)
  out$rlmFits<-robustList
  hypers<-c(a0,b0,swSq,w1,w2,a_psir,b_psir, mu_bstr,psi_bstr)
  names(hypers)<-c("a0","b0","swSq",'w1',"w2","a_psir","b_psir", 'mu_bstr', "psi_bstr")
  out$hypers<-hypers
  out
}       
       


# hierNormTheoryRestLm<-function(y,
#                                X,
#                                regEst='Huber',
#                                scaleEst='Huber',
#                                nkeep=1e4, nburn=1e3,
#                                mu0,
#                                Sigma0,
#                                a0, 
#                                b0,
#                                #alpha_mustr,
#                                #beta_mustr,
#                                # a_psib,
#                                # b_psib,
#                                mu_bstr,
#                                psi_bstr,
#                                swSq=1,
#                                w1,
#                                w2, 
#                                a_psir,
#                                b_psir,
#                                maxit=400)
# {
#   #y is a list of the responses for each group
#   #X is the design Matrix-a list of the design matrices Xi for each group
#   #X[[i]] is the design matrix for each group
#   
#   mu0<<-mu0
#   Sigma0<<-Sigma0
#   a0<<-a0
#   b0<<-b0
#   alpha_mustr<<-alpha_mustr
#   beta_mustr<<-beta_mustr
#   a_psib<<-a_psib
#   b_psib<<-b_psib
#   swSq<<-swSq
#   w1<<-w1; w2<<-w2
#   a_psir<<-a_psir
#   b_psir<<-b_psir
#   ############
#   projList=NULL #leave this as null; projection onto deviation space for each group
#   XtX<-lapply(X, FUN=function(X) t(X)%*%X)
#   p<<-length(mu0) #the number of reg. coefficients per group
#   pTot<<-length(X)*p
#   ni<<-sapply(y, length)
#   nGroups<<-pTot/p
#   Sigma0Inv<<-solve(Sigma0)
#   total<-nkeep+nburn
#   #initialize outputs
#   #list of betaSamples: each element of list is betaSamples from corresponding group
#   betaGroupSamples<-array(NA, c(p,nGroups, total))
#   BetaSamples<-matrix(NA,total,p)
#   sigma2Samples<-matrix(NA,total,nGroups)
#   #mu_bstrSamples<-numeric(total)
#   #psi_bstrSamples<-numeric(total)
#   bstarSamples<-numeric(total)
#   mu_rhoSamples<-numeric(total)
#   psi_rhoSamples<-numeric(total)
#   rhoSamples<-numeric(total)
#   yAccept<-matrix(NA,total,nGroups)
#   logprop.curr<-matrix(NA,total,nGroups)
#   yMhRatios<-matrix(NA,total,nGroups)
#   
#   #initial values
#   #starting values of parameters
#   #mu_bstr<-rbeta(1, alpha_mustr,beta_mustr)
#   #psi_bstr<-rgamma(1, a_psib,b_psib)  
#   v<-fn.compute.ab(mu_bstr,psi_bstr)
#   bstar<-mu_bstr#rbeta(1, v[1],v[2])
#   Beta<-mvrnorm(1,mu=mu0, Sigma=(1-bstar)*Sigma0)
#   betalMat<-matrix(Beta, p,nGroups)
#   mu_rho<-rbeta(1,w1,w2)
#   psi_rho<-rgamma(1, a_psir,b_psir)
#   ab<-fn.compute.ab(mu_rho,psi_rho)
#   #rho<-rbeta(1,ab[1],ab[2])
#   rho<-runif(1)
#   Sigma_rho<-(1-rho)*diag(nGroups)+matrix(1, nGroups,nGroups)*rho
#   Z<-mvrnorm(1, rep(0,nGroups), Sigma_rho)
#   
#   
#   
#   
#   ############################
#   #define the psi function
#   ############################
#   if(regEst=='Huber') {
#     psi<-get('psi.huber') #internal
#     fn.psi<-get('fn.psi.huber')
#     
#   } else { 
#     if(regEst=='Tukey'){
#       psi<-get('psi.bisquare') #internal
#       fn.psi<-get('fn.psi.bisquare')
#     } else {stop("only set up for Huber or Tukey regression estimates")}}
#   
#   if(scaleEst!='Huber'){
#     stop('scale estimate only set up for Hubers Prop2')
#   }
#   fn.chi<-fn.chi.prop2
#   
#   robustList<-list();length(robustList)<-nGroups
#   bHatObsList<-list();length(bHatObsList)<-nGroups
#   sigHatObsList<-list();length(sigHatObsList)<-nGroups
#   
#   #fit the robust linear model to each group separately
#   #modified for small groups?
#   for(groups in 1:nGroups){
#     robust<-rlm(X[[groups]],y[[groups]],psi=psi, scale.est=scaleEst, maxit=maxit)
#     robustList[[groups]]<-robust
#     bHatObsList[[groups]]<-robust$coef
#     sigHatObsList[[groups]]<-robust$s
#   }
#   
#   #####################
#   #choose starting value for yi
#   #####################
#   log.prop.den.curr<-numeric(nGroups)
#   yCurr<-list(); length(yCurr)<-nGroups
#   for(i in 1:nGroups){
#     yCurr[[i]]<-numeric(ni[i])
#   }
#   if(is.null(projList)) {projList<-list(); length(projList)<-nGroups
#                         
#                          for(i in 1:nGroups){
#                            Q<-qr.Q(qr(X[[i]]))
#                            projList[[i]]<-diag(ni[i])-tcrossprod(Q,Q)
#                          
#                          }
#   }
#   #finding starting yi's for each group 
#   for(i in 1:nGroups){
#     y.prop <- rnorm(ni[i])
#    Xi<-X[[i]]
#    bHatObsi<-bHatObsList[[i]] #observed stats in each group
#    sigHatObsi<-sigHatObsList[[i]]
#    proji<-projList[[i]]
#   yCurri <- fn.comp.ystst(y.prop,Xi,l1obs=bHatObsi,s1obs=sigHatObsi,psi,scale.est=scaleEst,maxit)
#    yCurri<-yCurri
#     yCurr[[i]]<-yCurri
#     log.prop.den.curr[i] <-log.prop.den(yCurri,Xi, proji,bHatObsi, sigHatObsi,fn.psi, fn.chi,ni[i],p)
#   }
#   
#   
#   ################## 
#   #MCMC 
#   #################
#     for(iter in 1:total){
#       
#       #[theta|everything] step
#       samp<-fn.hier.one.rep(yCurr,
#                             X,
#                             XtX,
#                             mu_bstr,
#                             psi_bstr,
#                             bstar,
#                             Beta,
#                             betalMat, # first to update
#                             Z,
#                             mu_rho,
#                             psi_rho,
#                             rho)
#       #update temp values
#       Beta<-samp$Beta
#       betalMat<-samp$betalMat
#       Z<-samp$Z
#       sigma2<-samp$sigma2
#       #sigma2<-unlist(sigHatObsList)^2
#       #mu_bstr<-samp$mu_bstr
#       #psi_bstr<-samp$psi_bstr
#       bstar<-samp$bstar
#       mu_rho<-samp$mu_rho
#       psi_rho<-samp$psi_rho
#       rho<-samp$rho
#       #update outputs
#       BetaSamples[iter,]<-Beta
#       betaGroupSamples[,,iter]<-betalMat
#       sigma2Samples[iter,]<-sigma2
#       #mu_bstrSamples[iter]<-mu_bstr
#       #psi_bstrSamples[iter]<-psi_bstr
#       bstarSamples[iter]<-bstar
#       mu_rhoSamples[iter]<-mu_rho
#       psi_rhoSamples[iter]<-psi_rho
#       rhoSamples[iter]<-rho
#       
#       #[y|everything + robust statistics] step; each group updated separately
#         for(gp in 1:nGroups){
#           yicurr<-yCurr[[gp]]
#           Xi<-X[[gp]]
#           betaCuri<-betalMat[,gp]
#           sigma2Curi<-sigma2[gp]
#           bHatObsi<-bHatObsList[[gp]]
#           sigHatObsi<-sigHatObsList[[gp]]
#           log.prop.den.curri<-log.prop.den.curr[gp]
#           proj<-projList[[gp]]
#           #Qt<-QtList[[gp]]
#           ySample<-fn.one.rep.y(yicurr,betaCuri,sqrt(sigma2Curi),bHatObsi, sigHatObsi,Xi, log.prop.den.curri, proj,fn.psi,fn.chi, psi,scaleEst,maxit)
#           yCurr[[gp]]<-ySample[[1]]
#           yAccept[iter,gp]<-ySample[[2]]
#           log.prop.den.curr[gp]<-ySample[[3]]
#           logprop.curr[iter,gp]<-ySample[[3]]
#           yMhRatios[iter,gp]<-ySample[[4]]
# }
#     }
#     
#     out<-list()
#     out$Beta<-BetaSamples[-c(1:nburn),]
#     out$betal<-betaGroupSamples[,,-c(1:nburn)]
#     out$sigma2s<-sigma2Samples[-c(1:nburn),]
#     #out$mu_bstr<-mu_bstrSamples[-c(1:nburn)]
#     #out$psi_bstr<-psi_bstrSamples[-c(1:nburn)]
#     out$bstar<-bstarSamples[-c(1:nburn)]
#     out$mu_rho<-mu_rhoSamples[-c(1:nburn)]
#     out$psi_rho<-psi_rhoSamples[-c(1:nburn)]
#     out$rho<-rhoSamples[-c(1:nburn)]
#     out$yAccept<-colMeans(yAccept)
#     out$rlmFits<-robustList
#     hypers<-c(a0,b0,alpha_mustr,beta_mustr,a_psib,b_psib,swSq,w1,w2,a_psir,b_psir, mu_bstr,psi_bstr)
#     names(hypers)<-c("a0","b0","alpha_mustr","beta_mustr",'a_psib','b_psib',"swSq",'w1',"w2","a_psir","b_psir", 'mu_bstr', "psi_bstr")
#     out$hypers<-hypers
#     out
#   } 
