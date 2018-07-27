
#Function to compute KL-divergence.
# prediction function is normal distribution
compute_pred_dist_normal<-function(ygrid, thetaSamples, sigma2Samples){
  #ygrid: grid of y values to evaluate the pred distribution
  #thetaSamples, sigma2Samples MCMC samples from the given group
  if(length(ygrid)>1000){
    sapply(ygrid, FUN=function(x) {
      mean(dnorm(x, thetaSamples, sqrt(sigma2Samples)))
    } 
    ) } else {
      ygridN<-matrix(ygrid, nrow=length(thetaSamples), ncol=length(ygrid), byrow = TRUE)
      colMeans(dnorm(ygridN, thetaSamples, sqrt(sigma2Samples)))
    }
}

#Function to compute KL-divergence.
# prediction function is normal distribution
compute_pred_dist_t<-function(ygrid, thetaSamples, sigma2Samples, nuSamples){
  #ygrid: grid of y values to evaluate the pred distribution
  #thetaSamples, sigma2Samples MCMC samples from the given group
  if(length(ygrid)>1000){
    sapply(ygrid, FUN=function(x) {
    mean(dt((x-thetaSamples)/sqrt(sigma2Samples), df = nuSamples)/sqrt(sigma2Samples))
    } 
    )  } else {
     ygridN<-matrix(ygrid, nrow=length(thetaSamples), ncol=length(ygrid), byrow = TRUE)
    y_grid_norm <- apply(ygridN, 2, function(x)  (x - thetaSamples)/sqrt(sigma2Samples))
     t_pdf <- apply(y_grid_norm, 2, function(x) dt(x, df = nuSamples))
    colMeans(apply(t_pdf, 2, function(x) x/sqrt(sigma2Samples)))
    }
}




# function to compute the K-L metric, for one group

compute_KL<-function(thetaSamples, sigma2Samples, nuSamples = NULL, ngridpts, prediction_dist = 'Normal', thetaTrue = 0, sig2True = 1){
  #thetaSamples: MCMC sample of the theta_i from the given group
  #sigma2Sample: MCMC sample of the sigma2_i from the given group
  #nuSamples: MCMC sample of df paramater if t model samples
  #thetaTrue: scalar of the true value of theta
  #sig2True: scalar of the true value of sigma2
  #prediction_dist: 'Normal' or 't' determines call to compute_pred_dist_normal or compute_pred_dist_t
  #ngridpts: number of grid points
  
  ygrid<-seq(thetaTrue-10*sqrt(sig2True),thetaTrue+10*sqrt(sig2True), length.out=ngridpts)
  predDist<-if(prediction_dist == 'Normal') {compute_pred_dist_normal(ygrid, thetaSamples, sigma2Samples)
  } else if(prediction_dist == 't'){
    compute_pred_dist_t(ygrid, thetaSamples, sigma2Samples, nuSamples)
  } else{
      stop('prediction_dist must be Normal or t')
    }
  trueDist<-dnorm(ygrid, thetaTrue, sqrt(sig2True))
  dif<-diff(ygrid)
  kl<-((log(trueDist)-log(predDist))*trueDist)[-1]
  KL<-sum((kl*dif)[kl!=Inf])
  KL
}
