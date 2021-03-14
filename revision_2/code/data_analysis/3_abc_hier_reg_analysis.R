# ABC Hierarchical regression - i.e., grouped by states. 

# download brlm package if needed
# library(devtools)
# install_github('jrlewi/brlm')
library(MCMCpack)
library(MASS)
library(parallel)
lib.loc = file.path("/Users", "john", "Dropbox",
                    "school", "osu", "dissertationResearch",
                    "snm", "rPackage", "brlm", "packrat", "lib",
                    "x86_64-apple-darwin15.6.0", "3.4.3")
library(brlm, lib.loc = lib.loc)
brlm:::update_group_abc
brlm:::proposal_group_abc
library(tidyverse)

nburn <- 20000
nkeep <- 20000
maxit <- 1000 #parameter in MASS::rlm

iter_check <- 200
bw_mult <- 1.2
min_accept_rate <- 0.1

#starting bandwidth for abc kernel for each group.
bandwidth <- .01 #sapply(models_lm, function(fit) sum(sqrt(diag(vcov(fit)))))

# aux funcitons ----

fn.compute.marginals.hierModelNormal<-function(betalsamples, sigma2lsamples, yhold,Xhold){
  
  
  fits<-function(betahats, X){X%*%betahats}
  
  betalSampsList <- lapply(1:dim(betalsamples)[3],function(x) betalsamples[,,x])
  names(betalSampsList) <- names(yhold)
  #holdout means
  
  nMuList<-mapply(fits, betalSampsList, Xhold)
  
  #sd's across samples for each holdout set
  nSigmaList<-split(sqrt(sigma2lsamples), rep(1:ncol(sigma2lsamples), each = nrow(sigma2lsamples)))
  names(nSigmaList)<-names(yhold)
  
  mapply(FUN=function(yh,mean,sd){
    pdfs <- dnorm(yh,mean,matrix(rep(sd, length(yh)), length(yh), length(sd),byrow=TRUE))
    mns <- rowMeans(pdfs)
    sds <- apply(pdfs, 1, sd)
    cbind(mns, sds)
  }
  ,yhold,nMuList,nSigmaList)
  
}


# load data and prior information ----
analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))

analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>%
  group_by(State) %>% 
  filter(n() >= 25) %>% ungroup() %>% 
  mutate(State = factor(State)) %>% #filter(!State %in% c(14, 30)) %>% 
  arrange(State)

analysis_data <- analysis_data %>% 
  filter(!State %in% c('14', '30')) %>% 
  mutate(State = droplevels(State))

# analysis_data %>%  group_by(State) %>% 
#   summarize(sd1 = sd(Associate_Count), u1 = length(unique(Associate_Count)),  sd2 = sd(Office_Employee_Count),
#             u2 = length(unique(Office_Employee_Count)),
#             n())

# center and scale X variables - to make scale of
# coefficients consistent for the ABC distance measure.
# columns_to_scale <- c("sqrt_count_2010",
#                        "Associate_Count", 
#                        "Office_Employee_Count")
# 
# analysis_data[, columns_to_scale] <- 
#   scale(analysis_data[, columns_to_scale])

# apply(analysis_data[, columns_to_scale], 2, mean)
# apply(analysis_data[, columns_to_scale], 2, sd)


#defined globally.
nGroups <<- length(unique(analysis_data$State))

state_sizes <- analysis_data %>% 
  dplyr::group_by(State) %>% 
  dplyr::summarise(n = n())

#Set prior parameters ----
parms_prior <- read_rds(file.path(here::here(), 'hier_parms_prior.rds'))
mu0 <<- parms_prior$mu_0
Sigma0 <<- parms_prior$Sigma_0
a0 <<- parms_prior$a_0
b0 <<- parms_prior$b_0
mu_bstr <<- parms_prior$mu_bstr
psi_bstr <<- parms_prior$psi_bstr
w1 <<- parms_prior$w1
w2 <<- parms_prior$w2
a_psir <<- parms_prior$a_psir
b_psir <<- parms_prior$b_psir
parms_prior$trend

trend <- sqrt_count_2012 ~ sqrt_count_2010  + -1 +
  Associate_Count +
  Office_Employee_Count


set.seed(123)
N <- nrow(analysis_data)
p <<- length(mu0)


#get the training/holdout indices used in the original analysis in 3_hier_reg_analysis
num_batch <- 5
rds_path <- file.path(here::here(), "hier_reg_n1547_sim_number_")
tmp1 <- readRDS(file.path(paste0(rds_path, 1, '.rds')))
holdIndicesMatrix <- tmp1$holdIndices
#combine results into one
for(i in 2:num_batch){
  tmp <-  readRDS(file.path(paste0(rds_path, i, '.rds')))
  holdIndicesMatrix <- rbind(holdIndicesMatrix, tmp$holdIndices)
}

trainIndicesMatrix <- t(apply(holdIndicesMatrix, 1, function(x) (1:N)[-x]))
dim(trainIndicesMatrix)
intersect(trainIndicesMatrix[1,], holdIndicesMatrix[1,])
n <- N - ncol(trainIndicesMatrix)

# estimates ----

#posterior means.
#order is beta_ls, sigma
reps <- nrow(trainIndicesMatrix) #number of training sets
group_estimates <- group_estimates_sds <- array(NA, c(p + 1,nGroups,reps)) 
#overall Beta
estimates <- estimates_sd <- array(NA, c(p, reps))
#marginals  and predictions of each y in holdout set  and each model ------
marginals <- predictions  <- array(NA, c(reps, N - n))
marginals_sd <- array(NA, c(reps, N - n))
# M-H acceptance rates -----
acceptY <- array(NA, c(nGroups, reps)) #acceptance rates in abc models

#auxilary functions and constants ---- 
#for preds on holdout set. 
fits <- function(betahats, X){X %*% betahats}


# simulation -----
strt <- Sys.time()

#set this up to follow same naming convention as in 3_hier_reg_analysis.R
# for the most part with some added info for abc specific anslysis. 


output_to_save <- function(abc_fit, yhold, Xhold, hold){
  
  out <- list()
  out$nkeep <- nrow(abc_fit$sigma2s)
  out$y_hold <- yhold
  out$y_open <- hold$Count_2012 > 0 
  out$y_type1 <- hold$Type == '1'
  
  betal <- aperm(abc_fit$betal, c(1,3,2)) #format expected for marginals computation
  
  
  mcmc_format <- 
    lapply(seq(dim(betal)[3]), function(x) betal[ , , x]) 
  mcmc_format <- do.call(rbind, mcmc_format) %>% t() %>% 
    mcmc()
  
  out$group_converge_betas <- 
    abs(geweke.diag(mcmc_format)$z)
  
  
  
  postMeansBetal <- apply(betal,c(1,3) , mean)
  #out$postMeansBetal <- postMeansBetal
  postSDsBetal<-apply(betal,c(1,3) , sd)
  #out$postSDsBetal <-postSDsBetal
  
  #sigma2s
  postMeansSigma2s  <- colMeans(abc_fit$sigma2s)
  #out$postMeansSigma2s <- postMeansSigma2s 
  
  
  postSDsSigma2s <- apply(abc_fit$sigma2s,2,sd)
  #out$postSDsSigma2s <- postSDsSigma2s
  
  out$group_estimates <- rbind(postMeansBetal, postMeansSigma2s)
  out$group_estimates_sds <- rbind(postSDsBetal, postSDsSigma2s)
  
  out$group_converge_sigma2 <- 
    abs(geweke.diag(mcmc(abc_fit$sigma2s))$z)
  
  #Beta
  postMeansBETA <- colMeans(abc_fit$Beta)
  out$estimates <- postMeansBETA
  postSDsBeta <- apply(abc_fit$Beta,2 , sd)
  out$estimates_sd <- postSDsBeta
  
  out$converge_beta <- abs(geweke.diag(mcmc(abc_fit$Beta))$z)
  
  out$converge_hyper <- c(
    #bstar converge
    abs(geweke.diag(mcmc(abc_fit$bstar))$z),
    #mu_rho converge
    abs(geweke.diag(mcmc(abc_fit$mu_rho))$z),
    #psi_rho_converge 
    abs(geweke.diag(mcmc(abc_fit$psi_rho))$z),
    #rho_converge
    abs(geweke.diag(mcmc(abc_fit$rho))$z))
  
  #Acceptance rates for new y's
  out$acceptY <- abc_fit$yAccept
  
  #predictionss on holdout set
  postMeansBetalList <- split(postMeansBetal, 
                              rep(1:ncol(postMeansBetal), 
                                  each = nrow(postMeansBetal)))
  restrictedPreds <- mapply(fits, postMeansBetalList, Xhold)
  
  predictions <- restrictedPreds %>% unlist()
  out$predictions <- predictions
  #computing marginal likelihoods for each element in holdout sample
  rest_marg_mn_sd <- 
    fn.compute.marginals.hierModelNormal(betal, abc_fit$sigma2s, yhold,Xhold)
  rest_marg <- lapply(rest_marg_mn_sd, function(x) x[,1])
  rest_marg <- rest_marg %>% unlist() 
  out$marginals <- rest_marg
  rest_marg_sd <- lapply(rest_marg_mn_sd, function(x) x[,2])
  rest_marg_sd <- rest_marg_sd %>% unlist() 
  out$marginals_sd <- rest_marg_sd
  out$bandwidth <- abc_fit$bandwidth
  out
  
  
} 
i <- 1
one_abc_fit <- function(i, holdIndicesMatrix, 
                         trainIndicesMatrix, 
                         analysis_data){
  holdIndices <- holdIndicesMatrix[i, ]
  trainIndices <- trainIndicesMatrix[i,]
  train <- analysis_data[trainIndices,]
  hold <- analysis_data[holdIndices,]
  
  # rlm on training -----      
  by_state <- train %>% 
    group_by(State) %>% 
    nest() 
  nis <- apply(by_state,1, function(dd) dd$data %>% nrow())
  state_lm <- function(df){
    MASS::rlm(trend, y.ret = TRUE, x.ret = TRUE, data = df, maxit = maxit)
  }
  
  models_lm <- by_state$data %>% 
    map(state_lm) 
  names(models_lm) <- by_state$State
  
  #used for tuning params below
  sigHats <- models_lm %>% 
    map(.f = function(m) m$s) %>% 
    unlist()
  
  # convert response and design matrix for training set to lists for brlm functions
  y <- models_lm %>% 
    map(.f = function(m) m$y)
  X <- models_lm %>% 
    map(.f = function(m) m$x)
  

  #prepare the holdout data for predictions
  #yhold is list of holdout responses from each group
  #Xhold is list of holfout design matrices for each group
  by_state_hold <- hold %>% 
    group_by(State) %>% 
    nest()
  
  models_lm_hold <- by_state_hold$data %>% 
    map(state_lm) 
  names(models_lm_hold) <- by_state_hold$State

  
  # convert response and design matrix for holdout set
  yhold <- models_lm_hold %>% 
    map(.f = function(m) m$y)
  
  Xhold <- models_lm_hold %>% 
    map(.f = function(m) m$x)

  ################################################
  # hier abc
  ################################################
  swSq <<- 1      
  
  #tunning parameters for MH step on bstar, mu_rho, psi_rho, and rho
  step_logbstar <- abs(log(mu_bstr/(sqrt(mu_bstr*(1 - mu_bstr)/(psi_bstr + 1))))) #abs log(mean/sd)
  mu_rho_step <- .3
  psi_rho_step <- (a_psir/b_psir^2)^.5 #mean/sd
  rho_step <- .1
  
  # nis <- unlist(lapply(y, length), use.names=FALSE)
  sigs2 <-  sigHats^2
  step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
  if (any(is.na(step_Z))){
    sigs2 <-  1
    step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))
  }  
  
  ################################################
  # ABC Versions -----
  ################################################
  #Tukey version ---
  
  abc_fit <- brlm::hierNormTheoryRestLm(y,
                                           X,
                                           regEst = 'Tukey',
                                           scaleEst = 'Huber',
                                           nkeep, 
                                           nburn,
                                           mu0,
                                           Sigma0,
                                           a0, 
                                           b0,
                                           mu_bstr,
                                           psi_bstr,
                                           swSq = 1,
                                           w1,
                                           w2, 
                                           a_psir,
                                           b_psir,
                                           maxit=maxit,
                                           step_logbstar, 
                                           mu_rho_step, 
                                           psi_rho_step, 
                                           rho_step,
                                           step_Z,
                                           abc = TRUE,
                                           bandwidth = bandwidth,
                                           iter_check = iter_check,
                                           min_accept_rate = min_accept_rate,
                                           bw_mult = bw_mult)
 out <-  output_to_save(abc_fit, yhold, Xhold, hold)
 out$holdIndices <- holdIndices
 out
  }



#cl <- parallel::makeCluster(2)
outs <- parallel::mclapply(X = seq(nrow(holdIndicesMatrix)), FUN = one_abc_fit, 
                            holdIndicesMatrix = holdIndicesMatrix,
                            trainIndicesMatrix = trainIndicesMatrix,
                            analysis_data = analysis_data,
                           mc.cores = parallel::detectCores())

print(
  lapply(outs, function(x) x$acceptY)
)

print(
  lapply(outs, function(x) x$bandwidth)
)


print(
  sapply(outs, function(x) x$nkeep)
)

write_rds(outs, file.path(here::here(), 
                          paste0('hier_abc_reg_n', n, "_n_keep", nkeep, '.rds' )))

Sys.time() - strt

