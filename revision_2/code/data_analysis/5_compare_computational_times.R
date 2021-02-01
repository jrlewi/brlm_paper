# simple comparision of computational timing differences
# between restricted and abc methods

#load locally
lib.loc = "/Users/john/Dropbox/school/osu/dissertationResearch/snm/rPackage/brlm/packrat/lib/x86_64-apple-darwin15.6.0/3.4.3"
library(brlm, lib.loc = lib.loc)
library(tidyverse)
#library(MASS)
# parameters
#Set prior parameters ----
parms_prior <- read_rds(file.path(here::here(), 'hier_parms_prior.rds'))
mu0 <- parms_prior$mu_0
Sigma0 <- parms_prior$Sigma_0
a0 <- parms_prior$a_0
b0 <- parms_prior$b_0
mu_bstr <- parms_prior$mu_bstr
psi_bstr <- parms_prior$psi_bstr
w1 <- parms_prior$w1
w2 <- parms_prior$w2
a_psir <- parms_prior$a_psir
b_psir <- parms_prior$b_psir
swSq <- 1


#additional for abc
iter_check <- 200
bw_mult <- 1.2
min_accept_rate <- 0.1
bandwidth <- 1


#sampling parameters
nburn <- 0 #set length of mcmc chains
nkeep <- 1000
maxit <- 1000




#Data
# load data and prior information ----
analysis_data <- read_rds(file.path(here::here(), 'data', 'analysis_data.rds'))



# parms_prior <- read_rds(file.path(here::here(), 'parms_prior.rds'))
analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), sqrt_count_2012 = sqrt(Count_2012)) %>%
  group_by(State) %>% 
  filter(n() >= 25) %>% ungroup() %>% 
  mutate(State = factor(State)) %>% #filter(!State %in% c(14, 30)) %>% 
  arrange(State)

analysis_data <- analysis_data %>% 
  filter(!State %in% c('14', '30')) %>% 
  mutate(State = droplevels(State))


# regressions by state
trend <- sqrt_count_2012 ~ sqrt_count_2010 - 1 + 
  Associate_Count +
  Office_Employee_Count

by_state <- analysis_data %>% 
  group_by(State) %>% 
  nest() 
nis <- apply(by_state,1, function(dd) dd$data %>% nrow())
state_lm <- function(df){
  lm(trend, y = TRUE, x = TRUE, data = df)
}

models_lm <- by_state$data %>% 
  map(state_lm) 
names(models_lm) <- by_state$State

#(nGroups <<- length(models_lm))
step_logbstar <- abs(log(mu_bstr/(sqrt(mu_bstr*(1-mu_bstr)/(psi_bstr+1))))) #abs log(mean/sd)
mu_rho_step <- .3
psi_rho_step <- (a_psir/b_psir^2)^.5 #mean/sd
rho_step <- .1
sigs2 <- unlist(models_lm %>% map(sigma))^2
step_Z <- abs(brlm:::fn.compute.Z(mean(sigs2), a0, b0)/(sqrt(nis)))

# convert response and design matrix for training set to lists for brlm functions
y <- models_lm %>% 
  map(.f = function(m) m$y)
X <- models_lm %>% 
  map(.f = function(m) m$x)


# restricted

(timing_bench <- rbenchmark::benchmark(
    "restricted" = {
      restricted <- brlm::hierNormTheoryRestLm(y,
                                               X,
                                               regEst = 'Tukey',
                                               scaleEst='Huber',
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
                                               step_Z)
  },
  
    "abc" = {
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
  }, replications = 5))

saveRDS(timing_bench, file = file.path(getwd(), "timing_bench.rds"))



### manual
nburn <- 0
nkeep <- 1000
strt_restricted <- Sys.time()
restricted <- brlm::hierNormTheoryRestLm(y,
                                         X,
                                         regEst = 'Tukey',
                                         scaleEst='Huber',
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
                                         step_Z)
end_restricted <- Sys.time() - strt_restricted



(strt_abc <- Sys.time())
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
(end_abc <- Sys.time() - strt_abc)

format(end_abc)
format(end_restricted)

