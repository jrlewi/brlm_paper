# Script to set up prior distributions ----

library(MASS)
library(tidyverse)
library(MCMCpack)

prior_data <- read_rds(file.path(here::here(), 'data', 'prior_data.rds'))
prior_data <- prior_data %>% 
  mutate(sqrt_count_2008 = sqrt(Count_2008), 
         sqrt_count_2010 = sqrt(Count_2010))  %>% filter(Type == 1)


# pooled regression analysis ----
prior_fit <- MASS::rlm(sqrt_count_2010 ~ sqrt_count_2008 - 1 + 
                         Associate_Count +
                         Office_Employee_Count, 
                       scale.est = 'Huber', data =  prior_data)


# MASS::rlm(sqrt_count_2012 ~ sqrt_count_2010 - 1, scale.est = 'Huber', data =  analysis_data )

theme_set(theme_bw(base_family = 'Times'))
#+ xlim(c(0,150)) +  ylim(c(0, 130)) +
ggplot(prior_data, aes(x = sqrt_count_2008, y = sqrt_count_2010, 
                       col = Type)) + geom_point(size = .3) + 
  stat_smooth(method = "rlm", formula = y ~ x -1, 
              method.args = list(scale.est = 'Huber'), 
              size = .5, lty = 2, se = FALSE, col = 1) + 
  guides(col = guide_legend(title = 'Agency Type'))


summary(prior_fit)
beta_0 <- coef(prior_fit)
sigma_beta_0 <- vcov(prior_fit)
beta_var_inflate <- nrow(prior_data)
sigma2_hat <- prior_fit$s^2

a_0 <- 5
b_0 <- sigma2_hat*(a_0 - 1)

b_0/(a_0 - 1)
b_0^2/((a_0-1)^2*(a_0-2))

curve(dinvgamma(x, a_0,b_0), from=0, to=.01)


parms_prior <- list(beta_0 = beta_0, sigma_beta_0 = sigma_beta_0, 
                    sigma2_hat = sigma2_hat, 
                    a_0 = a_0, b_0 = b_0, 
                    beta_var_inflate = beta_var_inflate )
write_rds(parms_prior, "parms_prior.rds")



#hierarchcical regression analysis -----
by_state <- prior_data %>% 
  group_by(State) %>% 
  filter(n() >= 40) %>% #use state with plenty of data
  nest() 
  
state_model <- function(df){
  MASS::rlm(sqrt_count_2010 ~ sqrt_count_2008 - 1, scale.est = 'Huber', data = df, maxit = 100)
}

models <- by_state$data %>% 
  map(state_model)

betaHats <-
  models %>% 
  map(coef) %>%
  unlist()

var_betaHats <-
  models %>% 
  map(.f = function(m) vcov(m)) %>%
  unlist()



sigma2Hats <-
  models %>% 
  map(.f = function(m) summary(m)$sigma^2) %>%
  unlist()

n_is <- by_state$data %>% 
  map(nrow) %>% unlist()


plot(betaHats)
abline(h = beta_0)



plot(sigma2Hats)
abline(h = sigma2_hat)
sigma2Hats <- sigma2Hats #[-which(sigma2Hats > .01)] #get rid of large estimates
plot(sigma2Hats)
abline(h = sigma2_hat)


X <- model.matrix(prior_fit)
p <- ncol(X)



nStates <- nrow(by_state)



# prior for beta_i's
wts <- n_is/sum(n_is)
p <- nrow(vcov(prior_fit))

delta_is <- sapply(betaHats, FUN=function(x) x - beta_0)
delta_is <- matrix(delta_is, p, nStates)
dList<-list()
for(i in 1:ncol(delta_is)){
  dList[[i]]<-delta_is[,i]%*%t(delta_is[,i])
}

SigmaDelta<-matrix(rep(0, p*p), p,p)
for(i in 1:length(dList)){
  SigmaDelta <- SigmaDelta+dList[[i]]*wts[i]
}

SigmaDelta  
K <-length(n_is)

vcov(prior_fit)

g <- (det(SigmaDelta)/(det(vcov(prior_fit))))^(1/p)
swSq <- sum(wts^2)



mu_b <- g/nrow(prior_data) #mean for b (i.e. mu_bstr)
swSq <- 1
psi_b <- 10 #( psi_bstr)

fn.compute.ab<-function(mu, psi){
  a<-mu*psi
  b<-psi-a
  c(a,b)
}
ab <- fn.compute.ab(mu_b, psi_b)
curve(dbeta(x, ab[1], ab[2]))


# -------------------------
# prior for sigma2_i's
# -------------------------


#given z~N(0,1) compute inverse gamma(ao, b0) random variable
# invGam<-function(z, a0, b0){
#   c<-pnorm(z)
#   x<-qgamma(c, a0, scale=1/b0)
#   1/x 
# }

curve(dinvgamma(x, a_0,b_0), from=0,to=.01)
points(sigma2Hats,rep(.1,length(sigma2Hats)))



#prior parameters for mu_rho ~beta(w1,w2)
#idea: compute the z_i's; calculate the MLE of rho, make this the mean of mu_rho

fn.compute.Z<-function(sigma2,a_0,b_0){
  X<-integrate(function(x) {dinvgamma(x, shape = a_0, scale = b_0)},0, sigma2)$value #rate=b0
  qnorm(X)
}
Z <- sapply(sort(sigma2Hats), FUN=fn.compute.Z, a_0 = a_0, b_0 = b_0)
#Z <- c(Z, rnorm(100, Z[1], 0.0001))
# Z <- (Z - mean(Z))/sd(Z)
K <- length(Z)
J <- matrix(1, K,K)

log.like.rho<-function(rho){
  Sigma_rho <- ((1-rho)*diag(K) + rho*J)
  Sigma_rho_inv <- solve(Sigma_rho)
  K*log(1-rho)+log(1+K*rho/(1-rho))+t(Z)%*%Sigma_rho_inv%*%Z
}

log.like.rho2 <-function(rho){
  Sigma_rho <- ((1-rho)*diag(K) + rho*J)
  Sigma_rho_inv <- solve(Sigma_rho)
  .5*log(det(2*pi*Sigma_rho)) + .5*t(Z)%*%Sigma_rho_inv%*%Z
}



#curve(-log.like.rho2(x),0.001,.999)
rho_seq <- seq(0, .99, by = .01)
#plot(rho_seq ,sapply(rho_seq , function(rho) log.like.rho(rho)))
plot(rho_seq ,sapply(rho_seq , function(rho) log.like.rho2(rho)))

curve(log.like.rho,0.001,.999)
#maximize
# op<-optim(.5, log.like.rho,control=list(fnscale=-1), lower=0, upper=.9999, method = c("L-BFGS-B"), hessian=TRUE)
#MLE of rho 
(op <- optim(.5, log.like.rho2, lower=0, upper=.999, method = c("L-BFGS-B"), hessian=TRUE))


log.like.rho<-function(rho){
  Sigma_rho <- ((1-rho)*diag(K) + rho*J)
  Sigma_rho_inv <- solve(Sigma_rho)
  K*log(1-rho)+log(1+K*rho/(1-rho))+t(Z)%*%Sigma_rho_inv%*%Z
}
(op<-optim(.5, log.like.rho, lower=0.001, upper=.9999, method = c("L-BFGS-B"), hessian=TRUE))


mean_mu_rho <- as.numeric(op$par)
var_mu_rho <- -2*as.numeric(op$hessian^-1)
#had trouble identifying any correlation in the sigma's in the prior data set

mean_mu_rho <- .2
var_mu_rho <-  .01


w1_plus_w2 <- mean_mu_rho*(1-mean_mu_rho)/var_mu_rho-1

var_mu_rho <- ((w1_plus_w2+1)/(mean_mu_rho*(1-mean_mu_rho)))^-1
w1 <- mean_mu_rho*w1_plus_w2
w2 <- w1_plus_w2-w1
w1;w2
curve(dbeta(x, w1,w2))
w1/(w1+w2)
w1*w2/((w1+w2)^2*(w1+w2+1))
quantile(rbeta(1e6, w1,w2), prob = c(.1, .9))
#prior parameters for  psi_rho ~ gamma(a_psir, b_psir)

mean_psi_rho <- (mean_mu_rho*(1-mean_mu_rho))/(var_mu_rho) - 1

a_psir <- mean_psi_rho 
b_psir<- 1 

curve(dgamma(x, a_psir, b_psir), from=0, to=100)
a_psir/b_psir
a_psir/b_psir^2

#samples of rho...

nsamps <- 10000
psi_rho_samp <- rgamma(nsamps, a_psir, b_psir)
mu_rho_samp <- rbeta(nsamps, w1,w2)
a_rho_samp <- mu_rho_samp*psi_rho_samp
b_rho_samp <- psi_rho_samp - a_rho_samp

rho_samp <- rbeta(nsamps, a_rho_samp, b_rho_samp)
hist(rho_samp)
quantile(rho_samp, probs = c(.1, .9))
#results in strong correlation amongst the sigmas. 

hier_parms_prior <- list(mu_0 = beta_0,
                         Sigma_0 = nrow(prior_data)*se_beta_0^2,
                         sigma2_hat = sigma2_hat,
                         a_0 = a_0, 
                         b_0 = b_0,
                         mu_bstr = mu_b,
                         psi_bstr = psi_b,
                         w1 = w1,
                         w2 = w2,
                         a_psir = a_psir,
                         b_psir = b_psir,
                         beta_var_inflate = nrow(prior_data), #beta_var_inflate,
                         prior_fit = prior_fit)
write_rds(hier_parms_prior, "hier_parms_prior.rds")

