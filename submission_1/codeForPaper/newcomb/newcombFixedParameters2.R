#Specify the prior Parameters for the newcomb Analysis
#this file is to ensure the various analysis use the same prior distributions
#using a second set of hyperparameters from MacEachern's paper. these are much more informative


library(MASS)
##############################################################
eta<-23.6 #value of the mean from previous estimates
tau<-2.04  #value of the standard deviations from previous estimates
#the prior for sigma^2 is  an inverse gamma
#inverse gamma parameters (Note R's parameterization)
alpha<-5   #shape parameter
beta<-10  #this gives a mean of 10/4=2.5
#####################################################
#beta^2/((alpha-1)^2*(alpha-2))

data(newcomb)
n<-length(newcomb)
