#Specify the prior Parameters for the newcomb Analysis
#this file is to ensure the various analysis use the same prior distributions


library(MASS)

##############################################################
eta<--222.32 #value of the mean from previous estimates
tau<-540.4248  #value of the standard deviations from previous estimates
#the prior for sigma^2 is  an inverse gamma
#inverse gamma parameters (Note R's parameterization)
alpha<-5   #shape parameter
beta<-50*(alpha-1)   #this gives a mean of 50
#####################################################
#beta^2/((alpha-1)^2*(alpha-2))

data(newcomb)
n<-length(newcomb)

