### PE Algorithm for Multiple Random Starts - Simulation 4 ###

# NB: This code is the same as the code given in 'Finding a Design Using the Point Exchange Algorithm.R'

## Clearing workspace and setting working directory ##

rm(list=ls()) # This clears the workspace, ensuring that any variables saved in the environment do not effect our current code. 
setwd("/Users/esr1g08/Documents/Shampoo - LiWei") # Sets working directory - will need to be amended to wherever the dataset 'all_data.csv' is stored.

## Installing/ loading libraries required to run this code ##

library(acebayes) # Package for utility function 
library(MASS) # Required for a function used in the utility function 
library(DiceKriging) # Required for a function used in the EGO algorithm
library(DiceOptim) # Required for a function used in the EGO algorithm
library(parallel) # Package for parallelisation
library(PEEGOAlgorithm)

## Initialising parameters used ##

rstart<-50 # Number of random starts
runs<-48 # Number of runs in the designs 
lb<-c(rep(0,6), 0.05, 0.05) # Lower bounds of the variables in the design, scaled to be 0/ 0.05 (as P1/T1 can't be zero)
ub<-rep(1,8) # Upper bounds of the variables in the design
s.indexes<-combn(1:5,3) # All possible combinations of 3 of the 5 surfactants considered in this design, given as indexes 
B<-1000 # Number of simulations in the Monte-Carlo estimate of the utility function 
Snum<-5 # Number of surfactants 
w<-0.5 # Weight on each part of the objective function
Xrow<-50 # Number of random starting rows used in the EGO algorithm
maxitsEGO<-100 # Maximum iterations for the EGO algorithm 
cov.type<-"powexp" # Covariance function used in EGO algorithm
sig<-0.05 # Significance level for KS test in point exchange algorithm
tol<-0.01 # Tolerance for point exchange algorithm
maxitsPE<-20 # Maximum number of iterations that can be considered in point exchange algorithm
# NB: These values have been chosen based on running this code on the University of Southampton cluser, Iridis 4. We would recommend testing/ adjusting for your own computer/ cluster. 
# NB: To run this for more random starts, you could just update this code and submit it to the cluster muliple times

## Random starting designs and indicator matrices ##

dands<-list() # Initialising list of random starting designs and indicator matrices
for(i in 1:rstart) {
  dands[[i]]<-matrix(0, nrow=runs, ncol=length(lb)+Snum) 
  for(j in 1:runs) {
    sindex<-s.indexes[, sample(1:ncol(s.indexes), 1)] # Randomly sampling which surfactants are active
    dands[[i]][j,length(lb)+sindex]<-rep(1,3) # These columns of this row of dands contains 0 if the surfactants are not active, and 1 if the surfactants are active
    dands[[i]][j, c(sindex, 6:8)]<-runif(6, min=lb[c(sindex, 6:8)], max=ub[c(sindex, 6:8)]) # Sampling the scaled variabels from a uniform distribution, but leaving non-active surfactants as zero
  }
  colnames(dands[[i]])<-c("S1", "S2", "S4", "S5", "S7", "Sum", "P1", "T1", "S1", "S2", "S4", "S5", "S7")
}

## Utility function ##

data<-read.csv("all_data.csv") # Reading in past data
data<-data[-which(data$Viscosity=="nv"),] # Removing rows which have no value for Viscosity
data$Sum<-apply(data[,1:5], 1, sum) # Creating dummy variable for sum
data<-data[-which(data$Sum<=13),] # Removing rows where the sum is less than 13, as we want to scale this variable
data<-data[-which(data$T1>2),] # Removing rows where T1 > 2
data.var<-(1:ncol(data))[-c(8,9)] # Setting columns with variables, that is removing P2 and T2
data.rescaled<-data # Rescaling data to be between 0 and 1, this is the data that has been used when fitting the model 
for(i in 1:length(data.var)) {
  data.rescaled[,data.var[i]]<-(data[,data.var[i]]-min(data[,data.var[i]]))/(max(data[,data.var[i]])-min(data[,data.var[i]]))
}

# Utility function #

initial.glm<-glm(Output~S1+S5+P1+T1+Sum-1, data=data.rescaled, family=binomial("logit")) # This is the model that has been found to best fit the past data for this experiment. This was found using a forward selection algorithm. 
p<-ncol(data.rescaled)-2 # Number of parameters in the model 
active.param<-c(1,4,6,7,8) # Which parameters are active (non-zero) in the model
beta.mu<-rep(0, p) # Initialising the mean of the prior for these parameters. 
beta.mu[active.param]<-initial.glm$coeff # Setting the mean to be the value of the coefficients for the model if the parameters are active in the momdel
beta.sigma<-diag(0.01, p) # Initialsing the variance-covariance matrix for the prior, assuming covariances are zero
for(i in 1:length(active.param)) {
  for(j in 1:length(active.param)) {
    beta.sigma[active.param[i], active.param[j]]<-vcov(initial.glm)[i,j] # Setting the variances for the active parameters to be ther variances found in the model 
  }
}
# NB: Non-active parameters are assumed to have a variance of 0.01. This is a parameter that could be amended in future designs. 

prior<-function(B) { # Setting up the prior function for the utility
  beta.mat<-matrix(0, nrow=B, ncol=p) 
  for(b in 1:B) {
    beta.mat[b,]<-mvrnorm(n=1, mu=beta.mu, Sigma=beta.sigma) # The prior for the parameters in each simulation of the utility function is the normal distribution with mean beta.mu and variance beta.var
    # NB: The mvrnorm function comes from the MASS package
  }
  beta.mat
}

util<-utilityglm(formula=~S1+S2+S4+S5+S7+Sum+P1+T1-1, family=binomial(link="logit"), prior=prior, criterion="D", method="MC") # Utility function using the D criterion, found using the utilityglm function in the acebayes package
# NB: This function is used in both the EGO and point exchange algorithms
#     Different utility functions could be considered by changing the criterion, model or prior. This was seen as a good starting point.

## Running the point exchange algorithm ##

PEAlgoutput<-mclapply(dands, PEAlg, s.indexes=s.indexes, B=B, Snum=Snum, w=w, Xrow=Xrow, maxitsEGO=maxitsEGO, lb=lb, ub=ub, cov.type=cov.type, sig=sig, tol=tol, maxitsPE=maxitsPE, mc.cores=12)
save(PEAlgoutput, file="PEAlgoutput4.Rdata") # Saving output

