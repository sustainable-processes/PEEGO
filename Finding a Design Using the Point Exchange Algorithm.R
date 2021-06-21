### Finding Design Using Point Exchange Algorithm ###

# We are going to use the point exchange algorithm to find an optimal design. 
# This design is based on the bridge designs discussed in Jones et al (2015), but has been adapted slightly to better suit this specific problem. 
# The objective function used for this design is a weighted sum of the log of the average simulated utility function and the average distance between all pairwise points in the design. 
# The utility function considers the binary response from the stability test, and the distance considers the continuous response from the viscosity test. 
# It is assumed, following the analysis of data from past experiments, that a logit model is suitable for the outcomes from the stability test. 
# As the form of the model for the response from the viscocisty test is unknown, a Gaussian Process model is assumed for these outcomes. 

# Further information about any of the functions can be found by typing ?functionname into R and pressing enter

## Clearing workspace and setting working directory ##

rm(list=ls()) # This clears the workspace, ensuring that any variables saved in the environment do not effect our current code. 
setwd("/Users/esr1g08/Documents/Shampoo - LiWei") # Setting working directory, will need to change this to whatever directory your inputs are currently saved in, and this will be the directory that any outputs will also be saved into

## Installing/ loading libraries required to run this code ##

# Installing libraries #
# If you are running this code for the first time, you will need to install various libraries to ensure this code will run. 
# NB: Once you have done this, you will not neeed to run this again. 
install.packages("PEEGOAlgorithm_0.1.tar.gz", repos=NULL) # This assumes that this file is in the working directory set above
install.packages(c("acebayes", "MASS", "DiceKriging", "DiceOptim", "parallel"))

# You may be asked to select a CRAN mirror for downloading the packages, I usually select Bristol or London, but it doesn't matter which UK one you choose. It is essentially just asking what online repositry you want to use. 

# Loading libraries #
# In order to use the functions in an R library, you need to load them. 
library(acebayes) # Package for utility function 
library(MASS) # Required for a function used in the utility function 
library(DiceKriging) # Required for a function used in the EGO algorithm
library(DiceOptim) # Required for a function used in the EGO algorithm
library(parallel) # Package for parallelisation 
library(PEEGOAlgorithm) # Package for point exchange algorithm 

## Initialising parameters used ##

rstart<-50 # Number of random starts
runs<-48 # Number of runs in the designs 
lb<-c(rep(0,6), 0.05, 0.05) # Lower bounds of the variables in the design, scaled to be 0/ 0.05 (as P1/T1 can't be zero)
ub<-rep(1,8) # Upper bounds of the variables in the design
s.indexes<-combn(1:5,3) # All possible combinations of 3 of the 5 surfactants considered in this design, given as indexes 
B<-1000 # Number of simulations in the Monte-Carlo estimate of the utility function 
Snum<-5 # Number of surfactants (currently S1, S2, S4, S5, S7)
w<-0.5 # Weight on each part of the objective function
Xrow<-50 # Number of random starting rows used in the EGO algorithm
maxitsEGO<-100 # Maximum iterations for the EGO algorithm 
cov.type<-"powexp" # Covariance function used in EGO algorithm
sig<-0.05 # Significance level for KS test in point exchange algorithm
tol<-0.01 # Tolerance for point exchange algorithm
maxitsPE<-20 # Maximum number of iterations that can be considered in point exchange algorithm
# NB: These values have been chosen based on running this code for the maximum 60 hour job time on the University of Southampton cluser, Iridis 4. We would recommend testing/ adjusting for your own computer/ cluster. 
# NB: To run this for more random starts, you could just update this code and submit it to the cluster muliple times

## Random starting designs and indicator matrices ##

# Random starting designs are used as the initial designs in the point exchange algorithm.
# Random matrices which indicate which of the surfactants are active are also generated. 
# These are combined into a single matrix so that the code can be parallelised

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

# NB: The following adjustments would need to be made to this code should the data in 'all_data.csv' be amended:
#     - If S3, S6, S8, P2 or T2 (or any combination of these are added back in), then the column indices used when reading in the data will have to be amended. 
#     - If this data is changed, a new model will need to be fitted using the 'Model Selection for Prior.R' code. 
#     - If a new model is chosen, the indexes used in the active.param object will need to be changed. 

# Reading in data for prior # 
# The prior for the utility function is based on the past data that has been collected for this experiment

data<-read.csv("all_data.csv") # Reading in past data
data<-data[-which(data$Viscosity=="nv"),] # Removing rows which have no value for Viscosity
data$Sum<-apply(data[,1:Snum], 1, sum) # Creating dummy variable for sum
data<-data[-which(data$Sum<=13),] # Removing rows where the sum is less than 13, as we want to scale this variable
data<-data[-which(data$T1>2),] # Removing rows where T1 > 2
data.var<-(1:ncol(data))[-c(8,9)] # Column indices for variables we want to consider, removing the columns for the two outputs
data.rescaled<-data # Rescaling data to be between 0 and 1, this is the data that has been used when fitting the model 
for(i in 1:length(data.var)) {
  data.rescaled[,data.var[i]]<-(data[,data.var[i]]-min(data[,data.var[i]]))/(max(data[,data.var[i]])-min(data[,data.var[i]]))
}

# Utility function #

initial.glm<-glm(Output~S1+S5+P1+T1+Sum-1, data=data.rescaled, family=binomial("logit")) # This is the model that has been found to best fit the past data for this experiment. This was found using the algorithm in 'Model Selection for Prior.R'
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

PEAlgoutput<-mclapply(dands, PEAlg, s.indexes=s.indexes, B=B, Snum=Snum, w=w, Xrow=Xrow, maxitsEGO=maxitsEGO, lb=lb, ub=ub, cov.type=cov.type, sig=sig, tol=tol, maxitsPE=maxitsPE, mc.cores=12) # mclapply runs the function over the number of cores specfied by mc.cores
# NB: mc.cores=12 in this case as that is the number of cores on Iridis. If you were running this on your own computer then you would need to reduce this dependent on your processor. 
#     mcapply sets mc.cores=1 if run on windows, so will not run any faster 
save(PEAlgoutput, file="PEAlgoutputtest.Rdata") # Saving output
# This took between 52 and 60 hours when run on Iridis

## Plotting results for optimal design from the point exchange algorithm ##
# The optimal design is the design which maximisies the objective function, we begin by finding this design and then plotting the results for this design

# Finding optimal design #

obj.func.vals<-rep(0, length(dands)) # Calculating the objective funciton value for each element of the list dands
for(i in 1:length(dands)) {
  obj.func.vals[i]<-(w*log(mean(PEAlgoutput[[i]]$final.utility)))+((1-w)*log(mean(PEAlgoutput[[i]]$final.dist)))
}
optPE<-PEAlgoutput[[which.max(obj.func.vals)]] # Selecting outputs for optimal design


# Plotting results #

plotname="OptimalDesign"
# A density estimate, which is found using the function density(), is used in all of these plots. 
# The function plot() produces a plot, and lines() adds lines to that plot. The 'main' argument sets the title of the plot, 'xlim' sets the limits on the x-axis, 'ylim' sets the limits on the y-axis and 'col' determines the colour of any lines in the plot.
# The pdf() function opens a pdf file, which the plot will be saved into, and the dev.off() function closes this file. 
# The plot name argument will be used in the file name for all these plots
# Plot of simulations of the utility function for the initial design and the design found using the point exchange algorithm
pdf(paste("UtilityPlot", plotname, ".pdf", sep=""))
plot(density(optPE$initial.utility), main="", xlim=c(min(density(optPE$initial.utility)$x,density(optPE$final.utility)$x), max(density(optPE$initial.utility)$x,density(optPE$final.utility)$x)), ylim=c(min(density(optPE$initial.utility)$y,density(optPE$final.utility)$y), max(density(optPE$initial.utility)$y,density(optPE$final.utility)$y)))
lines(density(optPE$final.utility), col="blue")
dev.off()
# Plot of the simulation of the objective function for the initial design and the design found using the point exchange algorithm
pdf(paste("ObjectivePlot", plotname, ".pdf", sep=""))
plot(density(optPE$initial.obj), main="", xlim=c(min(density(optPE$initial.obj)$x,density(optPE$final.obj)$x), max(density(optPE$initial.obj)$x,density(optPE$final.obj)$x)), ylim=c(min(density(optPE$initial.obj)$y,density(optPE$final.obj)$y), max(density(optPE$initial.obj)$y,density(optPE$final.obj)$y)))
lines(density(optPE$final.obj), col="blue")
dev.off()
# NB: If the design found using the point exchange is better than the original design, then the blue line in the plot will be shifted to the right of the black line
#     This plot is only to be used in conjunction with the p value from the KS test, as the shift needs to be big enough to say the difference in distributions is statistically significant
# Plot of the pairwise distances between rows in the initial design and the design found using the point exchange algorithm
pdf(paste("DistancesPlot", plotname, ".pdf", sep=""))
plot(density(optPE$initial.dist, na.rm=T), main="", xlim=c(min(density(optPE$initial.dist, na.rm=T)$x,density(optPE$final.dist, na.rm=T)$x), max(density(optPE$initial.dist, na.rm=T)$x,density(optPE$final.dist, na.rm=T)$x)), ylim=c(min(density(optPE$initial.dist, na.rm=T)$y,density(optPE$final.dist, na.rm=T)$y), max(density(optPE$initial.dist, na.rm=T)$y,density(optPE$final.dist, na.rm=T)$y)))
lines(density(optPE$final.dist, na.rm=T), col="blue")
dev.off()

# Plotting the design #

optimal.design<-matrix(0, nrow=48, ncol=7)
colnames(optimal.design)<-c("S1", "S2", "S4", "S5", "S7", "P1", "T1")
for(i in 1:nrow(optimal.design)) {
  optimal.design[i,]<-converttotruevals(optPE$final.design[i,], optPE$final.active.s[i,], 5)
}
plotname="OptimalDesign"
# Pairs plots for the surfactants, and P1 and T1
surfactantpairs<-combn(1:5,2) # All pairs of surfactants, to be plotted
pdf(paste("DesignPlot", plotname, ".pdf", sep=""))
par(mfrow=c(4,3))
for(i in 1:ncol(surfactantpairs)) {
  plot(x=optimal.design[,surfactantpairs[1,i]], y=optimal.design[,surfactantpairs[2,i]], main="", xlab=colnames(optimal.design)[surfactantpairs[1,i]], ylab=colnames(optimal.design)[surfactantpairs[2,i]], xlim=c(0,15), ylim=c(0,15))
}
plot(x=optimal.design[,6], y=optimal.design[,7], main="", xlab=colnames(optimal.design)[6], ylab=colnames(optimal.design)[7], xlim=c(0,2), ylim=c(0,2))
dev.off()

# Saving optimal design 

write.csv(optimal.design, file="OptimalDesign.csv")

##########################################################

# Bibliography 
# Jones, B., Silvestrini, R. T., Montgomery, D. C. and Steinberg, D. M. (2015) 'Bridge Designs for Modelling Systems with Low Noise.' Technometrics, 57, 155-163.