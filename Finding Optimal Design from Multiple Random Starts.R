### Finding Optimal Design from Multiple Random Starts ###

# This code demonstrates how an optimal design can be found from running the code in 'Finding a Design Using the Point Exchange Algorithm.R' multiple times. 

# Here we assume that the code has been run 10 times, so we have results for 500 random starting designs. 
# It is assumed that the output for these 10 pieces of code has been saved as 'PEAlgoutput1.R', 'PEAloutput2.R', ..., 'PEAlgoutput10.R.'
# Example code to perform these simulations and get these outputs is given in the 'Example Simulation - Multiple Random Starts' folder. 
# NB: The example code is for w=0.5. 

## Clearing workspace and setting working directory ##

rm(list=ls()) # Clears workspace
setwd("/Users/esr1g08/Documents/Shampoo - LiWei") # Sets working directory - will need to be amended to wherever the point exchange algorithm results are stored. 
library(PEEGOAlgorithm)

## Reading in algorithm results ##

PEAlgall<-list() # Setting up a list where all the results for all 500 random starts are going to be saved

for(i in 1:10) { # Looping over all 10 results - would need to be amended if the number of scripts run was changed. 
  load(paste("PEAlgoutput", i, ".Rdata", sep="")) # Loading the results from running the script - would need to be amended if file name changed. 
  # For example if i=2, then the file "PEAlgoutput2.Rdata" is loaded. 
  PEAlgall<-c(PEAlgall, PEAlgoutput) # Adds the current loaded PEAlgoutput to the list of all outputs. 
  # NB: This would need to be amended if the name of the object where the results are saved is changed from PEAlgoutput in the 'Finding a Design Using the Point Exchange Algorithm.R' script. 
}

## Finding and saving point exchange algorithm output and design matrix for optimal design ##

w<-0.5 # w used in simulations, set to 0.5 as this is what is used in the examples
filename<-"w0p5" # setting what we want to use in the filename 

obj.func.vals<-rep(0, length(PEAlgall)) # Initialising a vector to store the objective function values in
for(i in 1:length(PEAlgall)) {
  if(length(PEAlgall[[i]])!=12) { # If the PE algorithm hasn't successfully run for a random starting design...
    obj.func.vals[i]<-NA # ... then we return a NA result
    # The algorithm may successfully run for all the random starting designs, but this catch is needed just in case it doesn't. 
  } else {
    obj.func.vals[i]<-PEAlgall[[i]]$mean.final.obj # If the algorithm has successfully run, then we want the objective function for the design found by the algorithm. 
  }
}
summary(obj.func.vals) # Prints a summary of the objective function values, which will also give you the number of NA's (if any). 
# If the number of NA's is small (say less than 25), we are not too concerned as there are a number of numeric processes which could have errors in the point exchange algorithm. 
# If the number of NA's is large, then it would be worth double checking the code and seeing where the errors occur. 
optPE<-PEAlgall[[which.max(obj.func.vals)]] # Gives the point exchange algorithm output for the optimal design (the one which maximises the objective function). 
save(optPE, file=paste("optPE", filename, ".Rdata", sep="")) # Saves the point exchange algorithm for the optimal design as 'OptPEw0p5.Rdata' when filename="w0p5".  
optdesign<-matrix(0, nrow=48, ncol=7) # Initialising optimal design
for(i in 1:nrow(optdesign)) {
  optdesign[i,]<-converttotruevals(optPE$final.design[i,], optPE$final.active.s[i,], 5) # Converting each row in the final design to its true values and saving in optdesign
}
write.csv(optdesign, file=paste("OptimalDesign", filename, ".csv", sep="")) # Saves the matrix for the optimal design as 'OptimalDesignw0p5.csv' when filename="w0p5".

# These saved outputs could then be used to produce plots using the code in 'Comparative Plots for Optimal Designs with Different w.R.'