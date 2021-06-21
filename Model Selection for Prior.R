### Model Selection for Prior ###

# This code provides the method used to find a model for the results of the stability test based on the current data set. 
# The parameters in this model will be used in the prior for the utility function in our objective function.

# This code assumes that S1, S2, S4, S5, S7, P1 and T1 are varied in the experiments. If any other controllable variables are to be included, this code will need amending.

## Clearing workspace and setting working directory ##

rm(list=ls()) # This clears the workspace, ensuring that any variables saved in the environment do not effect our current code. 
setwd("/Users/esr1g08/Documents/Shampoo - LiWei") # Setting working directory, will need to change this to whatever directory your inputs are currently saved in, and this will be the directory that any outputs will also be saved into

## Reading in data ##

data<-read.csv("all_data.csv") # Reading in data
data<-data[-which(data$Viscosity=="nv"),] # Removing rows which have no value for the viscosity test
data$Sum<-apply(data[,1:5], 1, sum) # Creating dummy variable for the sum of the surfactants 
# NB: The columns selected, denoted by [,1:5], will need to be amended if more than 5 surfactants are considered. 
data<-data[-which(data$Sum<=13),] # Removing rows where the sum is less than 13, as these do not meet the current restrictions on the controllable variable
data<-data[-which(data$T1>2),] # Removing rows where T1 > 2, as this is outside the current boundary for this variable

## Rescaling the data ##

# The point exchange algorithm assumes that the data is scaled to be between 0 and 1, therefore the parameters in the model used for the prior need to based on scaled data

data.var<-(1:ncol(data))[-c(8,9)] # Removing the columns relating to the two outputs
# NB: The columns removed, denoted by -c(8,9), would need to be amended if more surfactants are considered or if P2 and T2 are to be included again.
data.rescaled<-data # Initialising matrix for rescaled data
for(i in 1:length(data.var)) { # Loop to rescale variables
  data.rescaled[,data.var[i]]<-(data[,data.var[i]]-min(data[,data.var[i]]))/(max(data[,data.var[i]])-min(data[,data.var[i]])) # Data rescaled to be between 0 and 1
}

## Model Selection Algorithm ## 

# This is a forward selection algorithm which has three steps:
# 1) Fit models to each variable individually and idenitfy the effects with p-values less than a given level of signficance (set to 5% in this code) when 
# 2) Fit models which contain at least two of the effects identified in step 1, and calculate the AIC, BIC and deviance for these models. 
# 3) Identify the models which minimise the AIC, BIC and deviance, select the model with the most terms (as we are using this to build a prior). 

# Here we consider the two most commonly used models for binary (pass/fail) data, logistic and probit regression. 
# The function used to fit models to the stability test data is glm, which fits a generalised linear model from the family stated using the argument family. 
# Logistic regression is fitted using the argument family=binomial("logit") and probit regression is fitted using the argument family=binomial("probit"). 

# The AIC, BIC and deviance are considered as model selection crtieria as they are commonly used in forward selection algorithms. They are by no means an exhasutive list. 

# The model without the intercept is considered in all of these cases, and no joint effects (or interactons) are considered. 

# Output is used as our response throughout, as this is the name given to the stability test response in 'all_data.csv.'

# Step 1 - Fitting models with one variable and identifying the effects with p-values < 0.05 

glm.logit.pval<-glm.probit.pval<-NULL # Initialising vector that p-values are going to be saved into. 
for(i in 1:length(data.var)) { # Only considering controllable variables, and not the outputs
  glm.logit<-glm(Output~data.rescaled[,data.var[i]]-1, data=data.rescaled, family=binomial("logit")) # Fitting the logistic (shortened to logit here) regression model with one parameter
  # The glm function fits a generalised linear model, and has the following arguments:
  # - The model, in the form response~inputs. The -1 means we are fitting the model without the intercept. 
  # - The data the response and inputs can be found in. 
  # - The type of generalised linear model to be fitted, given by the family argument.
  glm.logit.pval[i]<-summary(glm.logit)$coeff[1,4] # Retrieving p-value from the summary for the glm model object (through the $coeff[1,4] argument) and assigning as an element of glm.logit.pval
  glm.probit<-glm(Output~data.rescaled[,data.var[i]]-1, data=data.rescaled, family=binomial("probit")) # Fitting the probit regression model with one parameter
  glm.probit.pval[i]<-summary(glm.probit)$coeff[1,4] # As above, but for probit regression
  # NB: The objects glm.logit and glm.probit are written over after each iteration
}

chosen.vals.logit<-which(glm.logit.pval<=0.05) # Gives indexes for parameters that have a p-value of less than 0.05 (5% significance) when a logistic regression model for that individual effect is considered
chosen.vals.probit<-which(glm.probit.pval<=0.05) # As above, but for probit regression 
all(chosen.vals.logit==chosen.vals.probit) # Check to see if the same terms have been selected. This doesn't have to be the case, but it makes the code for step two easier if it is, hence why it has been checked. 
# [1] TRUE 
# In this case, the chosen values are both the same
colnames(data.rescaled)[data.var[chosen.vals.logit]] # These are the individual parameters that have been identified in step 1
# [1] "S1"  "S5"  "P1"  "T1"  "Sum"

# Step 2 - Fit models which contain at least two of the effects identified in step 1 and calculate model selection crtieria

# Creating a list of the combinations of two or more elements of the vector of individual effects identified in step 1
combinations<-list() # Initialising the list 
for(i in 1:(length(chosen.vals.logit)-1)) {
  combinations[[i]]<-combn(chosen.vals.logit, i+1) # Creates a matrix of all the possible combinations of i+1 elements of the vector chosen.vals.logit
  # For example, combinations [[3]] is 
  #       [,1] [,2] [,3] [,4] [,5]
  # [1,]    1    1    1    1    4
  # [2,]    4    4    4    6    6
  # [3,]    6    6    7    7    7
  # [4,]    7    8    8    8    8
  # which is a matrix of all possible combinations of 4 elements of the indices in chosen.vals.logit
  #NB: If chosen.vals.logit and chosen.vals.probit were not identical (as checked in step 1), we would need a seperate list for these two vectors. 
}

# Fitting the models and saving the AIC, BIC and deviance 
glm.AIC<-glm.BIC<-glm.deviance<-NULL # Initialisng vectors that model selection criteria will be stored in
glm.AIC.fit<-glm.BIC.fit<-glm.deviance.fit<-list() # Initialising list that models will be saved in
# NB: glm.AIC[i] will give the AIC for the best model with i+1 individual effects and glm.AIC.model[[i]] will be that model. The same is also true for BIC and deviance. 
for(i in 1:length(combinations)) { # For all possible elements of the list combinations
  glm.probit<-glm.logit<-list() # Initialising a list to save the models
  glm.probit.AIC<-glm.logit.AIC<-glm.probit.BIC<-glm.logit.BIC<-glm.probit.dev<-glm.logit.dev<-NULL # Initialising vector of AIC/BIC/deviance
  # These lists and vectors will be written over at each iteration of the loop. 
  # Fitting the models and saving the AIC, BIC and deviance
  for(j in 1:ncol(combinations[[i]])) {
    glm.logit[[j]]<-glm(Output~.-1, data=data.rescaled[,c(data.var[combinations[[i]][,j]], 8)], family=binomial("logit")) # Fitting the logistic regression model to the parameters indexed by data.var[combinations[[i]][,j]].
    # The notation .-1 in Output ~ .-1 means that all columns in the argument data are considered in the model. 
    # For example, when i=1 and j=1, the model will be fitted will be Output ~ S1 + S5 -1
    # NB: If the matrix data.rescaled changes due to the addition of new varaibles, the argument 8 in c(data.var[combinations[[i]][,j]], 8) may have to be changed, as this is the index for the stability test output
    glm.logit.AIC[j]<-AIC(glm.logit[[j]]) # Calculating AIC for the logistic regression model
    glm.logit.BIC[j]<-BIC(glm.logit[[j]]) # Calculating BIC for the logistic regression model 
    glm.logit.dev[j]<-summary(glm.logit[[j]])$deviance # Extracting deviance from the summary for the logistic regression model
    glm.probit[[j]]<-glm(Output~.-1, data=data.rescaled[,c(data.var[combinations[[i]][,j]], 8)], family=binomial("probit")) # As above but for the probit regression model.
    glm.probit.AIC[j]<-AIC(glm.probit[[j]]) # The following three lines are as above but for probit regression model
    glm.probit.BIC[j]<-BIC(glm.probit[[j]])
    glm.probit.dev[j]<-summary(glm.probit[[j]])$deviance
  }
  # Identifying the best model with respect to AIC with i+1 effects 
  if(min(glm.logit.AIC)<min(glm.probit.AIC)) { # If the minimum AIC value for all the logistic regression models is less than that for the probit regression models, then the logistic regression model with i+1 effects which minimises AIC should be saved
    glm.AIC[i]<-min(glm.logit.AIC) # Minimum AIC value 
    glm.AIC.fit[[i]]<-glm.logit[[which.min(glm.logit.AIC)]] # Model with minimum AIC value
  } else { # Otherwise, the probit regression model with i+1 effects which minimises AIC should be saved
    glm.AIC[i]<-min(glm.probit.AIC) # Minimum AIC value 
    glm.AIC.fit[[i]]<-glm.probit[[which.min(glm.probit.AIC)]] # Model with minimum AIC value 
  }
  # Identifying the best model with respect to BIC with i+1 effects
  # NB: The following two if/else statements are as above, but for BIC and deviance instead of AIC. 
  if(min(glm.logit.BIC)<min(glm.probit.BIC)) { 
    glm.BIC[i]<-min(glm.logit.BIC)
    glm.BIC.fit[[i]]<-glm.logit[[which.min(glm.logit.BIC)]]
  } else { 
    glm.BIC[i]<-min(glm.probit.BIC)
    glm.BIC.fit[[i]]<-glm.probit[[which.min(glm.probit.BIC)]]
  }
  # Identifying the best model with respect to deviance with i+1 effects
  if(min(glm.logit.dev)<min(glm.probit.dev)) { 
    glm.deviance[i]<-min(glm.logit.dev)
    glm.deviance.fit[[i]]<-glm.logit[[which.min(glm.logit.dev)]]
  } else {  
    glm.deviance[i]<-min(glm.probit.dev)
    glm.deviance.fit[[i]]<-glm.probit[[which.min(glm.probit.dev)]]
  }
}
# Hence the ith element of the vector glm.AIC gives you the AIC value for the best model with i+1 effects, and the ith element of the list glm.AIC.fit gives you the model fit for the best model with i+1 effects (and similarly for BIC and deviance)
# For example:
glm.AIC[2] # AIC for model with 3 terms
#[1] 103.6777
glm.AIC.fit[[2]] # Model fit for model with 3 terms
# 
# Call:  glm(formula = Output ~ . - 1, family = binomial("logit"), data = data.rescaled[, 
#                                                                                       c(data.var[combinations[[i]][, j]], 8)])
# 
# Coefficients:
#   S1       T1      Sum  
# -0.7142  -1.2264   0.9204  
# 
# Degrees of Freedom: 77 Total (i.e. Null);  74 Residual
# Null Deviance:	    106.7 
# Residual Deviance: 97.68 	AIC: 103.7

# Step 3 - Find the models with optimise the model selection crtieria

summary(glm.AIC.fit[[which.min(glm.AIC)]]) # Model which minimises AIC
# 
# Call:
#   glm(formula = Output ~ . - 1, family = binomial("logit"), data = data.rescaled[, 
#                                                                                  c(data.var[combinations[[i]][, j]], 8)])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3848  -0.9292  -0.8301   1.4015   1.6746  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# T1   -1.2825     0.6747  -1.901   0.0573 .
# Sum   0.7697     0.8598   0.895   0.3707  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 106.745  on 77  degrees of freedom
# Residual deviance:  98.201  on 75  degrees of freedom
# AIC: 102.2
# 
# Number of Fisher Scoring iterations: 4
# 
summary(glm.BIC.fit[[which.min(glm.BIC)]]) # Model which minimises BIC
# 
# Call:
#   glm(formula = Output ~ . - 1, family = binomial("logit"), data = data.rescaled[, 
#                                                                                  c(data.var[combinations[[i]][, j]], 8)])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3848  -0.9292  -0.8301   1.4015   1.6746  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# T1   -1.2825     0.6747  -1.901   0.0573 .
# Sum   0.7697     0.8598   0.895   0.3707  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 106.745  on 77  degrees of freedom
# Residual deviance:  98.201  on 75  degrees of freedom
# AIC: 102.2
# 
# Number of Fisher Scoring iterations: 4
# 
summary(glm.deviance.fit[[which.min(glm.deviance)]]) # Model which minimises deviance
# 
# Call:
#   glm(formula = Output ~ . - 1, family = binomial("logit"), data = data.rescaled[, 
#                                                                                  c(data.var[combinations[[i]][, j]], 8)])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3807  -0.8971  -0.7417   1.3153   1.9406  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# S1   -0.8739     1.0283  -0.850    0.395
# S5   -0.6944     0.8653  -0.802    0.422
# P1    0.1618     0.9162   0.177    0.860
# T1   -1.2168     1.0379  -1.172    0.241
# Sum   0.9968     0.9388   1.062    0.288
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 106.74  on 77  degrees of freedom
# Residual deviance:  96.99  on 72  degrees of freedom
# AIC: 106.99
# 
# Number of Fisher Scoring iterations: 4

# We select the model which has the most terms in it for our prior, which is the model which minimises the deviance. 
