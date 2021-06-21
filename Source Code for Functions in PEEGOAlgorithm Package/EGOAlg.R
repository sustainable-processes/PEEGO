#' Row optimisation using Efficient Global Optimisation (EGO) algorithm
#'
#' Optimisation method used to chose row which is added to the design in the point exchange algorithm, based on the EGO algorithm from Jones et al (1998)
#'
#' @param Xrow number of random starting rows considered
#' @param maxits maximum number of iterations for EGO algorithm
#' @param s.index index of active (non-zero) surfactants for random starting rows
#' @param Snum number of surfactants
#' @param rownum index of row to be swapped
#' @param d a design
#' @param s an indicator matrix for the surfactants which corresponds to d
#' @param lb lower bounds for controllable variables, on 0 to 1 scale
#' @param ub upper bounds for controllable variables, on 0 to 1 scale
#' @param B number of Monte-Carlo simulations in utility function
#' @param w weights in objective function
#' @param cov.type covariance for Gaussian process used in EGO algorithm
#' @param sig level of significance for Kolmogorov-Smirnov test for comparison of objective functions (Smirnov, 1939)
#'
#' @return A list with the following elements:
#' \item{optimal.row}{optimal row to swap for the row at rownum in d}
#' \item{optimisation}{optimisation method; "Original" if the row which maximises the objective function found using the EGO algorithm does not improve on the original row at rownum in d, and "EGO" if the row which maximises the objective function found using the EGO algorithm is an improvement on the row at rownum of d}
#' \item{optimal.EU}{Monte-Carlo simulations of the utility values for the design where the row at rownum is swapped with the optimal row}
#' \item{optimal.meandist}{mean distance for the design where the row at rownum is swapped with the optimal row}
#' \item{its}{number of iterations in the EGO algorithm}
#' @export
#'
#' @examples # None given as relies on the utility function, which has to be defined using the experimental data
#'
#' @references \itemize{
#' \item Jones, D. R., Schonlau, M. and Welch, W. J. (1998) 'Efficient Global Optimization of Expensive Black-Box Functions.' Journal of Global Optimization, 13, 455-492.
#' \item Smirnov, N. V. (1939) 'Estimate of Deviation Between Empriical Distribution Functions in Two Independent Samples.' (Russian). Bulletin Moscow Univerity, 2, 3-16.
#' }

EGOAlg<-function(Xrow, maxits, s.index, Snum, rownum, d, s, lb, ub, B, w, cov.type, sig) {
  # Setting up optimisation
  new.s<-rep(0, Snum) # Need to initialise the row which will be added to the s matrix for the active surfactants indexed by s.index
  new.s[s.index]<-rep(1, length(s.index)) # Setting the surfactants indexed by s.index to 1
  X<-randomrow(Xrow, s.index=s.index, lb=lb, ub=ub) # Generating Xrow random rows to be swapped with the current row at index rownum of d
  Z<-apply(X, 1, Zfunction, srow=new.s, rownum=rownum, d=d, s=s, B=B, Snum=Snum, w=w) # Calculating the reciprocal of the objective function (given by Zfunction) for all the designs with the random rows
  fit.mod<-km(~1, design=X[,c(s.index,6:8)], response=Z, covtype=cov.type, control=list(trace=FALSE)) # Fitting a Gaussian process between the designs with these random rows and the objective functions for designs.
  # NB: The function km is from DiceKriging package.
  # NB: control=list(trace=FALSE) argument stops this function printing out its progress.
  its<-0 # Initialising the number of iterations in the EGO algorithm.
  EIval<-0.1*min(Z) # Initialising the starting value for Expected Improvement (EI), so the while loop is entered.
  # Optimisation using EGO algorithm - adding rows which maximise the EI.
  while((abs(EIval)>abs(0.01*min(Z)))&(its<maxits)) { # This while loop defines the stopping rule for the algorithm.
    # Algorithm stops when the EI is less than 1% of current objective function value, or when the maximum number of iterations (maxits) is reached.
    # NB: absolute values used to avoid any issues with negative objective functions.
    sink(paste("EGOAlgB", B, "covtype", cov.type, ".txt", sep=""), append=T) # This saves output of this algorithm into a text file, and hence stops it showing up when running.
    # NB: This output file is not required, so is overwritten everytime the algorithm is run.
    max.EI<-max_EI(fit.mod, lower=lb[c(s.index,6:8)], upper=ub[c(s.index,6:8)]) # Using the Gaussian Process model to find a point which maximises the EI.
    # NB: The function max_EI is from the DiceOptim package.
    sink() # Closes file with unwanted output, with path given by first sink() function.
    Xopt<-rep(0, ncol(d)) # Initialising the new random row which will be added to the matrix of rows, X.
    Xopt[c(s.index, 6:8)]<-max.EI$par # Parameter values for point which maximises EI, added to the correct indicies.
    EIval<-max.EI$value # EI value for point which maximises EI
    X<-rbind(X, Xopt) # Adding the point which maximises EI to the X matrix
    Z<-apply(X, 1, Zfunction, srow=new.s, rownum=rownum, d=d, s=s, B=B, Snum=Snum, w=w) # Calculating the reciprocal of the objective function values for the updated X matrix
    fit.mod<-km(~1, design=X[,c(s.index,6:8)], response=Z, covtype=cov.type, control=list(trace=FALSE)) # Fitting the Gaussian Process for the updated X matrix
    its<-its+1 # Updating the number of iterations
  }
  utility.sim.origd<-util$utility(d=d, B=B) # Utility for initial design
  utility.sim.origd<-utility.sim.origd[which(utility.sim.origd>=0)] # Removing any utility values less than 0, as these cannot be logged
  meandist.origd<-mean(designdist(d=d, s=s, Snum=Snum)) # Average distance for initial design
  obj.origd<-(w*log(utility.sim.origd))+((1-w)*log(meandist.origd)) # Objective function for initial design
  # Getting the range of indices in X which have been added to the initial matrix of random rows.
  if(its==1) { # If the algorithm only did one iteration,
    EIrange<-(Xrow+1) # then the range is just Xrow+1,
  } else { # otherwise
    EIrange<-(Xrow+1):(Xrow+its) # the range is from Xrow+1 to the number of iterations, its.
  }
  # NB: this is done to ensure we do not get errors when finding the optimal point.
  opt.row.ei<-X[Xrow+which.min(Z[EIrange]),] # Row which maximised EI during the EGO algorithm and minimises the reciprocal of the objective function (and hence maximises the objective function).
  d.ei<-d # Initialising the design with the optimal row added
  d.ei[rownum,]<-opt.row.ei # Design with optimal row from the EGO algorithm swapped with original row
  utility.sim.eid<-util$utility(d=d.ei, B=B) # Utility for new design
  utility.sim.eid<-utility.sim.eid[which(utility.sim.eid>=0)] # Removing any utility values less than zero
  meandist.eid<-mean(designdist(d=d, s=s, Snum=Snum)) # Average distance for new design
  obj.eid<-(w*log(utility.sim.eid))+((1-w)*log(meandist.eid)) # Objective function for new design
  kst<-ks.test(obj.origd, obj.eid) # Using the Kolmogorov-Smirnov (KS) test (Smirnov, 1939) to compare the distributions of the objective functions for the two designs.
  # The null hypothesis of this test is the that the two samples being compared come from the same distribution.
  # Hence if the null is rejected in this case, there is evidence to suggest the objective function has changed by amedning this row.
  # Determining the output for the function (dependent on the results of KS test)
  if((kst$p.value<sig)&(mean(obj.eid)>mean(obj.origd))) { # If there is evidence to reject the null at significance level sig, then:
    opt.row<-opt.row.ei # Optimal row is that found in the EGO algorithm,
    opt.row.type<-"EGO" # Return a string stating that this row was found using the EGO algorithm,
    opt.EU<-utility.sim.eid # Expected utility is that for the design with opt.row.ei swapped for the original row in index rownum, and
    opt.mean.dist<-meandist.eid # Mean distance is that for the design with opt.row.ei swapped for the original row in index rownum.
  } else { # If the conditions in the if statement above are not met, then:
    opt.row<-d[rownum,] # Optimal row is the row at index rownum in the initial design,
    opt.row.type<-"Original" # Return a string stating this row was in the initial design,
    opt.EU<-utility.sim.origd # Expected utility is that for the initial design, and
    opt.mean.dist<-meandist.origd # Mean distance is that for the initial design.
  }
  return(list(optimal.row=opt.row, optimisation=opt.row.type, optimal.EU=opt.EU, optimal.meandist=opt.mean.dist, its=its))
}

