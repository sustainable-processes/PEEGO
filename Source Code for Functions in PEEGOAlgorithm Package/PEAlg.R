#' Point exchange algorithm
#'
#' Point exchange algorithm for finding optimal designs for the formulated product experiment. Points in algorithm optimised using EGO algorithm. For more information on point exchange algorithms see Federov (1972)
#'
#'
#' @param dands a design and the indiciator matrix corresponding to this design combined together, so that the first 8 columns are the design and the final 5 columns are the inidicator matrix
#' @param s.indexes the matrix of inidces for all possible combinations of active surfactants (found using combn() function)
#' @param B number of Monte-Carlo simulation in utility function
#' @param Snum number of surfactants
#' @param w weights in objective function
#' @param Xrow number of random starting rows considered in EGO algorithm
#' @param maxitsEGO maximum number of iterations in EGO algorithm
#' @param lb lower bounds for controllable variables, on 0 to 1 scale
#' @param ub upper bounds for controllable variables, on 0 to 1 scale
#' @param cov.type covariance for Gaussian process in EGO algorithm
#' @param sig level of significance in Kolmogorov-Smirnov test for comparison of objective functions
#' @param tol tolerance for average absolute difference between the initial and final design in each loop of the algorithm
#' @param maxitsPE maximum number of iterations allowed in the point exchange algorithm
#'
#' @return A list with the following elements:
#' \item{final.design}{the optimal design}
#' \item{final.active.s}{the inidicator matrix for active surfactants in the optimal design}
#' \item{ks.test.pval}{the p-value for the Kolmogorov-Smirnov test comparing the objective function for the initial design and the design found using the point exchange algorithm}
#' \item{mean.initial.obj}{the average of the objective function for the initial design}
#' \item{mean.final.obj}{the average of the objective function for the design found using the point exchange algorithm}
#' \item{initial.utility}{the Monte-Carlo simulations of the utility function for initial design}
#' \item{final.utility}{the Monte-Carlo simulations of the utility function for design found using the point exchange algorithm}
#' \item{initial.obj}{the simulations of the objective function for the initial design}
#' \item{final.obj}{the simulations of the objective function for the design found using the point exchange algorithm}
#' \item{initial.dist}{the Euclidean distances between all possible pairs of points in the initial design}
#' \item{final.dist}{the Euclidean distances between all possible pairs of points in the design found using the point exchange algorithm}
#' \item{its}{the number of iterations of the point exchange algorithm performed}
#' @export
#'
#' @seealso EGOAlg()
#'
#' @examples # None given as relies on the utility function, which has to be defined using experimental data
#'
#' @references Fedorov, V. V. (1972) 'Theory of Optimal Experiments' No. 12 in Probability and Mathematical Statistics. Academic Press, New York.
PEAlg<-function(dands, s.indexes, B, Snum, w, Xrow, maxitsEGO, lb, ub, cov.type, sig, tol, maxitsPE) {
  # Seperating d and s
  # NB: these are provided as one argument because we want to use the function mclapply to parallelise this code
  d<-dands[,1:length(lb)]
  s<-dands[,(length(lb)+1):ncol(dands)]
  # Calculating objective function for initial design
  initial.utility<-util$utility(d=d, B=B) # Utility for initial design
  initial.utility<-initial.utility[which(initial.utility>=0)] # Removing utility values which are less than zero
  initial.dist<-designdist(d=d, s=s, Snum=Snum) # Distances for initial design
  initial.meandist<-mean(initial.dist, na.rm=T) # Average distance for initial design
  initial.obj<-(w*log(initial.utility))+((1-w)*log(initial.meandist)) # Objective function for initial design
  # Initialising values for while loop, to ensure it is true for the first iteration
  abs.diff<-tol+1
  kstest.pval<-sig
  its<-0
  # Point exchange loop
  while(all((mean(abs.diff)>=tol),(kstest.pval<=sig),its<maxitsPE)) { # This while loop will continue whilst the mean absolute difference between all values in the start and end design in each loop is greater than or equal to the tolerance,
    # and the pvalue for the KS test comparing the objective functions for the design at the start and end of each loop is less than or equal to the level of significance, and the number of iterations has not exceed maxitsPE.
    if(its==0) { # Setting start design and indicator matrix for first iteration, which are d and s, respectively
      EGO.d<-d
      EGO.s<-s
    }
    EGO.d.start<-EGO.d # Setting the starting design for this loop
    EGO.s.start<-EGO.s # Setting the starting indicator matrix for this loop
    start.utility<-util$utility(d=EGO.d.start, B=B) # Utility at start of current loop
    start.utility<-start.utility[which(start.utility>=0)] # Removing negative expected utility values
    start.meandist<-mean(designdist(d=EGO.d.start, s=EGO.s.start, Snum=Snum)) # Average distance at start of current loop
    start.loop.obj<-(w*log(start.utility))+((1-w)*log(start.meandist)) # Objective at start of current loop
    for(i in 1:nrow(d)) { # Looping over all possible rows in the design
      start.utility<-util$utility(d=EGO.d, B=B) # Utility before current row is considered
      start.utility<-start.utility[which(start.utility>=0)] # Removing negative expected utility values
      start.meandist<-mean(designdist(d=EGO.d, s=EGO.s, Snum=Snum)) # Average distance before current row considered
      start.obj<-(w*log(start.utility))+((1-w)*log(start.meandist)) # Objective before current row considered
      # Run EGO algorithm for current row and all possible combinations of active surfactants and perform KS test
      opt.row.sindex<-list() # Initialise list of results from EGO algorithm for current row and surfactant combinations
      ks.test.pvals<-rep(0, ncol(s.indexes)) # Initialise vector of p-values from KS test
      for(j in 1:ncol(s.indexes)) { # Loop over all possible surfactant combinations
        opt.row.sindex[[j]]<-EGOAlg(Xrow=Xrow, maxits=maxitsEGO, s.index=s.indexes[,j], Snum=Snum, rownum=i, d=EGO.d, s=EGO.s, lb=lb, ub=ub, B=B, w=w, cov.type=cov.type, sig=sig) # Run EGO algorithm for current row and surfactant combination
        obj.sindex<-(w*log(opt.row.sindex[[j]]$optimal.EU))+((1-w)*log(opt.row.sindex[[j]]$optimal.meandist)) # Calculate objective function for design with current row swapped for optimal row from EGO algorithm for current row and surfactant combination
        ks.test.pvals[j]<-ks.test(obj.sindex, start.obj)$p.value # P-value for KS test between objective function for original design before current row swap considered, and design with swapped row
      }
      # Find row with smallest KS test p-value
      min.index<-which.min(ks.test.pvals) # Index of which surfactant combinations minimises the ks test results
      if(opt.row.sindex[[min.index]]$optimisation=="EGO") { # Only swap row in current design if it has actually been changed by the EGO algorithm, as if this argument is "Original", then row hasn't been amended by the algorithm
        opt.row<-opt.row.sindex[[min.index]]$optimal.row # Row which has smallest p-value for KS test comparison with design before swap considered, and hence best row to swap in the original design
        opt.sindex<-s.indexes[,min.index] # Index for active surfactants for the optimal row
        opt.srow<-rep(0, Snum) # Initialise row of inidicator matrix for active surfactants of optimal row
        opt.srow[opt.sindex]<-rep(1,3) # Set surfactants at opt.sindex to be 1
        opt.obj<-(w*log(opt.row.sindex[[min.index]]$optimal.EU))+((1-w)*log(opt.row.sindex[[min.index]]$optimal.meandist)) # Objective function for the optimal row
        # Consider swap
        if((min(ks.test.pvals)<sig)&(mean(opt.obj)>mean(start.obj))) { # If the KS test p-value is less than sig and the mean of the objective function for the design with the optimal row swapped for the row in the start design for this loop is greater than the mean of the objective function for the start design in this loop, then:
          EGO.d[i,]<-opt.row # The row at index i in the design from the start of this loop is swapped for the optimal row, and
          EGO.s[i,]<-opt.srow # The row at index i in the indicator matrix from the start of this loop is swapped for the optimal row of s.
        } # Otherwise no swap is performed and the design remains unchanged.
      }
    }
    end.utility<-util$utility(d=EGO.d, B=B) # Utility following one loop of the algorithm, after all rows in the design have been considered
    end.utility<-end.utility[which(end.utility>=0)] # Removing negative values from expected utility simulations
    end.dist<-designdist(d=EGO.d, s=EGO.s, Snum=Snum) # Distance for design following one loop of algorithm
    end.meandist<-mean(end.dist, na.rm=T) # Mean distance for design following one loop of algorithm
    end.loop.obj<-(w*log(end.utility))+((1-w)*log(end.meandist)) # Objective function following one loop of algorithm
    abs.diff<-abs(c(EGO.d.start-EGO.d)) # Absolute difference between all pairs of controllable variable values in the start and end design from this loop of the algorithm
    kstest.pval<-ks.test(start.loop.obj, end.loop.obj)$p.val # P-value of KS test of the objective function for the designs from the start and end of this loop
    its<-its+1 # Increasing the number of iterations
    print(its) # Printing the current number of iterations, as a method of tracking the progress of the algorithm
  }
  # Calculating objective function for final design found using point exchange algorithm
  final.utility<-util$utility(d=EGO.d, B=B) # Utility for final design
  final.utility<-final.utility[which(final.utility>=0)] # Removing any negative utility values
  final.dist<-c(designdist(d=EGO.d, s=EGO.s, Snum=Snum)) # Distances between points in final design
  final.meandist<-mean(final.dist, na.rm=T) # Average distance between points in final design
  final.obj<-(w*log(final.utility))+((1-w)*log(final.meandist)) # Objective function for final design
  kstest<-ks.test(initial.obj, final.obj) # KS test for comparing the initial design (d) and the final design found using the point exchange algorithm
  # Determining output for function
  if((kstest$p.val<sig)&(mean(final.obj)>mean(initial.obj))) { # If the p-value from the KS test is less than the level of significance, sig, and the mean of the objective function for the design found using the point exchange algorithm is less than the mean of the objective funciton for the original design, then:
    final.d<-EGO.d # the design found using the point exchange is the final design, and
    final.s<-EGO.s # the indicator matrix found using the point exchange is the final indicator matrix.
  } else { # Otherwise,
    final.d<-d # the initial design, and
    final.s<-s # the initial indicator matrix is returned.
  }
  return(list(final.design=final.d, final.active.s=final.s, ks.test.pval=kstest$p.val, mean.initial.obj=mean(initial.obj), mean.final.obj=mean(final.obj), initial.utility=initial.utility, final.utility=final.utility, initial.obj=initial.obj, final.obj=final.obj, initial.dist=initial.dist, final.dist=final.dist, its=its))
}

