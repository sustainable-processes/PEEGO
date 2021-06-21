#' Function to generate random rows
#'
#' Function used to generate the random rows considered in the EGO algorithm
#'
#' @param nrow number of random rows
#' @param s.index index of active (non-zero) surfactants for random starting rows
#' @param lb lower bounds for controllable variables, on 0 to 1 scale
#' @param ub upper bounds for controllable variables, on 0 to 1 scale
#'
#' @return A matrix of random rows which are considered for swapping for a certain row in the design.
#' @export
#'
#' @examples set.seed(123) # Set random seed for reproducible example
#' randomrow(nrow=5, s.index=c(1,0,1,0,1), lb=c(rep(0, 6), rep(0.05, 2)), ub=rep(1,8))
#'#          [,1] [,2] [,3] [,4] [,5]      [,6]      [,7]       [,8]
#'#[1,] 0.4089769    0    0    0    0 0.8830174 0.9434439 0.09327867
#'#[2,] 0.5514350    0    0    0    0 0.4566147 0.9589917 0.48066745
#'#[3,] 0.1029247    0    0    0    0 0.8998250 0.2837833 0.08995656
#'#[4,] 0.8895393    0    0    0    0 0.6928034 0.6584815 0.99455629
#'#[5,] 0.5440660    0    0    0    0 0.5941420 0.3247018 0.18975796
#' @seealso EGOAlg()
randomrow<-function(nrow, s.index, lb, ub) {
  random.row<-matrix(0, nrow=nrow, ncol=length(lb)) # Initialising the set of unique rows
  for(i in 1:nrow(random.row)) { # Loop over all added rows
    random.row[i,c(s.index, 6:8)]<-runif(6, lb[c(s.index, 6:8)], ub[c(s.index, 6:8)]) # Generating a random row, where the active surfactants are given by s.index, and the other surfactants are 0
  }
  return(random.row) # Returning these random rows
}
