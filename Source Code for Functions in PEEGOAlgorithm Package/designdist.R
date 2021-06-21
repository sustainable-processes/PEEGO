#' Function for the distance between all rows in a design
#'
#' Calculates the Euclidean distance between all possible pairings of rows in the design
#'
#' @param d a design
#' @param s an inidcator matrix for the surfactants which corresponds to d
#' @param Snum the number of surfactants
#'
#' @return A vector of Euclidean distances for all possible pairs of rows in the design
#' @export
#'
#' @examples set.seed(123) # Set random seed for reproducible example
#' lb<-c(rep(0,6), 0.05, 0.05) # Lower bounds of the variables in the design
#' ub<-rep(1,8) # Upper bounds of the variables in the design
#' d<-matrix(0, nrow=8, ncol=8) # Initialising an 8 run design
#' s<-matrix(0, nrow=8, ncol=5)
#' for(i in 1:nrow(d)) { # For all the rows in d and s
#' sindex<-sample(1:5,3) # Index for active surfactants
#' s[i,sindex]<-rep(1,3) # Setting active surfactants to 1 in indicator matrix
#' d[i,c(sindex,6:8)]<-runif(6, min=lb[c(sindex,6:8)], max=ub[c(sindex,6:8)])
#' # Sampling values in d matrix from uniform distribution
#' }
#' designdist(d=d, s=s, Snum=5)
#' # [1] 10.845509  9.174131 10.341590 12.607289 16.118427 13.045025  5.601129  4.835845
#' # [9]  8.954437 13.297260 16.487215  9.529606 10.814740  5.123037 12.874674 15.535015
#' # [17]  8.972709  8.972521 13.157941 15.212356  7.140476  9.584821  4.312826 12.942136
#' # [25]  7.354380 15.276135 10.596318 11.925301
designdist<-function(d, s, Snum) {
  # Convert d to true values
  d.trans<-matrix(0, nrow=nrow(d), ncol=ncol(d)-1) # Initialising new matrix, which will have one less column than d as the sum of the surfactants is removed
  for(i in 1:nrow(d)) { # Loop over all possible rows in d
    d.trans[i,]<-converttotruevals(d[i,], s[i,], Snum) # Each row of d.trans is the row of the d matrix converted to original values, using the convert.to.truevals function
  }
  # Calculating Euclidean distance between all pairwise rows of d
  rowcomb<-combn(1:nrow(d), 2) # The combn function with these arguments gives all possible combinations of two rows of d
  # Each column of the rowcomb matrix gives the indices for the two rows whos distance we want to measure
  distance<-rep(0, ncol(rowcomb)) # Initialising distance vector
  for(i in 1:ncol(rowcomb)) { # Loop over all columns in rowcomb
    distance[i]<-sqrt(sum((d.trans[rowcomb[1,i],]-d.trans[rowcomb[2,i],])^2)) # Each element of distance is the Euclidean distance between a pair of rows in d.trans.
  }
  return(distance) # Returns the vector of distances
}
