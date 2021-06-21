#' Function to convert to true values
#'
#' Converts a row in the design matrix from 0 to 1 scaling to true values
#'
#' @param drow row in the design matrix, on 0 to 1 scale
#' @param srow row in the indicator matrix for surfactants which corresponds to drow
#' @param Snum number of surfactants
#'
#' @return A row in the design matrix which has been converted to the true value
#' @export
#'
#' @examples converttotruevals(drow=c(0,0.5, 0.6, 0.7,0,0.8,0.9,1), srow=c(0,1,1,1,0), 5)
#' #[1] 0.000000 3.318182 4.778182 6.503636 0.000000 1.800000 2.000000
converttotruevals<-function(drow, srow, Snum) {
  drow.trans<-drow # Initialise a new vector
  drow.trans[(Snum+1)]<-13+(2*drow[(Snum+1)]) # Converting sum of the surfactants (which is in the Snum + 1 persition of the design matrix) to be between 13 and 15
  drow.trans[which(srow==1)]<-drow.trans[(Snum+1)]*converttoweights(drow[1:Snum], srow) # Converting surfactants to values which sum to the sum value, using the convert.to.weights function.
  drow.trans[(Snum+2):length(drow)]<-2*drow[(Snum+2):length(drow)] # Converting P1 & T1
  return(drow.trans[-(Snum+1)]) # Return drow without sum, as don't want this included in the distance calculation
}
