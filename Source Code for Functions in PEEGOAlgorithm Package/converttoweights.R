#' Function to convert surfactants to weights
#'
#' Converts the surfactants a row in the design matrix from 0 to 1 scaling to relative weights using transformation on page 130 of Atkinson et al (2007)
#'
#' @param drow row in the design matrix, on 0 to 1 scale
#' @param srow row in the indicator matrix for surfactants which corresponds to drow
#'
#' @return The weights for the three active (non-zero) surfactants
#' @export
#'
#' @examples converttoweights(drow=c(0.1, 0.2, 0.3, 0, 0, 1, 1, 1), srow=c(1,1,1,0,0))
#' # [1] 0.07142857 0.28571429 0.64285714
#'
#' @references Atkinson, A. C., Donev, A. N. and Tobias, R. D. (2007) 'Optimum Experimental Designs, with SAS.' Oxford University Press, Oxford.
converttoweights<-function(drow, srow) {
  drow<-drow[which(srow==1)] # Only selecting surfactants which are active (which elements of the row of s are non-zero)
  if(all(drow==0)) {  # If all the surfactants are zero, then we assume they have equal weighting.
    drow.conv<-rep(1/length(drow), length(drow))
  } else { # Otherwise, we use the conversion given in the aforementioned book.
    drow.conv<-(drow^2)/sum(drow^2)
  }
  return(drow.conv) # We return the converted drow, which will now have the same number of elements as there are active surfactants.
}

