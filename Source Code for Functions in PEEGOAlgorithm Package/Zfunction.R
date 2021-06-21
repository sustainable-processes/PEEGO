#' Objective function
#'
#' Objective function which is optimised in the EGO and point exchange algorithms. The weighted sum of the log of the estimated expected utility and the log of the average distance between all pairwise combinations of rows in the design. Loosley based on the objective function used to find bridge designs in Jones et al (2015)
#'
#' @param drow row to be swapped at rownum to d, on 0 to 1 scale
#' @param srow row to be swapped at rownum to s, which corresponds to drow
#' @param rownum index of row to be swapped
#' @param d a design
#' @param s an indicator matrix for the surfactants which corresponds to d
#' @param Snum number of surfactants
#' @param B number of Monte-Carlo simualtions in utility function
#' @param w weight on log of expected utility estimate (weight on log of average pairwise distances is 1-w)
#'
#' @return Reciprocal of objective function value, which is used in EGO algorithm. EGO algorithm is finds the minimum of a function, and we want to maximise the objective function, hence the use of the reciprocal.
#' @export
#'
#' @examples # None given as relies on the utility function, which has to be defined using the experimental data
#'
#' @references Jones, B., Silvestrini, R. T., Montgomery, D. C. and Steinberg, D. M. (2015) 'Bridge Designs for Modelling Systems with Low Noise.' Technometrics, 57, 155-163.
Zfunction<-function(drow, srow, rownum, d, s, Snum, B, w) {
  d[rownum,]<-drow # Swap the row at index rownum of d with drow
  s[rownum,]<-srow # Swap the row at index rownum of s with srow
  utility.sim<-util$utility(d=d, B=B) # Makes B Monte-Carlo samples of the utility function for design d
  meandist<-mean(designdist(d=d, s=s, Snum=Snum)) # Calculates the average difference between all possible pairs of points in the design
  output<-(w*log(mean(utility.sim)))+((1-w)*log(meandist)) # Objective function, the weighted sum of the mean of the Monte-Carlo samples from the utility, which is an estimate of the expected utility, and the log of the mean avereage distance.
  # NB: - Use the log so that any large differences in scaling do not affect the results.
  #     - w controls the weight on each component, the higher w is, the more weight there is on the expected utility.
  return(1/output)
}
