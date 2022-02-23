#' Calculation of frailty variance for the Gamma frailty distribution
#' 
#' 
#' This function calculates the frailty variance for any given value of Kendall's tau between [0,1]
#' 
#' @details 
#' The possible values of Kendall's tau under the Gamma distribution are [0,1]. 
#'
#'
#' @param tau the wanted tau value 
#'
#' @return Return the corresponds frailty variance, theta
#' @export
#'
#' @examples
theta.gamma<-function(tau){
  2*tau/(1-tau)
}
