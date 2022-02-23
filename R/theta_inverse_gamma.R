#' Calculation of frailty variance for the Inverse Gaussian frailty distribution
#' 
#' 
#' This function calculates the frailty variance for any given value of Kendall's tau between [0,0.5)
#' 
#' @details 
#' The possible values of Kendall's tau under the inverse Gaussian distribution are [0,0.5). 
#' Since the function that connects between the Kendall's tau and theta is not invertable, a file with two
#' column : tau and the corresponds theta values needs to be specified. 
#'
#' @param data - data of Kendall's tau and corresponds theta values
#' @param tau -the wanted tau value
#'
#' @return Return the corresponds frailty variance, theta
#' @export
#'
#' @examples
theta.inverse.gaussian<-function(data,tau){
  theta<-data$theta[which.min(abs(data$tau-tau))]
  tau<-data$tau[which.min(abs(data$tau-tau))]
  
  return(list(theta=theta,tau=tau))
}