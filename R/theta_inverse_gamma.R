
#' Title
#'
#' @param data 
#' @param tau 
#'
#' @return
#' @export
#'
#' @examples
theta.inverse.gaussian<-function(data,tau){
  theta<-data$theta[which.min(abs(data$tau-tau))]
  tau<-data$tau[which.min(abs(data$tau-tau))]
  
  return(list(theta=theta,tau=tau))
}