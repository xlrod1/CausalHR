#this the function that given the hazard function and the cumulative hazard calculate the causal HR
# This function can calculate the causal hazard ratio under two frailty distribution-
# 1=Gamma frailty distribution
# 2=Inverse Gaussian distribution

#' Calculation of the causal hazard ratio for different frailty distribution
#' 
#' 
#' This function calculate the causal hazard ratio for the Gamma or the Inverse Gaussian frailty distributions 
#' 
#' @details 
#' This function calculates the causal hazard ratio given the hazard rate and the cumulative hazard rate for treatment arm 1 and the hazard rate and the cumulative hazard in treatment arm 0.
#' The treatment arm 1, will be considered as the "treatment arm" and the treatment arm 0, will be considered as the "control arm".
#' The hazard ratio will be calculated as the hazard rate in the treatment arm vs. the hazard rate in the control arm.
#' The hazard rates and the cumulative hazards will be substitutes in the following formulas depending on the chosen frailty distribution.
#'\itemize{
#' \item{For the Gamma frailty}{ -hr=(haz1/haz0)*exp(theta*(cum.haz1-cum.haz0))}
#' \item{For the Inverse Gaussian frailty }{ -hr=(haz1/haz0)*((1+theta*cum.haz1)/(1+theta*cum.haz0))}
#' }
#'
#' @param frailty.type This function can calculate the causal hazard ratio under two frailty distribution-1=Gamma frailty distribution, 2=Inverse Gaussian distribution. The default value is the Gamma distribution.
#' @param haz1 This is the hazard rate in treatment arm 1 (usually the treatment group)
#' @param haz0 This is the hazard rate in treatment arm 0 (usually the control group)
#' @param cum.haz1 This is the cumulative hazard rate in treatment arm 1 (usually the treatment group)
#' @param cum.haz0 This is the cumulative hazard rate in treatment arm 0 (usually the control group)
#' @param theta This is the frailty variance
#'
#' @return The kernel-based causal hzard ratio estimator 
#' @references 
#' "Sensitivity analysis for the causal hazard ratio in randomized and observational studies" by Rachel Axelrod and Daniel Nevo (arXiv link pending)
#' @export

causal.hazard.frailty.type<-function(frailty.type=1,haz1,haz0,cum.haz1,cum.haz0,theta){
  
  #Gamma frailty
  if(frailty.type==1){
    hr<-(haz1/haz0)*exp(theta*(cum.haz1-cum.haz0))
    
  }
  #IG frailty
  if(frailty.type==2){
    hr<-(haz1/haz0)*((1+theta*cum.haz1)/(1+theta*cum.haz0))
    
  }
  return(hr)
}
