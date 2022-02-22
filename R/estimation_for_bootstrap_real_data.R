#' The wrapper function of the bootstrap estimation
#' 
#' 
#' Calculation of the kernel-based and the Cox-based estimators for the causal hazard ratio for each bootstrap sample
#' 
#' @details 
#'  This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for each bootstrap sample. T
#'  This function can be used for RCT or for observational studies. For the observational studies setting, a vector name "w",of the observations' weights needs to pre-specified and need to be of the columns' .
#'  The weights are used also in the kernel-based estimator and ahe Cox-based estimator.
#'  This function calls \code{\link{muhaz_hr_for_real_data_example}} function that estimate the kernel-based casual hazard ratio in each sample,  \code{\link{coxph}} and \code{\link{base.haz.estimator}} ( from the survival package) to calculate the Cox-based estimator.
#'  
#'  
#' @param data The data that the estimators will be estimated on. In case of RCT, a data with column named "time, treatment, status" needs to be specified. In case of observational studies a firth column named "W" that contains the weights need to be specified also. 
#' @param indices allows boot to select sample, dont need to be sppecified
#' @param frailty.type This function can calculate the causal hazard ratio under two frailty distribution-1=Gamma frailty distribution, 2=Inverse Gaussian distribution. The default value is the Gamma distribution.
#' @param tau The wanted Kendall's tau. Under the Gamma frailty the possible values of Kendall's tau are [0,1] and under the inverse Gaussian [0,0.5)
#' @param confounder "no"- for RCT setting, or in case a non-weighted estimator wanted to be calculated. "yes"-for observational studies setting, or in evey other case a weighted estimator need to be calculated. The defoult option is "no".
#' @param n.est.grid Number of points in the estimation grid, where hazard estimates are computed. Default value is 101.
#' @param b.cor Boundary correction type. Possible values are: "none" - no boundary correction "left" - left only correction "both" - left and right corrections Default value is "both". Only the first letter needs to be given (e.g. b.cor="n").
#' @param min.time Left bound of the time domain used in analysis. If missing, min.time is considered 0.
#' @param end.time Right bound of the time domain used in analysis.
#' @param kern Boundary kernel function to be used. Possible values are: "rectangle", "epanechnikov", "biquadratic", "triquadratic". Default value is "epanechnikov". Only the first letter needs to be given (e.g. kern="b").
#' @param data.tau a csv file of the taus file, need to be specified in case a causal hazard ratio under the inverse Gaussian need to be calculated
#'
#' @return Returns a vector the bootstrap sample causal hazard ratio estimator.

estimation.for.bootsrap.real.data<-function(data, indices,frailty.type,tau, confounder,n.est.grid ,b.cor,end.time,min.time,kern,data.tau=NA){
  
  d <- data[indices,] # allows boot to select sample
  
  if (frailty.type==1){
    theta<-theta.gamma(tau=tau)
  }
  
  if (frailty.type!=1){
    tau.theta<-theta.inverse.gaussian(data.tau,tau=tau)
    theta<-tau.theta$theta
    
  }
 
  
  if(confounder=="no"){
    cox.reg<-coxph(Surv(time,status)~treatment,data=d)
  }
  if(confounder=="yes"){
    cox.reg<-coxph(Surv(time,status)~treatment,data=d,weights = W)
  }
  #kernel
  fit <-muhaz_hr_for_real_data_example(frailty.type=frailty.type,d,theta=theta, confounder=confounder,
                                       n.est.grid = n.est.grid ,b.cor=b.cor,end.time=end.time,min.time =min.time,kern=kern)
  
  fit<-fit$muhaz.hr.df
  
  #cox
  time<-fit$time
  beta.hat<-as.numeric(cox.reg$coefficients[1])
  base.haz.estimator<-stepfun(x=basehaz(cox.reg,centered = FALSE)$time,y=c(0,basehaz(cox.reg,centered = FALSE)$hazard))
  base.line.cumsum<- base.haz.estimator(time) 

  if (frailty.type==1){
    HR.sp<-exp(beta.hat)*exp(theta*base.line.cumsum*(exp(beta.hat)-1))
  }
  
  if (frailty.type!=1){
    HR.sp<-exp(beta.hat)*((1+theta*exp(beta.hat)*base.line.cumsum)/(1+theta*base.line.cumsum))
    
  }
  
  
 
  return(c(fit$HR,HR.sp))
}

