#' Calculation of the kernel-based and the Cox-based estimators for the causal hazard ratio
#' 
#' 
#' This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for a given value of tau
#' 
#' @details 
#'  This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for a given value of tau. This function can be used for RCT ot for obaservational studies
#'  For the observational studies setting, a vector name "w",of the observations' weights needs to pre-specified and need to be of the columns' .
#'  The weights are used also in the kernel-based estimator and ahe Cox-based estimator.
#'  This function calls \code{\link{theta.gamma}} or\code{\link{theta.gamma}} to calculate the frailty variance, theta. When calculating 
#'  the causal hazard ratio under the inverse Gaussian, a data of tau needs to be specified. A tau data can be found in the first author Githab account.
#'  This function also calls  \code{\link{muhaz_hr_for_real_data_example}} to calculate the cox-based estimator and the 
#'  \code{\link{coxph}} and \code{\link{base.haz.estimator}} ( from the survival package to calculate the Cox-based estimator.
#'  If a calculation of standard errors  and 95\% pointwise percentile CIs  are wanted, please use the \code{\link{boot_ci_real_data}} function.
#'  
#'  
#'   
#'
#'
#' @param tau The wanted Kendall's tau. Under the Gamma frailty the possible values of Kendall's tau are [0,1] and under the inverse Gaussian [0,0.5)
#' @param data The data that the estimators will be estimated on. In case of RCT, a data with column named "time, treatment, status" needs to be specified. In case of observational studies a firth column named "W" that contains the weights need to be specified also. 
#' @param frailty.type This function can calculate the causal hazard ratio under two frailty distribution-1=Gamma frailty distribution, 2=Inverse Gaussian distribution. The default value is the Gamma distribution.
#' @param confounder "no"- for RCT setting, or in case a non-weighted estimator wanted to be calculated. "yes"-for observational studies setting, or in evey other case a weighted estimator need to be calculated.
#' @param n.est.grid Number of points in the estimation grid, where hazard estimates are computed. Default value is 101.
#' @param b.cor Boundary correction type. Possible values are: "none" - no boundary correction "left" - left only correction "both" - left and right corrections Default value is "both". Only the first letter needs to be given (e.g. b.cor="n").
#' @param end.time Right bound of the time domain used in analysis.
#' @param min.time Left bound of the time domain used in analysis. If missing, min.time is considered 0.
#' @param kern Boundary kernel function to be used. Possible values are: "rectangle", "epanechnikov", "biquadratic", "triquadratic". Default value is "epanechnikov". Only the first letter needs to be given (e.g. kern="b").
#' @param data.tau a csv file of the taus file, need to be specified in case a causal hazard ratio under the inverse Gaussian need to be calculated
#'
#' @return Returns dataframe containing the following columns:time=the estimation grid points, HR=the kernel-based estimator, HR.sp=the Cox-based estimator , "tau"=the selected tau value 
#' 

CausalHR.without.bootstrap<-function(tau,data,frailty.type,confounder,n.est.grid=101,b.cor,end.time,min.time=0,kern,data.tau=NA){
  if (frailty.type==1){
    theta01<-theta.gamma(tau=tau)
  }
  
  if (frailty.type!=1){
    tau.theta<-theta.inverse.gaussian(data.tau,tau=tau)
    theta01<-tau.theta$theta
    
  }
  
  bb.sample<-muhaz_hr_for_real_data_example(frailty.type,data,theta01, confounder,n.est.grid =n.est.grid,b.cor=b.cor,end.time=end.time,min.time =min.time,kern=kern)
  muhaz.hr.df.list<-bb.sample$muhaz.hr.df
  muhaz.hr.df.list$tau<-as.character(tau)
  
  if(confounder=="yes"){
    cox.reg<-coxph(Surv(time,status)~treatment,data=data,weights = W)
  }
  if(confounder=="no"){
    cox.reg<-coxph(Surv(time,status)~treatment,data=data)
  }
  beta.hat<-as.numeric(cox.reg$coefficients[1])
  time<-muhaz.hr.df.list$time
  base.haz.estimator<-stepfun(x=basehaz(cox.reg,centered = FALSE)$time,y=c(0,basehaz(cox.reg,centered = FALSE)$hazard))
  base.line.cumsum<- base.haz.estimator(time) 
  
  if(frailty.type==1){
    HR.sp<-exp(beta.hat)*exp(theta01*base.line.cumsum*(exp(beta.hat)-1))
  }
  if(frailty.type!=1){
    
    HR.sp<-exp(beta.hat)*((1+theta01*exp(beta.hat)*base.line.cumsum)/(1+theta01*base.line.cumsum))
  }
  muhaz.hr.df.list$HR.sp<-HR.sp
  
  
  
  muhaz.hr.df.list
}