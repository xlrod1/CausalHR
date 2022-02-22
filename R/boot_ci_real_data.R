#' Calculation of the kernel-based and the Cox-based estimators 
#' 
#' 
#' Calculation of the kernel-based and the Cox-based estimators for the causal hazard ratio along with standard errors and 95\% pointwise percentile bootstrap CIs
#' 
#' @details 
#'  This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for a given value of tau. This function can be used for RCT ot for observational studies
#'  For the observational studies setting, a vector name "w",of the observations' weights needs to pre-specified and need to be of the columns' .
#'  The weights are used also in the kernel-based estimator and ahe Cox-based estimator.
#'  This function calls \code{\link{boot}} function that estimate the csual hazard ratio in each smple using the function \code{\link{estimation.for.bootsrap.real.data}} to preform bootstrap sampling. Then it calculates the empirical standard errors and the 95\% pointwise CIs by the percentile method.
#'  
#'  
#' @param data The data that the estimators will be estimated on. In case of RCT, a data with column named "time, treatment, status" needs to be specified. In case of observational studies a firth column named "W" that contains the weights need to be specified also. 
#' @param R The number of bootstrap replicates. Usually this will be a single positive integer. For importance resampling, some resamples may use one set of weights and others use a different set of weights. In this case R would be a vector of integers where each component gives the number of resamples from each of the rows of weights.Defoult value is 300.
#' @param frailty.type This function can calculate the causal hazard ratio under two frailty distribution-1=Gamma frailty distribution, 2=Inverse Gaussian distribution. The default value is the Gamma distribution.
#' @param confounder "no"- for RCT setting, or in case a non-weighted estimator wanted to be calculated. "yes"-for observational studies setting, or in evey other case a weighted estimator need to be calculated. The default option is "no".
#' @param tau The wanted Kendall's tau. Under the Gamma frailty the possible values of Kendall's tau are [0,1] and under the inverse Gaussian [0,0.5)
#' @param min.time Left bound of the time domain used in analysis. If missing, min.time is considered 0.
#' @param end.time Right bound of the time domain used in analysis.
#' @param n.est.grid Number of points in the estimation grid, where hazard estimates are computed. Default value is 101.
#' @param b.cor Boundary correction type. Possible values are: "none" - no boundary correction "left" - left only correction "both" - left and right corrections Default value is "both". Only the first letter needs to be given (e.g. b.cor="n").
#' @param kern Boundary kernel function to be used. Possible values are: "rectangle", "epanechnikov", "biquadratic", "triquadratic". Default value is "epanechnikov". Only the first letter needs to be given (e.g. kern="b").
#' @param max.HR a correction factor for the bootstrap results. if the bootstrap result is max.HR times bigger than the dataset result than the result is NA.
#' @param data.tau a csv file of the taus file, need to be specified in case a causal hazard ratio under the inverse Gaussian need to be calculated
#'
#' @return Returns a list of two dataframes. Each dataframe contains the following columns are :A dataframe with the kernel based estimator results. The columns are :time=the estimation grid points,
#'  HR= The causal hazard ratio estimator , "tau"=the selected tau value ,"dist"-the selected frailty distribution ,"method"=the estimation method,
#'  "BSE"=the emprical standard erorrs , "X025"=the lower bound of the 95% CI, "X095"-the upper bound of the 95% CI
#' \itemize{
#' \item{CI.Kernel}{ The results for the kernel-based estimator}
#' \item{{CI.Cox}{ The results for the Cox-based estimator}
#' }
#'
#' @export
#' "Sensitivity analysis for the causal hazard ratio in randomized and observational studies" by Rachel Axelrod and Daniel Nevo (arXiv link pending)
#' @examples
boot.ci.tau<-function(data,R=300,frailty.type=1,confounder="no",tau,min.time=0,end.time,n.est.grid=101,kern="e",b.cor="b",max.HR=5,data.tau=NA){

  
  results<- boot(data=data, statistic=estimation.for.bootsrap.real.data,R=R, 
                 frailty.type=frailty.type,tau=tau,confounder=confounder,
                 n.est.grid=n.est.grid ,b.cor="b",
                 end.time=end.time,min.time=min.time,kern=kern,data.tau=data.tau)  
  
  results$t<-apply(results$t,c(1,2),function(x){ifelse(x%in%c(Inf,-Inf),NA,x)})
  t0<-results$t0
  #yy<-results$t
  #yyy<-results$t
  for( j in 1:ncol(results$t)){
    for( i in 1:nrow(results$t)){
      if(!is.na( results$t[i,j])&&results$t[i,j]>=max.HR*t0[j]){ results$t[i,j]=NA}
    }
    
  }
  
  
  #CI kernel 
  ci.Kernel<-data.frame(t(apply(results$t[(1:R),(1:n.est.grid)], 2, function(u){quantile(u, probs = c(.025, .975), type = 6,na.rm = TRUE)})))
  ci.Kernel$BSE<-apply(results$t[(1:R),(1:n.est.grid)], 2, function(u){sqrt(var(u,na.rm = TRUE))})
  ci.Kernel$HR<-t0[(1:n.est.grid)]
  
  #CI cox
 ci.Cox<-data.frame(t(apply(results$t[(1:R),((n.est.grid+1):(2*n.est.grid))], 2, function(u){quantile(u, probs = c(.025, .975), type = 6,na.rm = TRUE)})))
 ci.Cox$BSE<-apply(results$t[(1:R),((n.est.grid+1):(2*n.est.grid))], 2, function(u){sqrt(var(u,na.rm = TRUE))})
 ci.Cox$HR<-t0[((n.est.grid+1):(2*n.est.grid))]


  
 ci.Cox$time<- ci.Kernel$time<-  time.vec<-seq(from=min.time,to = end.time,length.out = 51)
  
  if(frailty.type==1){ci.Cox$dist<- ci.Kernel$dist<-"Gamma"}
  else{ci.Cox$dist<- ci.Kernel$dist<-"IG"}
  
 ci.Cox$tau<- ci.Kernel$tau<-as.character(tau)
 
  # ci.Kernel$HR<-hr.data$HR
  # ci.Cox$HR<-hr.data$HR.sp
 ci.Cox$method<-"Cox"
 ci.Kernel$method="Kernel"
  #write.csv( ci,paste("ci_","ft_",frailty.type,"_tau_",tau,sep = "",".csv"),row.names = F)
  return(list(CI.Kernel=ci.Kernel,CI.Cox=ci.Cox))
}
