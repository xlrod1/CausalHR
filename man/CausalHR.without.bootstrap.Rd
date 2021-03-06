% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/muhaz_by_tau.R
\name{CausalHR.without.bootstrap}
\alias{CausalHR.without.bootstrap}
\title{Calculation of the kernel-based and the Cox-based estimators for the causal hazard ratio}
\usage{
CausalHR.without.bootstrap(
  tau,
  data,
  frailty.type,
  confounder,
  n.est.grid = 101,
  b.cor,
  end.time,
  min.time = 0,
  kern,
  data.tau = NA
)
}
\arguments{
\item{tau}{The wanted Kendall's tau. Under the Gamma frailty the possible values of Kendall's tau are [0,1] and under the inverse Gaussian [0,0.5)}

\item{data}{The data that the estimators will be estimated on. In case of RCT, a data with column named "time, treatment, status" needs to be specified. In case of observational studies a firth column named "W" that contains the weights need to be specified also.}

\item{frailty.type}{This function can calculate the causal hazard ratio under two frailty distribution-1=Gamma frailty distribution, 2=Inverse Gaussian distribution. The default value is the Gamma distribution.}

\item{confounder}{"no"- for RCT setting, or in case a non-weighted estimator wanted to be calculated. "yes"-for observational studies setting, or in evey other case a weighted estimator need to be calculated.}

\item{n.est.grid}{Number of points in the estimation grid, where hazard estimates are computed. Default value is 101.}

\item{b.cor}{Boundary correction type. Possible values are: "none" - no boundary correction "left" - left only correction "both" - left and right corrections Default value is "both". Only the first letter needs to be given (e.g. b.cor="n").}

\item{end.time}{Right bound of the time domain used in analysis.}

\item{min.time}{Left bound of the time domain used in analysis. If missing, min.time is considered 0.}

\item{kern}{Boundary kernel function to be used. Possible values are: "rectangle", "epanechnikov", "biquadratic", "triquadratic". Default value is "epanechnikov". Only the first letter needs to be given (e.g. kern="b").}

\item{data.tau}{a csv file of the taus file, need to be specified in case a causal hazard ratio under the inverse Gaussian need to be calculated}
}
\value{
Returns dataframe containing the following columns:time=the estimation grid points, HR=the kernel-based estimator, HR.sp=the Cox-based estimator , "tau"=the selected tau value
}
\description{
This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for a given value of tau
}
\details{
This function calculate the kernel-based and the Cox-based estimators for the causal hazard ratio for a given value of tau. This function can be used for RCT ot for obaservational studies
 For the observational studies setting, a vector name "w",of the observations' weights needs to pre-specified and need to be of the columns' .
 The weights are used also in the kernel-based estimator and ahe Cox-based estimator.
 This function calls \code{\link{theta.gamma}} or\code{\link{theta.gamma}} to calculate the frailty variance, theta. When calculating 
 the causal hazard ratio under the inverse Gaussian, a data of tau needs to be specified. A tau data can be found in the first author Githab account.
 This function also calls  \code{\link{muhaz_hr_for_real_data_example}} to calculate the cox-based estimator and the 
 \code{\link{coxph}} and \code{\link{base.haz.estimator}} ( from the survival package to calculate the Cox-based estimator.
 If a calculation of standard errors  and 95\% pointwise percentile CIs  are wanted, please use the \code{\link{boot_ci_real_data}} function.
}
