% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta_inverse_gamma.R
\name{theta.inverse.gaussian}
\alias{theta.inverse.gaussian}
\title{Calculation of frailty variance for the Inverse Gaussian frailty distribution}
\usage{
theta.inverse.gaussian(data, tau)
}
\arguments{
\item{data}{- data of Kendall's tau and corresponds theta values}

\item{tau}{-the wanted tau value}
}
\value{
Return the corresponds frailty variance, theta
}
\description{
This function calculates the frailty variance for any given value of Kendall's tau between [0,0.5)
}
\details{
The possible values of Kendall's tau under the inverse Gaussian distribution are [0,0.5). 
Since the function that connects between the Kendall's tau and theta is not invertable, a file with two
column : tau and the corresponds theta values needs to be specified.
}
