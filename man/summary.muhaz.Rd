\name{summary.muhaz}
\alias{summary.muhaz}
\title{Display the most important input parameters used in calling the
  `muhaz' function.}
\description{It also displays some of the output data.
Common to all three methods:
* number of observations
* number of censored observations
* bandwidth method used (global, local or nearest neighbor)
* boundary correction type (none, left only, both left and right)
* kernel type (rectangle, Epanechnikov, biquadradic, triquadratic)
* minimum time
* maximum time
* number of points in MSE minimization grid
* number of points in estimation grid
* pilot bandwidth
* estimated IMSE for optimal bandwidth
For bw.method="global" also reports optimal global bandwidth.
For bw.method="knn" also reports optimal number of nearest neighbors.
For bw.method="local" and bw.method="knn" also reports smoothing bandwidth
used to smooth the optimal local bandwidths.
}
\usage{
\method{summary}{muhaz}(object, ...)
}
\arguments{
  \item{object}{
    Object of class \code{muhaz} (output from calling \code{muhaz}
    function)}
  \item{...}{Ignored for now.}
}
\seealso{
\code{\link{muhaz}}, \code{\link{muhaz.object}}
}
\examples{
data(cancer, package="survival")
attach(ovarian)
fit <- muhaz(futime, fustat)
summary(fit)
}
\keyword{survival}

