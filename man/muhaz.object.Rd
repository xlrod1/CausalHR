\name{muhaz.object}
\alias{muhaz.object}
\title{
Estimated Hazard Rate Object
}
\description{
This class of objects is returned by \code{muhaz}, which estimates the
hazard function from censored data.
}
\section{METHODS}{
Objects of this class have methods for the functions \code{summary}, 
\code{plot}, and \code{lines}.
}
\section{STRUCTURE}{
Common components of a \code{muhaz} object:
}
\arguments{
\item{pin.common}{
list containing the input parameters.  Its components are:
* times - the survival times vector
* delta - the censoring vector
* nobs  - the number of observations
* min.time - the minimum time used in analysis
* max.time - the maximum time used in analysis
* n.min.grid - number of points in the minimization grid
* min.grid - the minimization grid
* n.est.grid - number of points in the estimation grid
* bw.pilot - the pilot bandwidth
* bw.smooth - the smoothing bandwidth for the local optimal bandwidths
* bw.method - the method used to estimate the hazard rates
* b.cor - boundary correction used
* kernel.type - kernel function used in the smoothing
}
\item{est.grid}{
the estimation grid where the hazard rates are computed.
}
\item{haz.est}{
the hazard estimates calculated at the estimation grid points.
}
\item{imse.opt}{
IMSE for the optimal bandwidth.
}
\item{bw.glob}{
optimal global bandwidth.  For bw.method="global" only.
}
\item{glob.imse}{
vector of IMSE, computed at each point in bw.grid.
For bw.method="global" only.
}
\item{bw.grid}{
Grid of bandwidth values used in the minimization.
For bw.method="global" or "local".
}
\item{bw.loc}{
vector of optimal local bandwidths computed by minimizing the MSE at each
point in the minimization grid. Not used for bw.method="local".
}
\item{bw.loc.sm}{
vector of smoothed local bandwidths, computed at each point in the estimation
grid by smoothing \code{bw.loc} using \code{bw.smooth} as the smoothing bandwidth.
Not used for bw.method="global".
}
\item{bias.min}{
vector of minimized bias, computed at each point in the minimization grid
for the optimal local bandwidth.  For bw.method="local" only.
}
\item{var.min}{
vector of minimized variance, computed at each point in the minimization grid
for the optimal local bandwidth.  For bw.method="local" only.
}
\item{k.grid}{
grid of nearest neighbor numbers used in the minimization.
For bw.method="knn" only.
}
\item{k.imse}{
vector of IMSE, computed at each of the points of \code{k.grid}.
For bw.method="knn" only.
}
\item{k.opt}{
optimum number of nearest neighbors.  For bw.method="knn" only.
}}
\seealso{
\code{\link{muhaz}}, \code{\link{plot.muhaz}}, \code{\link{summary.muhaz}}.
}
\keyword{survival}
% Converted by Sd2Rd version 1.10.
