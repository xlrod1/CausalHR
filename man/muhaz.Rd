\name{muhaz}
\alias{muhaz}
\title{
Estimate hazard function from right-censored data.
}
\description{
Estimates the hazard function from right-censored data using kernel-based
methods. Options include three types of bandwidth functions, three types of
boundary correction, and four shapes for the kernel function. Uses the
global and local bandwidth selection algorithms and the boundary kernel
formulations described in Mueller and Wang (1994). The nearest neighbor
bandwidth formulation is based on that described in Gefeller and Dette
(1992). The statistical properties of many of these estimators are reported
and compared in Hess et al (1999). Based on the HADES program developed by
H.G. Mueller. Returns an object of class 'muhaz.'
NOTE: For comparison to the smoothed hazard function estimates,
we have also made available a set of functions based on piecewise
exponential estimation.  These estimates are similar in concept to
the histogram estimator of the density function.  They give a feel for
the features of the data without the manipulations involved in
smoothing. They also help to confirm that muhaz is generating
realistic estimates of the underlying hazard function.  These
functions are called: pehaz, plot.pehaz, lines.pehaz, print.pehaz.
}
\usage{
muhaz(times, delta, subset, min.time, max.time, bw.grid, bw.pilot,
      bw.smooth, bw.method="local", b.cor="both", n.min.grid=51,
      n.est.grid=101, kern="epanechnikov")
}
\arguments{
\item{times}{
Vector of survival times.  Does not need to be sorted.
}
\item{delta}{
Vector indicating censoring
      0 - censored (alive)
      1 - uncensored (dead)
If delta is missing, all the observations are assumed uncensored.
}
\item{subset}{
Logical vector, indicating the observations used in analysis.
      T - observation is used
      F - observation is not used
If missing, all the observations will be used.
}
\item{min.time}{
Left bound of the time domain used in analysis.
If missing, min.time is considered 0.
}
\item{max.time}{
Right bound of the time domain used in analysis.
If missing, max.time is the time at which ten patients remain at risk.
}
\item{bw.grid}{
Bandwidth grid used in the MSE minimization.
If bw.method="global" and bw.grid has one component only,
no MSE minimization is performed.  The hazard estimates are computed
for the value of bw.grid.
If bw.grid is missing, then a bandwidth grid of 21 components is built,
having as bounds:
      \deqn{[0.2*bw.pilot, 20*bw.pilot]}
}
\item{bw.pilot}{
Pilot bandwidth used in the MSE minimization.
If missing, the default value is the one recommended by Mueller and
Wang (1994):
      \deqn{bw.pilot = (max.time-min.time) / (8*nz^0.2)}
where nz is the number of uncensored observations
}
\item{bw.smooth}{
Bandwidth used in smoothing the local bandwidths.
Not used if bw.method="global"
If missing:
      \deqn{bw.smooth = 5 * bw.pilot}
}
\item{bw.method}{
Algorithm to be used.  Possible values are:
"global" - same bandwidth for all grid points.  The optimal
bandwidth is obtained by minimizing the IMSE.
"local"  - different bandwidths at each grid point.  The optimal
bandwidth at a grid point is obtained by minimizing the local MSE.
"knn"    - k nearest neighbors distance bandwidth.  The optimal number of
neighbors is obtained by minimizing the IMSE.
Default value is "local".
Only the first letter needs to be given (e.g. "g", instead of "global").
}
\item{b.cor}{
Boundary correction type.  Possible values are:
"none" - no boundary correction
"left" - left only correction
"both" - left and right corrections
Default value is "both". Only the first letter needs to be given
(e.g. b.cor="n").
}
\item{n.min.grid}{
Number of points in the minimization grid.  This value greatly influences
the computing time.
Default value is 51.
}
\item{n.est.grid}{
Number of points in the estimation grid, where hazard estimates are
computed.
Default value is 101.
}
\item{kern}{
Boundary kernel function to be used.  Possible values are:
"rectangle", "epanechnikov", "biquadratic", "triquadratic".
Default value is "epanechnikov". Only the first letter needs to be given
(e.g. kern="b").
}}
\value{
Returns an object of class 'muhaz', containing input and output values.
Methods working on such an object are: plot, lines, summary.  For a
detailed description of its components, see object.muhaz.
}
\details{
The muhaz object contains a list of the input data and parameter values as
well as a variety of output data. The hazard function estimate is contained
in the haz.est element and the corresponding time points are in est.grid.
The unsmoothed local bandwidths are in bw.loc and the smoothed local
bandwidths are in bw.loc.sm.


For bw.method = 'local' or 'knn', to check the shape of the bandwidth
function used in the estimation, use
\code{plot(fit$pin$min.grid, fit$bw.loc)} to
plot the unsmoothed bandwidths and use
\code{lines(fit$est.grid, fit$bw.loc.sm)}
to superimpose the smoothed bandwidth function. Use bw.smooth to change the
amount of smoothing used on the bandwidth function.


For bw.method='global', to check the minimization process, plot the
estimated IMSE values over the bandwidth search grid. Use
\code{plot(fit$bw.grid, fit$globlmse)}.
Use k.grid and k.imse for bw.method='k'. You may want to
repeat the search using a finer grid over a shorter interval to fine-tune
the optimization or if the observed minimum is at the extreme of the grid
you should specify a different grid.


}
\references{
1. H.G. Mueller and J.L. Wang - Hazard Rates Estimation Under Random Censoring
with Varying Kernels and Bandwidths, Biometrics 50, 61-76, March 1994


2. O. Gefeller and H. Dette - Nearest Neighbour Kernel Estimation of the
Hazard Function From Censored Data, J. Statist. Comput. Simul., Vol.43, 1992,
93-101


3. K.R. Hess, D.M. Serachitopol and B.W. Brown - Hazard Function Estimators: A Simulation Study, Statistics in Medicine (in press).


}
\seealso{
summary.muhaz, plot.muhaz, lines.muhaz, muhaz.object
}
\examples{
# to compute a locally optimal estimate
data(cancer, package="survival")
attach(ovarian)
fit1 <- muhaz(futime, fustat)
plot(fit1)
summary(fit1)
# to compute a globally optimal estimate
fit2 <- muhaz(futime, fustat, bw.method="g")
# to compute an estimate with global bandwidth set to 5
fit3 <- muhaz(futime, fustat, bw.method="g", bw.grid=5)
}
\keyword{survival}
\keyword{smooth}
