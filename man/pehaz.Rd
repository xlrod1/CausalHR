\name{pehaz}
\alias{pehaz}
\title{
Estimates piecewise exponential hazard function from
right-censored data.
}
\description{
Divides the time domain into bins of equal width, and then
estimates the hazard in each bin as the number of events in
that bin divided by the total follow-up time in that bin.
}
\usage{
pehaz(times, delta=NA, width=NA, min.time=0, max.time=NA)
}
\arguments{
\item{times}{
Vector of survival times.  Does not need to be sorted.
}
\item{delta}{
Vector indicating censoring
0 - censored (alive)
1 - uncensored (dead)
If status is missing, all observations are assumed
uncensored.
}
\item{width}{
Bin width.  Default value is that recommended by
Mueller, \eqn{width = (max.time-min.time) / (8*(nu)^0.2)} where
nu is the number of uncensored observations.
}
\item{min.time}{
Left bound of the time domain used in the
analysis.  If missing, min.time is considered 0.
}
\item{max.time}{
Right bound of the time domain used in the
analysis.  If missing, max.time is considered
max(times).
}}
\value{
Returns  an  object of class 'pehaz', containing input and
output values.  Methods working on such an object are:
plot, lines, print.  For a detailed description of its
components, see \code{object.pehaz}.
}
\seealso{
\code{\link{pehaz.object}}
}
\examples{
data(cancer, package="survival")
attach(ovarian)
fit <- pehaz(futime, fustat)
plot(fit)
}
\keyword{survival}

