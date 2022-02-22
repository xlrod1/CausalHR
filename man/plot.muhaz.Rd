\name{plot.muhaz}
\alias{plot.muhaz}
\alias{lines.muhaz}
\title{Plots estimated hazard function from an object of class muhaz.}
\description{
Default time limits are those provided to \code{muhaz}, which default
to zero and the time corresponding to when ten patients remain
at risk.  Default y-axis limits are 0 and the maximum estimated
hazard rate.
Additional lines can be added to the same set of axes using the \code{lines}
method.
}
\usage{
\method{plot}{muhaz}(x, ylim, type, xlab, ylab, ...)
\method{lines}{muhaz}(x, ...)
}
\arguments{
  \item{x}{
    Object of class \code{muhaz} (output from calling \code{muhaz}
    function)}
  \item{ylim}{Limits for the y axis.}
  \item{type}{type argument for \code{plot}.}
  \item{xlab}{Label for the x axis.}
  \item{ylab}{Label for the y axis.}
  \item{...}{Additional arguments to be passed along.}
}
\seealso{
\code{\link{muhaz}}
}
\keyword{survival}
