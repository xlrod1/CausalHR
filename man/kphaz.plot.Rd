\name{kphaz.plot}
\alias{kphaz.plot}
\title{
kphaz.plot
}
\description{
Plots a Kaplan-Meier type hazard estimate.
}
\usage{
kphaz.plot(fit, ...)
}
\arguments{
\item{fit}{
A list representing the results of a call to "kphaz.fit".
}
\item{...}{
Any legal argument for the plot function.
}}
\section{Side Effects}{
A plot with multiple hazard curves. One for each unique
strata with 1 or more point.
}
\seealso{
\code{\link{kphaz.fit}}
}
\examples{
# Use "kphaz.fit" to generate a hazard estimate
data(cancer, package="survival")
attach(ovarian)
kpfit <- kphaz.fit(futime, fustat)
# Use "kphaz.plot" to plot the estimate
kphaz.plot(kpfit)
}
\keyword{survival}
