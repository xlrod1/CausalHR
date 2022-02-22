\name{pehaz.object}
\alias{pehaz.object}
\title{
Estimated Piecewise Exponential Hazard Rate Object
}
\description{
This class of objects is returned by \code{pehaz}, which estimates the
hazard function from censored data.
}
\section{METHODS}{
Objects of this class have methods for the functions \code{summary}, 
\code{plot}, and \code{lines}.
}
\section{STRUCTURE}{
Common components of a \code{pehaz} object:
\describe{
  A list containing the following components:
\item{call}{the call to pehaz}
\item{Width}{the width of the bins}
\item{Cuts}{the cutpoints used for the bins}
\item{Hazard}{the estimated hazard for each bin}
\item{Events}{the number of events in each bin}
\item{At.Risk}{the number at risk in each bin}
\item{F.U.Time}{the followup time (for the bin?)}
}
}

\seealso{ \code{\link{pehaz}}}
\keyword{survival}

