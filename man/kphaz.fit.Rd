\name{kphaz.fit}
\alias{kphaz.fit}
\title{
kphaz.fit
}
\description{
Calculates Kaplan-Meier type hazard estimates.
}
\usage{
kphaz.fit(time,status,strata,q=1,method="nelson")
}
\arguments{
\item{time}{
A vector of time values; all values must be greater than
or equal to zero. Missing values (NAs) are not allowed.
}
\item{status}{
A vector of status values.  The values are 0 for
censored or 1 for uncensored (dead). Missing values
(NAs) are not allowed. Must have the same length as time.
}
\item{strata}{
An optional vector that will be used to divide the
subjects into disjoint groups. Each group generates a
hazard curve. If missing, all subjects are assumed to be in the
same strata. Missing values (NAs) are allowed.
}
\item{q}{
Number of failure times combined for estimatingthe hazard at their 
midpoint. Default is 1.
}
\item{method}{
Type of hazard estimation made. Must be one of "nelson"
or "product-limit". The default is "nelson".
}
}
\value{
A list representing the results of the hazard estimation,
with the following components:

\item{time}{
A vector containing the times at which hazard estimations were made.
}
\item{haz}{
A vector containing the hazard estimate at each time.
}
\item{var}{
A vector containing variance estimates for each hazard estimate.
}
\item{strata}{
A vector which divides the hazard estimate into
disjoint groups. This vector is returned only if 'strata' 
is defined when 'kphaz.fit' is called.
}}
\details{
Let \deqn{t[1] < t[2] < \cdots < t[m]} denote the m "distinct" death times.


1. Estimate the cumulative hazard, H[t[j]], and the variance of
the cumulative hazard, Var(H[t[j]]), at each of the m distinct
death times according to the method selected.


a. For the "nelson" method:
   \deqn{H[t[j]] = sum(t[i] <= t[j]) status[i]/(n-i+1)}
   \deqn{Var(H[t[j]]) = \sum(t[i] <= t[j]) status[i]/((n-i+1)^2)}


b. For the "product-limit" metod:
   \deqn{H[t[j]] = sum(t[i] <= t[j]) -log(1 - status[i]/(n-i+1))}
   \deqn{Var(H[t[j]]) = sum(t[i] <= t[j]) 
                     status[i]/((n-i+1)*(n-i))}


2. For k=1,\dots{},(m-q), define the hazard estimate and variance at 
time[k] = (t[q+j]+t[j])/2 to be
   \deqn{haz[time[k]] = (H[t[q+j]]-H[t[j]])/(t[q+j]-t[j])}
   \deqn{var[time[k]] = (Var(H[t[q+j]])-Var(H[t[j]]))/
                     (t[q+j]-t[j])^2}


Note that if the final time is a death time rather than a censoring
time, the "product-limit" estimate will be Inf for the final hazard
and variance estimates.
}
\references{
Jarjoura, David (1988).
Smoothing Hazard Rates with Cubic Splines.
\emph{Commun. Statist. -Simula.}
\bold{17(2)}, 377-392.
}
\seealso{
\code{\link{kphaz.plot}}
}
\examples{
time <- 1:10
status <- rep(1,10)
kphaz.fit(time,status)
}
\keyword{survival}

