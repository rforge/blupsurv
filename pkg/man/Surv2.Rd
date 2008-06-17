\name{Surv2}
\alias{Surv2}
\title{ Bivariate survival object }
\description{
 Creates a survival object for bivariate data, to be used in \code{\link{bivrec}}.
}
\usage{
Surv2(start, stop, status1, status2)
}
\arguments{
  \item{start}{ starting time for the interval. }
  \item{stop}{ ending time for the interval. }
  \item{status1}{ status indicator for the first event process, 1=event, 0=no event. }
  \item{status2}{ analogous to \code{status1}, for the second event process. }
}
\value{
    an object of class \code{Surv2}, implemented as a data frame of 4 columns.
}
\seealso{ \code{\link{Surv}}, \code{\link{bivrec}}}
\keyword{ survival }
