\name{summary.bivrec}
\alias{summary.bivrec}
\alias{print.summary.bivrec}
\title{Summary method for bivrec objects}
\description{Prints a summary of a fit from \code{\link{bivrec}}, including
    pretty-printed p-values and significance indicators.}
\usage{
\method{summary}{bivrec}(object,digits=4,...)
}
\arguments{
    \item{object}{an object of type \code{\link{bivrec}}.}
    \item{digits}{number of digits to be used in pretty-printing output.}
    \item{...}{additional parameters for \code{\link{print}}.}
}
\value{
    An item of type \code{summary.bivrec} with components \code{call, summary.reg,
    summary.disp}
    as documented in \code{\link{bivrec}}.
}
\seealso{
\code{\link{bivrec}}
}
\keyword{survival}