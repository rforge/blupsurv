\name{summary.unirec}
\alias{summary.unirec}
\alias{print.summary.unirec}
\title{Summary method for unirec objects}
\description{Prints a summary of a fit from \code{\link{unirec}}, including
    pretty-printed p-values and significance indicators.}
\usage{
\method{summary}{unirec}(object,digits=4,...)
}
\arguments{
    \item{object}{an object of type \code{\link{unirec}}.}
    \item{digits}{number of digits to be used in pretty-printing output.}
    \item{...}{additional parameters for \code{\link{print}}.}
}
\value{
    An item of type \code{summary.unirec} with components \code{call, summary.reg,
    summary.disp}
    as documented in \code{\link{unirec}}.
}
\seealso{
\code{\link{unirec}}
}
\keyword{survival}