\name{plot.bivrec}
\alias{plot.bivrec}
\alias{plot.unirec}
\title{Plot of survivor function for bivrec or unirec objects}
\description{Function to plot the the survivor function
    estimated by \code{\link{bivrec}} or \code{\link{unirec}}.
} 
\usage{
\method{plot}{bivrec}(x, which = c(0,1,2), main=NULL,
    xscale = 1, hazscale = 1, add = FALSE, legend = NULL, ...)
    
\method{plot}{unirec}(x, main=NULL, xscale = 1, hazscale = 1, add = FALSE,...)
}
\arguments{
    \item{x}{an object of type \code{\link{bivrec}} or \code{\link{unirec}}.}
    \item{which}{an integer that determines which plot to draw. If \code{0},
        plots are drawn for both processes, if \code{1} or \code{2}, a plot is
        drawn for process 1 or 2 respectively.}
    \item{main}{for \code{bivrec}, a character vector of length 2, with two plot
        titles, for \code{unirec}, a string giving the main title.}
    \item{xscale}{amount by which to scale the x axis (for example to convert time
        from days to years). The hazard is automatically scaled accordingly.}
    \item{hazscale}{additional scaling for the hazard.}
    \item{add}{logical, whether to add the lines to the current plot or create
        a new plot.}
    \item{legend}{logical, whether to draw a legend.}
    \item{...}{additional parameters passed on to \code{\link{plot}}}
}
\seealso{\code{\link{bivrec}}, \code{\link{unirec}}}
\examples{
    data(m10Ji5.fit) # Example fitted model (see example(bivrec))
    
    # Default plot of the survivor function of the first process
    plot(m10Ji5.fit, which=1, main="Example Plot")
    
    # Add a line showing the effect of the covariate
    plot(m10Ji5.fit, which=1, 
        hazscale=exp(m10Ji5.fit$regression$coefficients1["Z1"]), 
        add=TRUE, lty=2)
        
    legend("topright", c("Z1 = 0", "Z1 = 1"), lty=c(1,2)) 
}
\keyword{survival}