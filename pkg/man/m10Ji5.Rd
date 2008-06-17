 \name{m10Ji5}
  \docType{data}
  \alias{m10Ji5}
  \alias{m10Ji5.fit}
  \title{Simulated data set with 10 clusters of 5 patients}
  \description{
    A simulated data set generated using the following settings:
   \itemize{
        \item 10 clusters
        \item 5 subjects per cluster
        \item 1 stratum
        \item lognormal cluster frailties with mean 1, variance 0.25
        \item lognormal subject frailties with variance .25 and covariance .125
        \item one time-fixed covariate, generated as Normal(0,.5)
        \item true regression coefficients are 1 for both processes
        \item both baseline hazards are Weibull with lambda=10, gamma=1.8
        \item censoring times are Weibull with lambda=1, gamma=1.8
   }
    used for package examples. The object \code{m10Ji5.fit} contains the
    fitted bivariate model for this data with default settings.
  }
  \usage{m10Ji5}
  \format{\code{m10Ji5} is a data frame of 9 columns and 298 rows, 
    \code{m10Ji5.fit} is an object of class \code{\link{bivrec}}.}
  \seealso{\code{\link{bivrec}}}
  \keyword{datasets}