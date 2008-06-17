 \name{vigndata}
  \docType{data}
  \alias{vigndata}
  \alias{vigndata.fit}
  \alias{vigndata.quickfit}
  \alias{vigndata.quickfit2}
  \title{Simulated data set used in the blupsurv package vignette}
  \description{
    A simulated data set for illustrative use in the package vignette. To make
    the discussion less abstract, covariates were given real-world names. The
    data represent a 10-year study during which patients in the 50 US states
    were monitored for events of severe pain or fever.
    The data were generated as follows:
   \itemize{
        \item 50 clusters, labeled with US state names
        \item between 5 and 25 subjects per cluster
        \item 1 stratum
        \item frailties for the pain process are lognormal with 
              cluster and subject- level variance 0.25
        \item frailties for the fever process are lognormal with
              cluster and subject- level variance 0.5
        \item subject-level frailty covariance is 0.15
        \item age is a rounded Normal variable with mean 50 and variance 10
        \item sex is a Bernoulli variable with mean 0.5
        \item true regression coefficients for the pain process are
              0 for age and 1 for sex
        \item true regression coefficient for the fever process are
              0.025 for age and 2 for sex
        \item both baseline hazards are Weibull with lambda=0., gamma=1.8
        \item censoring times are fixed at 10
   }
    used for package examples. The object \code{vigndata.fit} contains the
    fitted bivariate model for this data with full discretization.
  }
  \usage{vigndata}
  \format{\code{vigndata} is a data frame of 8 columns and 1982 rows, 
    \code{vigndata.fit} and \code{vigndata.quickfit}
    are objects of class \code{\link{bivrec}}.}
  \seealso{\code{\link{bivrec}}}
  \keyword{datasets}