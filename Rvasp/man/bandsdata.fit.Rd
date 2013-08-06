\name{bandsdata.fit}
\alias{bandsdata.fit}
\title{Fit a function to a band}
\usage{
  bandsdata.fit(bandsdata, bandnr, kpoints,
    fitname = c("dirac", "quadratic"), startingparameters,
    constants = NULL)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bandnr}{index of band to use}

  \item{kpoints}{indices of kpoints to use}

  \item{fitname}{used fitfunction}

  \item{startingparamters}{for use in fitfunction}

  \item{constants}{for use in fitfunction}
}
\description{
  \code{bandsdata.fit} fit a function to a band. will
  return a object of class bandsfit.
}
\seealso{
  \code{\link{bandsdata.fit.dirac.function}} and
  \code{\link{bandsdata.fit.quadratic.function}}
}

