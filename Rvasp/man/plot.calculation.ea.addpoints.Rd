\name{plot.calculation.ea.addpoints}
\alias{plot.calculation.ea.addpoints}
\title{Adds a curve to an existing e over a curve}
\usage{
  plot.calculation.ea.addpoints(calculation, afactor = 1,
    energyfactor = 1, energyshift = 0, fit = F, ...)
}
\arguments{
  \item{calculation}{object of type calculation}

  \item{afactor}{will be multiplied with lattice constant}

  \item{energyfactor}{will be multiplied with energy}

  \item{energyshift}{will be substracted from energy}

  \item{fit}{if \code{TRUE} will fit with
  \code{\link{ea.fitEOS}} plot and return result}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.calculation.ea.addpoints} adds a curve to an
  existing e over a curve.
}

