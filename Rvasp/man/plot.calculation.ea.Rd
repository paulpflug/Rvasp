\name{plot.calculation.ea}
\alias{plot.calculation.ea}
\title{Will plot a e over a curve}
\usage{
  plot.calculation.ea(calculation, energyfactor = 1,
    energyshift = 0, fit = F, type = "p", ...)
}
\arguments{
  \item{calculation}{object of type calculation}

  \item{energyfactor}{will be multiplied with energy}

  \item{energyshift}{will be substracted from energy}

  \item{fit}{if \code{TRUE} will fit with
  \code{\link{ea.fitEOS}} plot and return result}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.calculation.ea} Will plot a e over a curve.
}

