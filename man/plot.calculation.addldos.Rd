\name{plot.calculation.addldos}
\alias{plot.calculation.addldos}
\title{adds a ldos to an existing plot. either positions or dataset must be given.}
\usage{
  plot.calculation.addldos(calculation, positions = NULL,
    efermi = 0, smearing = NULL, dataset = NULL,
    col = "black", lty = 1, ...)
}
\arguments{
  \item{calculation}{object of class calculation}

  \item{positions}{object of class positions}

  \item{dataset}{list of ldos data objects}

  \item{efermi}{(optional) alignment to the Fermi level}

  \item{smearing}{(optional) adds gaussian smearing}
}
\description{
  \code{plot.calculation.addldos} adds a ldos to an
  existing plot
}

