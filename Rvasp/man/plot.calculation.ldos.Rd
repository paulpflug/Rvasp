\name{plot.calculation.ldos}
\alias{plot.calculation.ldos}
\title{Creates a ldos image.}
\usage{
  plot.calculation.ldos(calculation, positions,
    col = rainbow(length(positions)), typ = "l",
    efermi = 0, smearing = NULL,
    xlab = expression("Energy" ~ epsilon ~ "(eV)"),
    ylab = expression(rho * "(x," * epsilon * ")" ~ ("Å"^-3)),
    legend = F, lty = 1, ...)
}
\arguments{
  \item{calculation}{object of class calculation}

  \item{positions}{object of class positions}

  \item{smearing}{(optional) adds gaussian smearing}

  \item{efermi}{(optional) alignment to the Fermi level}

  \item{legend}{(optional) should an automatic legend
  should be created?}
}
\description{
  \code{plot.calculation.ldos} Creates a ldos image Based
  on base package.
}

