\name{chgcar.calcdensity}
\alias{chgcar.calcdensity}
\title{Calculates the charge density}
\usage{
  chgcar.calcdensity(chgcar,
    parameter = c("r", "phi", "theta"),
    nmean = pchg$n[[1]] * 10, smearing = NA)
}
\arguments{
  \item{chgcar}{object of chgcar class}

  \item{paramter}{spherical paramter which is used}

  \item{nsum}{count of adjacent density values, which will
  be averaged}

  \item{smearing}{(optional) sigma for gaussian smearing}
}
\description{
  \code{chgcar.calcdensity} calculates the charge density
  in dependence of a spherical parameter of choice. Will
  return a dataframe with the spherical parameter and the
  calculated charge density.
}

