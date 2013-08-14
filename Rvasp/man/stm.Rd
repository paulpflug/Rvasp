\name{stm}
\alias{stm}
\title{Calculates a constant-current stm}
\usage{
  stm(chgcar, emax, direction = 3, cpus = 4,
    interpolation = c("linear"))
}
\arguments{
  \item{chgcar}{object of class chgcar}

  \item{emax}{cut-off electron density}

  \item{direction}{direction of stm creation (negativ will
  go from small values to big, positiv vice versa)}

  \item{cpus}{uses snowfall package to parallelize
  calculation}

  \item{interpolation}{only linear implemented}
}
\description{
  \code{stm} Calculates a constant-current stm by a given
  chgcar. preferred orientation: z-direction
}

