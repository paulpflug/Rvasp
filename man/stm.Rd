\name{stm}
\alias{stm}
\title{Calculates a constant-current stm by a given chgcar.
preferred orientation: z-direction
searches bottom-up (Will give the bottom side of a slab)}
\usage{
  stm(chgcar, emax, cpus = 4)
}
\arguments{
  \item{chgcar}{object of class chgcar}

  \item{emax}{cut-off electron density}

  \item{cpus}{uses snowfall package to parallelize
  calculation}
}
\description{
  \code{stm} Calculates a constant-current stm
}

