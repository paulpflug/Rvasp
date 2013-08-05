\name{plot.dosdata.add}
\alias{plot.dosdata.add}
\title{Will add dosdata to existing plot.}
\usage{
  plot.dosdata.add(dosdata, factor = 1, smearing = 0,
    orbitals = NULL, atomindices = NULL, border = c(NA),
    col = c("grey"), type = c("line", "polygon"), ...)
}
\arguments{
  \item{dosdata}{object of type dosdata}

  \item{factor}{is multiplyed with density}

  \item{smearing}{if greater zero will add a gaussian
  smearing with \code{smearing} as standard deviation}

  \item{orbitals}{vector of indices of orbitals which are
  seperatly plotted. Use \code{all} to sum all orbitals}

  \item{atomindices}{vector of indices of atoms which are
  summed}

  \item{type}{\code{line} or \code{polygon} mode}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.dosdata.add} will add dosdata to existing
  plot. Needs partial dos. Offers two modes. Will return
  object of type dosdata where plotting parameters are
  saved.
}

