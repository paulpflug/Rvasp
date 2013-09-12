\name{plot.bandsdata.persp}
\alias{plot.bandsdata.persp}
\title{Plots the perspective bandsdata}
\usage{
  plot.bandsdata.persp(bandsdata, bands,
    maxdistancetobz = 1e-04, sym = NA, n = 201,
    projected = F, projected.atoms = 1:bandsdata$natoms,
    projected.energyintervall = NULL,
    projected.bands = 1:bandsdata$nbands,
    projected.orbitals = list(1, 2, 3, 4),
    projected.colorpalette = colorRampPalette(c("red", "blue", "green")),
    xlab = expression(k[x]), ylab = expression(k[y]), ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bands}{indices of bands to plot}

  \item{maxdistancetobz}{(optional) allowed kpoint distance
  to the first brillouinzone}

  \item{n}{resolution of datalayers}

  \item{sym}{See \code{\link{dataframe.applysymoperations}}
  for usage and \code{\link{plot.bandsdata.grid}} for
  testing.}

  \item{projected}{activate for additional datalayer with
  projected states.}

  \item{projected.atoms}{indices of atoms over which will
  be summed}

  \item{projected.energyintervall}{in which bands will be
  included (overwrites \code{projected.bands})}

  \item{projected.bands}{used for normating of color.
  searches for maximum projected value in these bands}

  \item{projected.orbitals}{list of orbitals to plot. To
  sum over orbitals 2 and 3 use \code{list(1,c(2,3),4)}}

  \item{projected.colorpalette}{color palette for orbitals}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata.persp} plots the perspective
  bandsdata. Projected states can be added as additional
  datalayer. First make sure, you have a statisfying grid.
  See \code{\link{plot.bandsdata.grid}}.
}

