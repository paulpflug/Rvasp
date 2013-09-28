\name{plot.bandsdata.contour}
\alias{plot.bandsdata.contour}
\title{Plots the contour bandsdata}
\usage{
  plot.bandsdata.contour(bandsdata, band, sym = NA,
    n = 201, breaks = 12,
    colorpalette = colorRampPalette(c("red", "blue", "green")),
    zone = c("bz", "basis"), linear = F, projected = F,
    projected.atoms = 1:bandsdata$natoms,
    projected.energyintervall = NULL,
    projected.bands = 1:bandsdata$nbands,
    projected.orbitals = list(1, 2, 3, 4),
    projected.colors = colorRampPalette(c("red", "blue", "green"))(length(projected.orbitals)),
    xlab = expression(k[x]), ylab = expression(k[y]), ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{band}{index of band to plot}

  \item{n}{resolution of datalayers}

  \item{sym}{See \code{\link{dataframe.applysymoperations}}
  for usage and \code{\link{plot.bandsdata.grid}} for
  testing.}

  \item{breaks}{number of steps in color coding}

  \item{colorpalette}{colors which are used if
  \code{projected} is \code{False}}

  \item{zone}{see
  \code{\link{reciprocalbasis.getbrillouinzone}} for
  usage.}

  \item{linear}{see \code{\link{akima::interp}} for usage.}

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
  \code{plot.bandsdata.contour} plots the contour bandsdata
  of a single band. Projected states can be added as
  additional datalayer. First make sure, you have a
  statisfying grid. See \code{\link{plot.bandsdata.grid}}.
}

