\name{plot.bandsdata.3d}
\alias{plot.bandsdata.3d}
\title{Plots 3d bandsdata}
\usage{
  plot.bandsdata.3d(bandsdata, bands, sym = NA, n = 201,
    colorpalette = colorRampPalette(c("red", "white", "blue")),
    asp = 0.1, breaks = 256, zone = c("bz", "basis"),
    linear = T, projected = F,
    projected.atoms = 1:bandsdata$natoms,
    projected.orbitals = list(1, 2, 3, 4),
    projected.bands = 1:bandsdata$nbands,
    projected.energyintervall = NULL,
    projected.colors = colorRampPalette(c("red", "blue", "green"))(length(projected.orbitals)),
    xlab = expression(k[x]), ylab = expression(k[y]),
    zlab = "energy (eV)", back = "filled",
    front = "filled", box = F, ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bands}{indices of bands to plot}

  \item{n}{resolution of datalayers}

  \item{sym}{See \code{\link{dataframe.applysymoperations}}
  for usage and \code{\link{plot.bandsdata.grid}} for
  testing.}

  \item{colorpalette}{colors which are used if
  \code{projected} is \code{False}}

  \item{asp}{aspect ration between energy and k-space}

  \item{breaks}{number of steps in color coding}

  \item{zone}{see
  \code{\link{reciprocalbasis.getbrillouinzone}} for
  usage.}

  \item{linear}{see \code{\link{akima::interp}} for usage.}

  \item{projected}{activate for additional datalayer with
  projected states.}

  \item{projected.atoms}{indices of atoms over which will
  be summed}

  \item{projected.bands}{(optional) used for normating of
  color. searches for maximum projected value in these
  bands. Make sure your desired bands are included.}

  \item{projected.energyintervall}{(optional) in which
  bands will be included for norming projected. Make sure
  your desired bands are in this intervall.}

  \item{projected.orbitals}{list of orbitals to plot. To
  sum over orbitals 2 and 3 use \code{list(1,c(2,3),4)}}

  \item{projected.colorpalette}{color palette for orbitals}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata.3d} plots 3d bandsdata. Projected
  states can be added as additional datalayer. First make
  sure, you have a statisfying grid. See
  \code{\link{plot.bandsdata.grid}}.
}
\seealso{
  \code{\link{rgl}}
}

