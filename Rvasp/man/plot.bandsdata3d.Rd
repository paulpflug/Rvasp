\name{plot.bandsdata3d}
\alias{plot.bandsdata3d}
\title{Plots the 3d bandsdata}
\usage{
  plot.bandsdata3d(bandsdata, band, projecttobz = T,
    sym = NA, n = 201, projected = F,
    projected.atoms = 1:bandsdata$natoms,
    projected.energyintervall = NULL,
    projected.bands = 1:bandsdata$nbands,
    projected.orbitals = list(1, 2, 3, 4),
    projected.colorpalette = colorRampPalette(c("red", "blue", "green")),
    xlab = expression(k[x]), ylab = expression(k[y]), ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{band}{index of band to plot}

  \item{projecttobz}{(optional) projects kpoints to the
  first brillouinzone}

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

  \item{energyintervall}{}
}
\description{
  \code{plot.bandsdata3d} plots the 3d bandsdata. Projected
  states can be added as additional datalayer. First make
  sure, you have a statisfying grid. See
  \code{\link{plot.bandsdata.grid}}.
}

