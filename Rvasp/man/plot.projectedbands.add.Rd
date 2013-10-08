\name{plot.projectedbands.add}
\alias{plot.projectedbands.add}
\title{Adds projected orbitals to existing plot}
\usage{
  plot.projectedbands.add(projectedbands,
    bands = 1:length(projectedbands$bands),
    orbitals = list(1, 2, 3, 4), col = NULL,
    col.palette = colorRampPalette(c("red", "blue", "green")),
    pch = 15:(14 + length(orbitals)), cex = 0.8,
    cex.cutoff = 0.2, legend = "topright", legendcex = 0.8,
    energyoffset = NULL, usetransparent = T, ...)
}
\arguments{
  \item{projectedbands}{object of class projectedbands}

  \item{bands}{limits plotting to specified bands}

  \item{orbitals}{list of orbitals to plot. To sum over
  orbitals 2 and 3 use \code{list(1,c(2,3),4)}}

  \item{col}{colors for orbitals, will be recycled by the
  length of orbitals, if not provided, will use col.platte
  for generation.}

  \item{col.palette}{color palette for orbitals}

  \item{legend}{position of legend, \code{NULL} will
  supress plotting}

  \item{legendcex}{size of legend}

  \item{energyoffset}{will be added to energy of all bands}

  \item{usetransparent}{determines usage of transparency}

  \item{cex.cutoff}{will only be used if usetransparent
  equals \code{FALSE}. Prevents plotting of smaller than
  cutoff points.}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.projectedbands.add} adds projected orbitals to
  existing plot.
}

