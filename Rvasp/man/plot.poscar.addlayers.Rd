\name{plot.poscar.addlayers}
\alias{plot.poscar.addlayers}
\title{Adds layers of atoms to existing plot}
\usage{
  plot.poscar.addlayers(poscar, layer, layers,
    direction = 3, color = rainbow(length(layer)),
    size = rep(1, length(layer)),
    lwd = rep(1, length(layer)), ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{direction}{of projection}

  \item{layer}{vector of indices of layers to plot}

  \item{layers}{total layer count in poscar}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.poscar.addlayers} adds layers of atoms to
  existing plot.
}

