\name{plot.poscar.addlayers}
\alias{plot.poscar.addlayers}
\title{Adds layers of atoms to existing plot}
\usage{
  plot.poscar.addlayers(poscar, layer, layers,
    direction = 3, color = rainbow(length(layer)),
    size = rep(1, length(layer)),
    lwd = rep(1, length(layer)), overwritestyle = list(),
    darkenlayer = NA, darkenalpha = 0.33, ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{direction}{of projection}

  \item{layer}{vector of indices of layers to plot}

  \item{layers}{total layer count in poscar}

  \item{color}{vector of colors, will be recycled along
  \code{layer}}

  \item{size}{vector of sizes, will be recycled along
  \code{layer}}

  \item{lwd}{vector of linewides, will be recycled along
  \code{layer}}

  \item{overwritestyle}{list of elements which style will
  be overwritten (see description)}

  \item{darkenlayer}{indices of layers which will be
  additionaly darkend}

  \item{darkenalpha}{transparancy which will be used to
  darken layers in \code{darkenlayer}}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.poscar.addlayers} adds layers of atoms to
  existing plot. All styleinformation will be the same in
  one layer. However, to highligt specific atom types you
  can use \code{overwritestyle} parameter. For example
  \code{list(H=list(col="black"),He=list(size=5))} will
  color all H atoms in black and all He atoms in size 5,
  regardless of layer.
}

