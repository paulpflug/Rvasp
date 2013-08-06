\name{plot.poscar.addnumbers}
\alias{plot.poscar.addnumbers}
\title{Adds numbers on atoms to existing plot}
\usage{
  plot.poscar.addnumbers(poscar, layers = 1, layer = 1,
    direction = 3, absolutenumber = T, ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{direction}{of projection}

  \item{layer}{vector of indices of layers to plot}

  \item{layers}{total layer count in poscar}

  \item{absolutenumber}{determines if index of atom in
  poscar(absolute) or index of atom in layer is printed}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.poscar.addnumbers} adds numbers on atoms,
  based on layers, to existing plot.
}

