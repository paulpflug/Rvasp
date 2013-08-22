\name{poscar.removelayers}
\alias{poscar.removelayers}
\title{Removes atoms of object of class poscar}
\usage{
  poscar.removelayers(poscar, layer, layers,
    vacuum = c(0, 0, 0), center = T)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layer}{indices of layers which will be removed}

  \item{layers}{total layer count}

  \item{vacuum}{sets vacuum, see \code{poscar.setvacuum}
  for futher informations}

  \item{center}{centers new poscar,see
  \code{poscar.setvacuum} for futher informations}
}
\description{
  \code{poscar.removelayers} removes atoms, based on
  layers, of object of class poscar. Will give object of
  class poscar. Remaining layers will be moved. The
  distance between two layers next to the removed layer is
  set to the mean of the two distances of the removed
  layer. For top / bottom layer, the distance to zero /
  zmax is used. Works best for centered cells.
}

