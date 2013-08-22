\name{poscar.extractlayers}
\alias{poscar.extractlayers}
\title{Extracts atoms of object of class poscar}
\usage{
  poscar.extractlayers(poscar, layer, layers,
    vacuum = c(0, 0, 0), center = T)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layer}{indices of layers which will be extracted}

  \item{layers}{total layer count}

  \item{vacuum}{sets vacuum, see \code{poscar.setvacuum}
  for futher informations}

  \item{center}{centers new poscar,see
  \code{poscar.setvacuum} for futher informations}
}
\description{
  \code{poscar.extractlayers} extracts atoms, based on
  layers, of object of class poscar. Will give object of
  class poscar. Will not change position of atoms.
}

