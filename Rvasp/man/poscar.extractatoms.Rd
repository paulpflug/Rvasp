\name{poscar.extractatoms}
\alias{poscar.extractatoms}
\title{Extracts atoms of object of class poscar}
\usage{
  poscar.extractatoms(poscar, atomindices,
    vacuum = c(0, 0, 0), center = F)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{atomindices}{indices of atoms which will be
  selected}

  \item{vacuum}{sets vacuum, see \code{poscar.setvacuum}
  for futher informations}

  \item{center}{centers new poscar,see
  \code{poscar.setvacuum} for futher informations}
}
\description{
  \code{poscar.extractatoms} extracts atoms, based on
  atomindices, of object of class poscar. Will give object
  of class poscar
}

