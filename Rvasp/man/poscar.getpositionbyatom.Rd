\name{poscar.getpositionbyatom}
\alias{poscar.getpositionbyatom}
\title{Gives positions object based on atoms}
\usage{
  poscar.getpositionbyatom(poscar, atomselector,
    zdist = c(1), layername = -1)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{atomselector}{indices of atoms which will be used
  for position calculation}

  \item{zdist}{distance of position to atom in z-direction}

  \item{layername}{name for automated legend creation}
}
\description{
  \code{poscar.getpositionbyatom} calculates positions
  based on atom positions. Will give object of class
  positions, which is used for LDOS calculations.
}

