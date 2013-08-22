\name{poscar.getpositionbylayer}
\alias{poscar.getpositionbylayer}
\title{Gives positions object based on atoms}
\usage{
  poscar.getpositionbylayer(poscar, layers, zdist = c(1),
    layer = layers)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layer}{indices of layers which will be used to
  search atom positions}

  \item{layers}{total layer count}

  \item{zdist}{distance of position to atom in z-direction}
}
\description{
  \code{poscar.getpositionbylayer} calculates positions
  based on atom positions. Will give object of class
  positions, which is used for LDOS calculations.
}

