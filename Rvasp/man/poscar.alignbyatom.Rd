\name{poscar.alignbyatom}
\alias{poscar.alignbyatom}
\title{aligns the atoms by a certain position}
\usage{
  poscar.alignbyatom(poscar, atomindex, alignto = c(0, 0),
    direction = 3)
}
\arguments{
  \item{poscar}{object of type poscar}

  \item{atomindex}{index of atom which will be used}

  \item{alignto}{will move all atoms, so that given atom
  will be positioned at \code{alignto}}

  \item{direction}{will determine in which plane the
  translation will be performed}
}
\description{
  \code{poscar.alignbyatom} aligns the atoms by a certain
  position. The position will be searched by the index of a
  atom. Will make a 2d translation
}

