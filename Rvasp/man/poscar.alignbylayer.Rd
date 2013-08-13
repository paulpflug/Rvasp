\name{poscar.alignbylayer}
\alias{poscar.alignbylayer}
\title{aligns the atoms by a certain position}
\usage{
  poscar.alignbylayer(poscar, layer, layers,
    atominlayer = 1, alignto = c(0, 0), direction = 3)
}
\arguments{
  \item{poscar}{object of type poscar}

  \item{layer}{the layer in which the atom will be
  searched}

  \item{layers}{total layer count}

  \item{atominlayer}{index of atom in given layer}

  \item{alignto}{will move all atoms, so that given atom
  will be positioned at \code{alignto}}

  \item{direction}{will determine in which plane the
  translation will be performed}
}
\description{
  \code{poscar.alignbylayer} aligns the atoms by a certain
  position. The position will be searched by a atom in a
  given layer. Will make a 2d translation
}

