\name{poscar.calculatewoodnotation}
\alias{poscar.calculatewoodnotation}
\title{Gives you the wood notation}
\usage{
  poscar.calculatewoodnotation(poscar1, poscar2,
    transmatrix1, transmatrix2, round = 1, latex = F)
}
\arguments{
  \item{poscar1}{object of type poscar of the substrate}

  \item{poscar2}{object of type poscar of the overlayer}

  \item{transmatrix1}{2x2 matrix which will transform the
  basis of the substrate into the basis of the supercell}

  \item{transmatrix2}{2x2 matrix which will transform the
  basis of the overlayer into the basis of the supercell}

  \item{latex}{if \code{TRUE} will give the output in latex
  format}
}
\description{
  \code{poscar.calculatewoodnotation} calculates the
  woodnotation for a slab system
}

