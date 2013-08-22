\name{poscar.getlayertranslation}
\alias{poscar.getlayertranslation}
\title{Gives the 2d translation of two layers}
\usage{
  poscar.getlayertranslation(poscar, fromlayer, tolayer,
    layers)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{fromlayer}{index of base layer}

  \item{tolayer}{index of layer which will be used to
  calculate translation}

  \item{layers}{total layer count}
}
\description{
  \code{poscar.getlayertranslation} gives the 2d
  translation of two layers. Only useful for equivalent
  layers. Searches for next neighbor, only works if next
  neighbor can be found easily (e.g. supercell): Will give
  2d vector, which could be used to move \code{fromlayer}
  to \code{tolayer}
}

