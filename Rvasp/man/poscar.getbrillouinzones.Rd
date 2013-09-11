\name{poscar.getbrillouinzones}
\alias{poscar.getbrillouinzones}
\title{Gets brillouinzone vector}
\usage{
  poscar.getbrillouinzones(poscar, rotate = 0, extend = 1,
    strain = 0)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{rotate}{rotates brillouinzone (in degrees)}

  \item{extend}{creates supercell of brillouinzone}

  \item{strain}{applied to brillouinzone}
}
\description{
  \code{poscar.getbrillouinzones} gets vector containing
  the vertices of the brillouinzones.
}

