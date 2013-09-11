\name{poscar.getbrillouinzone}
\alias{poscar.getbrillouinzone}
\title{Gets  brillouinzone vector}
\usage{
  poscar.getbrillouinzone(poscar, rotate = 0, extend = 1,
    strain = 0)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{rotate}{rotates brillouinzone (in degrees)}

  \item{extend}{creates supercell of brillouinzone}

  \item{strain}{applied to brillouinzone}
}
\description{
  \code{poscar.getbrillouinzone.hexagonal} gets vector
  containing the vertices of the brillouinzone.
}

