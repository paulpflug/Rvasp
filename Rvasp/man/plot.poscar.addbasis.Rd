\name{plot.poscar.addbasis}
\alias{plot.poscar.addbasis}
\title{Adds basis to existing plot}
\usage{
  plot.poscar.addbasis(poscar, xoffin = 0.1, yoffin = 0.1,
    direction = 3, basisnames = c("a", "b", "c"),
    arrowlength = 0, arrowsize = 0.1, fullcell = F, ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{direction}{of projection}

  \item{xoffin}{offset in x direction in inch}

  \item{yoffin}{offset in y direction in inch}

  \item{basisnames}{names of basis vectors}

  \item{fullcell}{if activated, draws full cell, rather
  than only the volume where atoms are}

  \item{arrowlength}{length of basis arrows in user
  coordinates}

  \item{arrowsize}{size of arrow tip}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.poscar.addbasis} adds basis to existing plot.
}

