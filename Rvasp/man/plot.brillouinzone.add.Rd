\name{plot.brillouinzone.add}
\alias{plot.brillouinzone.add}
\title{Adds brilloinzone to existing plot}
\usage{
  plot.brillouinzone.add(poscar, type = "hexagonal",
    rotate = 0, extend = 1, strain = 0, col = "blue", ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{type}{of brillouinzone (currently only hexagonal)}

  \item{rotate}{rotates brillouinzone (in degrees)}

  \item{extend}{creates supercell of brillouinzone}

  \item{strain}{applied to brillouinzone}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.brillouinzone.add} adds hexagonal brilloinzone
  to existing plot based on poscar.
}

