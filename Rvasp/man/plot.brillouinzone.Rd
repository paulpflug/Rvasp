\name{plot.brillouinzone}
\alias{plot.brillouinzone}
\title{Plots brillouinzone}
\usage{
  plot.brillouinzone(poscar, type = "hexagonal",
    rotate = 0, extend = 1, strain = 0, xlab = NA,
    ylab = NA, xaxt = "n", yaxt = "n", col = "blue",
    typ = "l", ...)
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
  \code{plot.brillouinzone} plots brillouinzone based on
  poscar.
}

