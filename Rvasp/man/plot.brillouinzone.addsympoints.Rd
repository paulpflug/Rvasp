\name{plot.brillouinzone.addsympoints}
\alias{plot.brillouinzone.addsympoints}
\title{Adds high symmetry points to existing plot}
\usage{
  plot.brillouinzone.addsympoints(brillouinzone,
    directcoordinates = list(c(0, 0)),
    labels = 1:length(directcoordinates), vectors = 1:2,
    textpos = 2, col = "blue", pch = 16, typ = "p",
    xoffset = 0, yoffset = 0, ...)
}
\arguments{
  \item{brillouinzone}{object of class brillouinzone}

  \item{directcoordinates}{positions of highsymmetry points
  in direct coordinates}

  \item{labels}{labels of highsymmetry points in order
  given by \code{directcoordinates}}

  \item{vectors}{indices of two vectors from brilloinzone,
  which will be used as basis for calculation}

  \item{textpos}{position of labels, see \code{\link{text}}
  for futher information}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.brillouinzone.addsympoints} adds high symmetry
  points to existing plot based on poscar.
}

