\name{plot.poscar}
\alias{plot.poscar}
\title{Custom plot for object of class poscar}
\usage{
  plot.poscar(poscar, direction = 3, xlab = "x",
    ylab = "y", basis = F, unitcell = F, fullcell = F, ...)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{direction}{of projection}

  \item{basis}{if activated, draws basis}

  \item{unitcell}{if activated, draws basis}

  \item{fullcell}{if activated, draws full cell, rather
  than only the volume where atoms are}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.poscar} custom plot for object of class
  poscar. Creates only the plotbox for 2d projection plot.
}

