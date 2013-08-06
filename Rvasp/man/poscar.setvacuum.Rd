\name{poscar.setvacuum}
\alias{poscar.setvacuum}
\title{Sets the vacuum}
\usage{
  poscar.setvacuum(poscar, vacuum = c(0, 0, 0), center = T)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{vacuum}{vector of new vacuum in all three
  directions. \code{c(0,0,10)} will change the vacuum in
  third direction to 10 (and leave the others).}

  \item{center}{\code{TRUE} will center the cell after
  change}
}
\description{
  \code{poscar.setvacuum} changes the vacuum. Vacuum is
  defined by the distance between atoms across unitcell
  borders.
}

