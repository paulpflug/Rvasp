\name{poscar.getshapedposcar}
\alias{poscar.getshapedposcar}
\title{Gives you a poscar in a given area}
\usage{
  poscar.getshapedposcar(poscar, x, y,
    shape = c("rectangular"))
}
\arguments{
  \item{poscar}{object of type poscar}

  \item{x}{max x (range will be c(0,x))}

  \item{y}{max y (range will be c(0,y))}

  \item{shape}{of the area in which atoms are selected}
}
\description{
  \code{poscar.getshapedposcar} gives you a poscar with
  atoms in area of given size and shape. Creates atoms
  beyond the border of the unitcell.
}

