\name{poscar.centeratoms}
\alias{poscar.centeratoms}
\title{Centers atoms in object of class poscar}
\usage{
  poscar.centeratoms(poscar, directions = 1:3,
    position = rep(0.5, 3))
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{directions}{subset of 1,2,3 determines which
  directions should be centered}

  \item{position}{relativ to which the atoms should be
  aranged (3d vector)}
}
\description{
  \code{poscar.centeratoms} centers the atoms relativ to a
  specific position in object of class poscar.
}

