\name{poscar.createsupercell}
\alias{poscar.createsupercell}
\title{Creates supercell}
\usage{
  poscar.createsupercell(poscar, A = diag(3), center = T,
    center.directions = 1:3, center.position = rep(0.5, 3))
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{A}{a 3x3 matrix which will connect new with old
  basis \code{B_new = A x B_old}}

  \item{center}{determines if new cell should be centered}

  \item{center.directions}{subset of 1,2,3, determines
  which directions should be centered (see
  \code{atoms.centeratoms})}

  \item{center.position}{relativ to which the atoms should
  be aranged (3d vector) (see \code{atoms.centeratoms})}
}
\description{
  \code{poscar.createsupercell} creates supercell.
}

