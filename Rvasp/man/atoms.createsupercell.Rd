\name{atoms.createsupercell}
\alias{atoms.createsupercell}
\title{Creates simple supercell}
\usage{
  atoms.createsupercell(atomsdirect, super = c(1, 1, 1),
    center = T, center.directions = 1:3,
    center.position = rep(0, 3))
}
\arguments{
  \item{atomsdirect}{atoms in direct coordinates}

  \item{super}{3d whole-number vector which defines, size
  of supercell in 3 directions}

  \item{center}{determines if new cell should be centered}

  \item{center.directions}{subset of 1,2,3 determines which
  directions should be centered (see
  \code{atoms.centeratoms})}

  \item{center.position}{relativ to which the atoms should
  be aranged (3d vector) (see \code{atoms.centeratoms})}
}
\description{
  \code{atoms.createsupercell} creates simple supercell.
  Only whole-number supercells.
}

