\name{atoms.centeratoms}
\alias{atoms.centeratoms}
\title{Centers atoms}
\usage{
  atoms.centeratoms(atomsdirect, directions = 1:3,
    position = rep(0, 3))
}
\arguments{
  \item{atomsdirect}{atoms in direct coordinates}

  \item{directions}{subset of 1,2,3 determines which
  directions should be centered}

  \item{position}{relativ to which the atoms should be
  aranged (3d vector)}
}
\description{
  \code{atoms.centeratoms} centers the atoms relativ to a
  specific position.
}

