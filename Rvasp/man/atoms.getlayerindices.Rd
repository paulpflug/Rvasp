\name{atoms.getlayerindices}
\alias{atoms.getlayerindices}
\title{Gives layerindices for all atoms}
\usage{
  atoms.getlayerindices(atoms, basis = NULL, layers)
}
\arguments{
  \item{atoms}{dataframe of atoms}

  \item{basis}{needed if atoms are in direct coordinates}

  \item{layers}{total layer count}
}
\description{
  \code{atoms.getlayerindices} will give a vector containig
  indices for all atoms, correspondig to the related layer.
  Will only work in z-direction.
}

