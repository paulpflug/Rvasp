\name{poscar.getatomlayerindices}
\alias{poscar.getatomlayerindices}
\title{Gives layerindices for all atoms}
\usage{
  poscar.getatomlayerindices(poscar, layers)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layers}{total layer count}
}
\description{
  \code{poscar.getatomlayerindices} will give a vector
  containig indices for all atoms, correspondig to the
  related layer. Will only work in z-direction.
}

