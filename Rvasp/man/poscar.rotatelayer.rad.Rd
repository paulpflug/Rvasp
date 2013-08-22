\name{poscar.rotatelayer.rad}
\alias{poscar.rotatelayer.rad}
\title{Rotates layers (2d)}
\usage{
  poscar.rotatelayer.rad(poscar, layer, layers, angle,
    offset = NULL)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layer}{indices of layers which will be rotated}

  \item{layers}{total layer count}

  \item{angle}{in radians}

  \item{offset}{point around which is rotated. If not
  provied will use atom nearest to zero (x-y projected)}
}
\description{
  \code{poscar.rotatelayer.deg} rotates layers around point
  (2d). Will rotate in x-y plane
}

