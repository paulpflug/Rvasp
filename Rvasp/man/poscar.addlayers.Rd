\name{poscar.addlayers}
\alias{poscar.addlayers}
\title{Insert atoms in object of class poscar}
\usage{
  poscar.addlayers(poscar, layer, layers, atomcart,
    type = c("above", "below"), distabove = NULL,
    distbelow = NULL)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layer}{index of layer which will be used for
  positioning}

  \item{layers}{total layer count}

  \item{atomcart}{atoms in cartesian coordinates}

  \item{type}{either above or below}

  \item{distabove}{(optional) custom distance to layer
  above}

  \item{distbelow}{(optional) custom distance to layer
  below}
}
\description{
  \code{poscar.addlayers} adds atoms in object of class
  poscar.
}

