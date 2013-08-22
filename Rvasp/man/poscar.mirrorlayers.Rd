\name{poscar.mirrorlayers}
\alias{poscar.mirrorlayers}
\title{mirrores a layer}
\usage{
  poscar.mirrorlayers(poscar, layers, baselayer,
    mirrorlayers)
}
\arguments{
  \item{poscar}{object of class poscar}

  \item{layers}{total layer count}

  \item{baselayer}{index of layer in which the
  mirror-diagonal will be based}

  \item{mirrorlayers}{index of layers which will be
  mirrored by the diagonal}
}
\description{
  \code{poscar.mirrorlayers} gives you a poscar with all
  atoms in mirrorlayer mirrored by the diagonal going
  trough a atom of the baselayer
}

