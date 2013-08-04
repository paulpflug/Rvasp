\name{chgcar.sumoverlayer}
\alias{chgcar.sumoverlayer}
\title{sums the charge over layers}
\usage{
  chgcar.sumoverlayer(chgcar, layer, layers)
}
\arguments{
  \item{chgcar}{object of chgcar class}

  \item{layer}{layers which should be summed}

  \item{layers}{total layer count}
}
\description{
  \code{chgcar.sumoverlayer} sums the charge over layers.
  Only for coherent layers. Only in z-direction. As a
  cutoff the middle between highest / lowest selected and
  next unselected layer is used. If highest / lowest
  selected layer is top / bottom layer of the slab, the
  cut-off is set to z=0 / z=max(z)
}

