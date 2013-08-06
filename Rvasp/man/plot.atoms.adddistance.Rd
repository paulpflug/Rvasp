\name{plot.atoms.adddistance}
\alias{plot.atoms.adddistance}
\title{Adds distance between atoms to existing plot}
\usage{
  plot.atoms.adddistance(atoms, basis = NULL,
    direction = 3, length = 0.1, lwd = 1, col = "black",
    selectedalpha = 150, selectedsize = 1.5, ...)
}
\arguments{
  \item{atoms}{dataframe of atoms, will only use first two
  atoms}

  \item{direction}{of projection}

  \item{basis}{basis if atoms are in direct coordinates}

  \item{length}{arrow tip length}

  \item{lwd}{arrow line width}

  \item{col}{arrow col}

  \item{selectedalpha}{transparency level of highlighted
  point over selected atom}

  \item{selectedsize}{if larger than 0 will draw white
  transparent point over selected atom for highlighting
  purpose}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.atoms.adddistance} adds distance between atoms
  to existing plot.
}

