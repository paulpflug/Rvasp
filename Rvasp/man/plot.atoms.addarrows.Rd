\name{plot.atoms.addarrows}
\alias{plot.atoms.addarrows}
\title{Adds arrows from specific atoms to other atoms to existing plot}
\usage{
  plot.atoms.addarrows(atomsold, atomsnew, basisold = NULL,
    basisnew = NULL, direction = 3, length = 0.1, ...)
}
\arguments{
  \item{atomsold}{dataframe of atoms, used for arrow start}

  \item{atomsnew}{dataframe of atoms, used for arrow end}

  \item{direction}{of projection}

  \item{basisold}{basis if atomsold are in direct
  coordinates}

  \item{basisnew}{basis if atomsold are in direct
  coordinates}

  \item{length}{size of arrow tips}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.atoms.addarrows} adds arrows from specific
  atoms to other atoms to existing plot. Arrows will be
  drawn pairwise.
}

