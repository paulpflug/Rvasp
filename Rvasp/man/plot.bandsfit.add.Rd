\name{plot.bandsfit.add}
\alias{plot.bandsfit.add}
\title{Adds a bandsfit to an existing plot}
\usage{
  plot.bandsfit.add(bandsfit, kpoints = bandsfit$knumbers,
    n = 201, energyoffset = NULL, ...)
}
\arguments{
  \item{bandsfit}{object of class bandsfit}

  \item{kpoints}{indices of kpoints, where the fit is
  plotted. All in between points will also be used.}

  \item{n}{count of points used for drawing}

  \item{energyoffset}{will be added to energy the fit}
}
\description{
  \code{plot.bandsfit.add} adds a bandsfit to an existing
  plot
}

