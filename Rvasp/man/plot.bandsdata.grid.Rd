\name{plot.bandsdata.grid}
\alias{plot.bandsdata.grid}
\title{Plots the 2d grid used in bandsdata}
\usage{
  plot.bandsdata.grid(bandsdata, sym = NA, foldback = T,
    maxdistancetobz = 1e-04, xlim = NA, ylim = NA, ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{sym}{See \code{\link{dataframe.applysymoperations}}
  for usage.}

  \item{foldback}{folding of points back in the BZ can be
  disabled}

  \item{maxdistancetobz}{(optional) allowed kpoint distance
  to the first brillouinzone (only if foldback is
  \code{True})}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata.grid} plots the 2d grid used in
  bandsdata. Main purpose is to test sym operations, to
  come to a statisfying grid.
}

