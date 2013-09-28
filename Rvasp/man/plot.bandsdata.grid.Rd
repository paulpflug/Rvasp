\name{plot.bandsdata.grid}
\alias{plot.bandsdata.grid}
\title{Plots the 2d grid used in bandsdata}
\usage{
  plot.bandsdata.grid(bandsdata, sym = NA, foldback = T,
    xlim = NA, ylim = NA, ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{sym}{See \code{\link{dataframe.applysymoperations}}
  for usage.}

  \item{foldback}{folding of points back in the BZ can be
  disabled}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata.grid} plots the 2d grid used in
  bandsdata. Main purpose is to test sym operations, to
  come to a statisfying grid.
}

