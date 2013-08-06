\name{plot.bandsdata}
\alias{plot.bandsdata}
\title{Plots bandsdata}
\usage{
  plot.bandsdata(bandsdata,
    bands = 1:length(bandsdata$bands), sympointpath = NULL,
    col = "black", type = "l", fermi = F, symnames = NULL,
    symcolor = "red", symlty = 3, xlim = NULL, xaxs = "i",
    yaxs = "i", energyoffset = 0, ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bands}{limits plotting to specified bands}

  \item{sympointpath}{calls
  \code{\link{bandsdata.calcsympointpath}}}

  \item{fermi}{adds blue line at Fermi level}

  \item{symnames}{adds labels at high symmetry points}

  \item{symcolor}{color of high symmetry point lines}

  \item{symlty}{line typ of high symmetry point lines}

  \item{energyoffset}{will be added to energy of all bands}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata} plots bandsdata. returns class of
  type bandsdata
}

