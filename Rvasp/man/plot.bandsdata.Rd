\name{plot.bandsdata}
\alias{plot.bandsdata}
\title{Plots bandsdata}
\usage{
  plot.bandsdata(bandsdata,
    bands = 1:length(bandsdata$bands), sympointpath = NULL,
    col.bands = "black", type = "l", fermi = F,
    fermi.col = "blue", sym.labels = NULL,
    sym.col.labels = "black", sym.lty = 3,
    sym.col.lines = "red", xlim = NULL, xaxs = "i",
    yaxs = "i", xlab = "Wave vector", ylab = "Energy (eV)",
    energyoffset = 0, ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bands}{limits plotting to specified bands}

  \item{col.bands}{color of bands}

  \item{sympointpath}{calls
  \code{\link{bandsdata.calcsympointpath}}}

  \item{fermi}{adds line at Fermi level}

  \item{fermi.col}{color of line at Fermi level}

  \item{sym.labels}{adds labels at high symmetry points}

  \item{sym.col.labels}{color of high symmetry point labes}

  \item{sym.lty}{line typ of high symmetry point lines}

  \item{sym.col.lines}{color of high symmetry point lines}

  \item{energyoffset}{will be added to energy of all bands}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata} plots bandsdata. returns class of
  type bandsdata
}

