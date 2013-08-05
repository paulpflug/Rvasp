\name{plot.dosdata}
\alias{plot.dosdata}
\title{Creates a plot for dosdata}
\usage{
  plot.dosdata(dosdata, smearing = 0, flip = F, norm = F,
    fermi = F, col.fermi = "blue", xaxs = "i", yaxs = "i",
    xlab = "Energy [eV]", ylab = "Electron density", ...)
}
\arguments{
  \item{dosdata}{object of type dosdata}

  \item{smearing}{if greater zero will add a gaussian
  smearing with \code{smearing} as standard deviation}

  \item{flip}{will exchange x and y}

  \item{norm}{will divide density by atom count}

  \item{fermi}{will draw fermi level}

  \item{col.fermi}{color of fermi level}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.dosdata} creates a plot for dosdata. Will
  return object of type dosdata where plotting parameters
  are saved. For adding curves see
  \code{\link{plot.dosdata.add}}
}

