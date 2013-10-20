\name{plot.bandsdata.addbands}
\alias{plot.bandsdata.addbands}
\title{Adds bands to existing plot}
\usage{
  plot.bandsdata.addbands(bandsdata,
    bands = 1:length(bandsdata$bands), energyoffset = NULL,
    ...)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{bands}{limits plotting to specified bands}

  \item{energyoffset}{will be added to energy of all bands}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.bandsdata.addbands} adds bands to existing
  plot
}

