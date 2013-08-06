\name{bandsdata.getprojecteddata}
\alias{bandsdata.getprojecteddata}
\title{Will give a object of class projectedbands}
\usage{
  bandsdata.getprojecteddata(bandsdata,
    atomindices = 1:bandsdata$natoms,
    bands = 1:bandsdata$nbands, energyintervall = NULL,
    cpus = 1)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{atomindices}{indices of atoms over which will be
  summed}

  \item{bands}{which will be included}

  \item{energyintervall}{in which bands will be included
  (overwrites \code{bands})}
}
\description{
  \code{bandsdata.getprojecteddata} will give a object of
  class projectedbands
}

