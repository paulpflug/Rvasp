\name{bandsdata.calcsympointpath}
\alias{bandsdata.calcsympointpath}
\title{Adds a path of high symmetry points to bandsdata}
\usage{
  bandsdata.calcsympointpath(bandsdata, sympointpath)
}
\arguments{
  \item{bandsdata}{object of class bandsdata}

  \item{sympointpath}{list of vectors containing indices of
  high symmetry points}
}
\description{
  \code{bandsdata.calcsympointpath} adds a custom path of
  high symmetry points to bandsdata.  ex.
  sympointpath=list(c(1,2),c(3,4)) goes from 1 to 2 and 3
  to 4 returns class of type bandsdata
}

