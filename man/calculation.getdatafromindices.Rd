\name{calculation.getdatafromindices}
\alias{calculation.getdatafromindices}
\title{gives the CHGCAR data by a list of given dataindices}
\usage{
  calculation.getdatafromindices(calculation, indices,
    efermi = 0, smearing = NULL)
}
\arguments{
  \item{calculation}{object of class calculation}

  \item{indices}{dataindices based on positions}

  \item{efermi}{(optional) alignment to the Fermi level}

  \item{smearing}{(optional) adds gaussian smearing}
}
\description{
  \code{calculation.getdatafromindices} gives the CHGCAR
  data by a given dataindex
}

