\name{dataframe.applysymoperations}
\alias{dataframe.applysymoperations}
\title{Applies chain of symmetric operation on a dataframe}
\usage{
  dataframe.applysymoperations(dataframe, symoperations,
    center = NULL)
}
\arguments{
  \item{dataframe}{dataframe with at least two columns}

  \item{symoperations}{list of symmetric operations.
  \code{list(c("reflection",-30),c("rotation",120,60))} for
  example will first reflect across a line with the slope
  of -30°, afterwards the new dataset (old+reflected) will
  be rotated by 120° and 60° the resulting dataset will
  consist of 6 combined datasets. The old, the reflected
  and the old+reflected, rotated by 120° and by 60°.}

  \item{center}{(optional) 2d point, which will be used as
  center for all symoperations}
}
\description{
  \code{dataframe.applysymoperations} applied a chain of 2d
  symmetric operations on the first two columns of the
  dataframe
}

