\name{plot.periodictable.addelements}
\alias{plot.periodictable.addelements}
\title{Adds elements to an existing periodic table}
\usage{
  plot.periodictable.addelements(elements,
    text = periodictable.getelementnames(elements),
    rows = NULL, rectmar = 0.05, col = "black",
    fill = "white", rows.col = "black", rows.cex = 0.8,
    text.cex = 1, round = NULL)
}
\arguments{
  \item{elements}{vector of (element) numbers or names}

  \item{col}{color of selected elements}

  \item{text}{additional text for selected elements list of
  vectors each containing one or two entries. List is
  recycled along selected elements}

  \item{rows.col}{color of selected texts}

  \item{fill}{color of the filling of the box}

  \item{rows.cex}{size of texts}

  \item{text.cex}{size of element name}

  \item{round}{if highlight texts are numeric type, value
  used for rounding}
}
\description{
  \code{plot.periodictable.addelements} Adds elements to an
  existing periodic table
}
\seealso{
  \code{\link{plot.periodictable.addunderlay}} and
  \code{\link{periodictable.getremainingelements}}
}

