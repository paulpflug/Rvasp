\name{plot.periodictable}
\alias{plot.periodictable}
\title{Draws a periodic table}
\usage{
  plot.periodictable(highlights = NULL,
    highlighttexts = NULL, underlay = NULL, col = "grey80",
    col.highlight = "black", col.underlay = "red",
    col.Lanthanoid = "grey60", col.Actinoid = "grey90",
    rectmar = 0.05, xlab = "", ylab = "", xaxs = "i",
    yaxs = "i", text.cex = 1, rows.cex = 0.8,
    col.rows = "black", round = NULL, typ = "", ...)
}
\arguments{
  \item{highlights}{vector of (element) numbers or names}

  \item{col.highlight}{color of highlighted elements}

  \item{highlighttexts}{text for highlighted elements list
  of vectors each containing one or two entries. List is
  recycled along highlights.}

  \item{col.rows}{color of highlight texts}

  \item{underlay}{vector of (element) numbers or names}

  \item{col.underlay}{color of rectangle underlaying
  elements}

  \item{col}{color of all not highlighted elements}

  \item{col.Lanthanoid}{underlay color of rectangle
  underlaying Lanthanoids}

  \item{col.Actinoid}{underlay color of rectangle
  underlaying Actinoids}

  \item{rectmar}{margin of boxes}

  \item{text.cex}{size of element names}

  \item{rows.cex}{size of highlight texts}

  \item{round}{if highlight texts are numeric type, value
  used for rounding}

  \item{typ}{"n" for suppress plotting}
}
\description{
  \code{plot.periodictable} Draws a periodic table. There
  are two ways of changing an entry. Highlight and
  underlay. Highlight colors the box of an entry and allows
  additional text in that box. Underlay puts a color
  beneath that box.
}
\details{
  For both parameters a vector of (element) numbers or
  names can be used.

  For very individual ploting use \code{typ="n"} and see
  \code{\link{plot.periodictable.addelements}},
  \code{\link{plot.periodictable.addunderlay}} and
  \code{\link{periodictable.getremainingelements}}
}

