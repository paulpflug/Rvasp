\name{zoomplot}
\alias{zoomplot}
\title{creates a new zoomplot within a plot
uses base package}
\usage{
  zoomplot(xpos, ypos, xmean, xwidth, ymean, ywidth = NULL,
    lines = T, box = T, mar = c(0, 0, 0, 0), xlab = "",
    ylab = "", bg = "", asp = NULL, xticks = 2, yticks = 2,
    ...)
}
\arguments{
  \item{xpos}{x position of the plot}

  \item{ypos}{y position of the plot}

  \item{xmean}{new plotting region will be centered around
  xmean}

  \item{xwidth}{width of the new plotting region in x
  direction}

  \item{ymean}{new plotting region will be centered around
  ymean}

  \item{ywidth}{width of the new plotting region in y
  direction (optional if asp is given)}

  \item{lines}{will draw lines from old plot to new
  plotting region}

  \item{box}{will draw box around old plot}

  \item{asp}{(optional) is used to calculate ywidth if not
  given. if both are not given the asp of original plot is
  used.}

  \item{mar}{margin for new plot}

  \item{...}{further plotting parameters}
}
\description{
  \code{zoomplot} creates a new plot within a plot
}

