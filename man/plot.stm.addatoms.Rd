\name{plot.stm.addatoms}
\alias{plot.stm.addatoms}
\title{Adds atoms to a stm image}
\usage{
  plot.stm.addatoms(stm, atomselector, col = "black",
    atomsize = 2 * super, super = 1, xlim = NULL,
    ylim = NULL, ...)
}
\arguments{
  \item{stm}{object of class stm}

  \item{atomselector}{selector to filter atoms which should
  be added}

  \item{super}{(optional) makes a supercell of filtered
  atoms (symmetric in x and y)}

  \item{...}{further plotting parameters}
}
\description{
  \code{plot.stm.addatoms} Adds atoms
}

