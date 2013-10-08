\name{makeTransparent}
\alias{makeTransparent}
\title{makes a vector of colors transparent
snippet from stackoverflow.com}
\usage{
  makeTransparent(someColor, alpha = 100, real = T,
    real.backgroundcolor = "white")
}
\arguments{
  \item{someColor}{colors}

  \item{alpha}{transparency}

  \item{real}{if \code{FALSE} will return a mixture with
  \code{real.backgroundcolor}, instead real transparent
  color.}

  \item{real.backgroundcolor}{color which will be used for
  mixing, if no real transparency can be used.}
}
\description{
  \code{makeTransparent} makes a vector of colors
  transparent
}

