\name{load.calculations}
\alias{load.calculations}
\title{Loads a calculations object in RData format}
\usage{
  load.calculations(name, folder = "", update = F)
}
\arguments{
  \item{name}{calculations object to load}

  \item{folder}{optional subfolder of working directory}

  \item{update}{if object is already loaded determines, if
  object will be loaded again}
}
\description{
  \code{load.calculations} loads a calculations object in
  RData format. Will search current working directory for
  \code{name}.RData
}

