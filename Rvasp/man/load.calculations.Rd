\name{load.calculations}
\alias{load.calculations}
\title{Loads a calculations object in RData format}
\usage{
  load.calculations(file, objectname = NULL, update = F)
}
\arguments{
  \item{file}{with calculations object to load}

  \item{objectname}{will check if already loaded}

  \item{update}{if object is already loaded determines, if
  object will be loaded again}
}
\description{
  \code{load.calculations} loads a calculations object in
  RData format. Will search current working directory for
  \code{name}.RData
}

