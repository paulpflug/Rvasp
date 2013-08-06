\name{read.calculations}
\alias{read.calculations}
\title{Reads severaly calculations}
\usage{
  read.calculations(folders, name = "calculations",
    calculations = list(), chgcarname = "",
    bandsxmlname = "", dosxmlname = "", update = F)
}
\arguments{
  \item{folders}{list/vector of folders containing a
  calculation each}

  \item{name}{name}

  \item{calculations}{provide an existing calculations
  object for updating or extending purpose}

  \item{chgcarname}{if name is provided, calls
  \code{\link{read.chgcar}}}

  \item{bandsxmlname}{if name is provided, calls
  \code{\link{read.bandsdata}}}

  \item{dosxmlname}{if name is provided, calls
  \code{\link{read.dosdata}}}

  \item{update}{if \code{TRUE} will override existing
  folders in calculations object provided by
  \code{calculations}}
}
\description{
  \code{read.calculations} reads severaly calculations.
  Needed folder structure
  /FolderYouProvide/Parameter_Value/POSCAR For two
  parameters:
  /FolderYouProvide/Parameter_Value+Parameter2_Value/POSCAR
  Will read CONTCAR / POSCAR / (total) Energy (from
  OSZICAR) To read CHGCAR / Bands / Dos provide CHGCAR-type
  filename / vasprun.xml / vasprun.xml
}

