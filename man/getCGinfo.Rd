\name{getCGinfo}
\alias{getCGinfo}
\title{
CpG probe annotation inforamtion
}
\description{
Extract CpG probe annotation inforamtion from a rgDataSet
}
\usage{
getCGinfo(rgSet, type="IandII")
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}}
  \item{type}{One of the following options "I","II","IandII","ctrl". "all" }
}
\value{
An object of class dataframe}
\author{
Zongli Xu
}
\examples{
\dontrun{
require(minfidata)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
info=getCGinfo(rgSet,type="all")
}}

