\name{getCGinfo}
\alias{getCGinfo}
\title{
CpG probe annotation inforamtion
}
\description{
Extract CpG probe annotation inforamtion from an rgDataSet
}
\usage{
getCGinfo(rgSet, type="IandII")
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}}
  \item{type}{One of the following options "I","II","IandII","ctrl", 
    indicating type I, type II type I & II or control probes type}
}
\value{
An object of data frame class}
\author{
Zongli Xu
}
\examples{
\dontrun{
require(minfidata)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
cginfo=getCGinfo(rgSet,type="IandII")
}}

