\name{getmeth}
\alias{getmeth}
\title{
Create a methDataSet
}
\description{
Create a methDataSet based on a rgDataset.
}
\usage{
getmeth(rgSet)
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}}
}
\value{
An object of class \code{methDataSet}
}
\author{
Zongli Xu
}
\examples{
\dontrun{
require(minfidata)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
meth=getmeth(rgSet)
}}

