\name{getmeth}
\alias{getmeth}
\title{
Create a methDataSet
}
\description{
To create a methDataSet based on a rgDataset
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

\details{
CpG annotation information is incorporated in the output methDataSet 
object, and can be extracted using command \code{rowData().}
}
\author{
Zongli Xu
}
\examples{
\donttest{

require(minfiData)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
meth=getmeth(rgSet)
meth
cginfo=rowData(meth)
}
}
