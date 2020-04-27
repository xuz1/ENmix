\name{predSex}
\alias{predSex}
\title{
  Estimating sample sex
}
\description{
  Estimating sample sex based on methylation data
}
\usage{
predSex(mdat, cutoff = -2)
}
\arguments{
  \item{mdat}{An object of class \code{MethDataSet} or \code{rgDataSet}.}
  \item{cutoff}{The difference in log2 total intensity between X and Y 
chromosomes}
}
\details{
  Estimation of sex is based on the difference of log2 median total 
intensity measures on the   X and Y chromosomes. 
}
\author{
Zongli Xu
}
\examples{
\dontrun{
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
sex=predSex(rgSet)
}
}
