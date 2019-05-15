\name{methyAge}
\alias{methyAge}
\title{
Calculate Methylation Age
}
\description{
Calculation of Methylation Age using Hovath, Hannum or PhenoAge methods
}
\usage{
methyAge(beta,type="all",fastImputation=FALSE,normalize=TRUE,nCores=2)
}
\arguments{
  \item{beta}{Methylation beta value matrix with CpG names(row names) and sample ids(column names).}
  \item{type}{Which method will be used for methylation age calculation: "hovath", "hannum", "phenoAge" or "all"}
  \item{fastImputation}{If "TRUE" reference methylation values will used for imputation, if "FALSE", KNN nearest neighbor method will be used.}
  \item{normalize}{TRUE or FALSE, if TRUE, Hovath modified BMIQ method will be used to perform normalization.}
  \item{nCores}{Number of cores will be used for normalization}
}
\value{
  A data frame with sample ids and estimates methylation age.
}

\references{
  Hovath 
  \emph{Title}
  Journal 
  Hannum 
  \emph{Title}
  Journal 
  PhenoAge 
  \emph{Title}
  Journal 
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
beta=getB(meth)
mage=methyAge(beta)
}
}

