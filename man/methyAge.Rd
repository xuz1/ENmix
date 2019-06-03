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
Horvath S. \emph{DNA methylation age of human tissues and cell types}. Genome biology 2013 14:R115.
Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, et al. \emph{Genome-wide methylation profiles reveal quantitative views of human aging rates}. Molecular cell 2013 49:359-367.
Levine ME, Lu AT, Quach A, Chen BH, Assimes TL, Bandinelli S, et al. \emph{An epigenetic biomarker of aging for lifespan and healthspan}. Aging (Albany NY) 2018 10:573-591.
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

