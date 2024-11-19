\name{methyAge}
\alias{methyAge}
\title{
Methylation Age estimator
}
\description{
To calculate Methylation Age using Hovath, Hannum or PhenoAge methods and 
pace of aging DunedinPACE.
}
\usage{
methyAge(beta,fastImputation=FALSE,normalize=TRUE,nCores=2)
}
\arguments{
  \item{beta}{Methylation beta value matrix with CpG names(row names) and
       sample ids(column names).}
  \item{fastImputation}{If "TRUE" reference methylation values will used
      for imputation, if "FALSE", KNN nearest neighbor method will be used.}
  \item{normalize}{TRUE or FALSE, if TRUE, Hovath modified BMIQ method will
       be used to perform normalization.}
  \item{nCores}{Number of cores will be used for normalization}
}
\value{
  A data frame with rows for sample and columns for types of methylation age.
}

\references{
Horvath S. \emph{DNA methylation age of human tissues and cell types}. Genome biology 2013 14:R115.

Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, et al. \emph{Genome-wide methylation profiles reveal quantitative views of human aging rates}. Molecular cell 2013 49:359-367.

Levine ME, Lu AT, Quach A, Chen BH, Assimes TL, Bandinelli S, et al.  \emph{An epigenetic biomarker of aging for lifespan and healthspan}. Aging (Albany NY) 2018 10:573-591.

Daniel W Belsky, Avshalom Caspi, David L Corcoran,et al. \emph{DunedinPACE, a DNA methylation biomarker of the pace of aging}. eLife, 2022
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
beta=getB(meth)
mage=methyAge(beta)
}
}
