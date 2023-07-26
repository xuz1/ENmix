\name{rcp2}
\alias{rcp2}
\title{
Modified RCP method
}
\description{
Modified RCP method to normalize user's data to a list of reference values
}
\usage{
        rcp2(datMeth, reference, quantile.grid=seq(0.001,0.999,by=0.001))
}
\arguments{
  \item{datMeth}{A matrix with row for probes and column for samples}
  \item{reference}{A data frame with two columns, "cg" for CpG names and "goldstandard" for reference values}
  \item{quantile.grid}{Quantile grid used in linear regression}
}
\details{
The function will normalize user data to a reference value distribution based on 
common set of probe between user data and reference data.
}
\value{
A matrix with same dimension of user data
}
\references{
Liang Niu, Zongli Xu and Jack A. Taylor
\emph{RCP: a novel probe design bias correction method for Illumina Methylation 
BeadChip, Bioinformatics 2016}
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
reference=data.frame(cg=rownames(beta),goldstandard=beta[,3])
reference=reference[sample(1:nrow(reference),2000),]
beta2=rcp2(beta,reference)
}}
