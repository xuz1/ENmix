\name{bmiq.mc}
\alias{bmiq.mc}
\title{
A multi-processor wrapper of BMIQ method
}
\description{
A multi-processor wrapper of BMIQ method. BMIQ is an intra-sample normalization procedure to correct the bias of Infinium 2 probe methylation beta values.
}
\usage{
bmiq.mc(mdat, nCores = 1,...)
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}.}
  \item{nCores}{Number of cores used for computation.}
  \item{...}{See BMIQ in R package wateRmelon for more options.}
}
\value{
A data matrix of Methylation beta value.
}
\references{
Teschendorff AE et. al. 
\emph{A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data}. Bioinformatics. 2013
}
\author{
Zongli Xu
}
\seealso{
See \code{BMIQ} in R package \code{watermelon} for model details
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
mdatq1=norm.quantile(mdat,method="quantile1")
beta=bmiq.mc(mdatq1,nCores=10)
}}
}


