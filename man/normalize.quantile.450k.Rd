\name{normalize.quantile.450k}
\alias{normalize.quantile.450k}
\title{
Quantile normalization.
}
\description{
Quantile normalization of methylation intensity data across samples for Illumina Infinium HumanMethylation 450 BeadChip.
}
\usage{
normalize.quantile.450k(mdat, method = "quantile1")
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}.}
  \item{method}{Quantile normalization method. This should be one of the following strings: "quantile1", "quantile2", or "quantile3".}
}
\details{
By default,  method = "quantile1" will separately quantile normalize Methylated or Unmethylated intensities for Infinium I or II probes. The "quantile2" will quantile normalize combined Methylated or Unmethylated intensities for Infinium I or II probes. The "quantile3" will quantile normalize combined Methylated or Unmethylated intensities for Infinium I and II probes together.
}
\value{
An object of class \code{MethylSet}.
}
\references{
Pidsley, R., CC, Y.W., Volta, M., Lunnon, K., Mill, J. and Schalkwyk, L.C. (2013) A data-driven approach to preprocessing Illumina 450K methylation array data. BMC genomics, 14, 293.
}

\author{
Zongli Xu
}

\examples{
if(FALSE){
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
mdatq1=normalize.quantile.450k(mdat,method="quantile1")
}}
}

