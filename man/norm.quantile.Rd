\name{norm.quantile}
\alias{norm.quantile}
\title{
Quantile normalization.
}
\description{
Quantile normalization of methylation intensity data across samples for Illumina 
Infinium HumanMethylation 450 and MethylationEPIC BeadChip.
}
\usage{
norm.quantile(mdat, method = "quantile1")
}
\arguments{
  \item{mdat}{An object of class \code{methDataSet} or \code{MethylSet}.}
  \item{method}{Quantile normalization method: "quantile1", "quantile2", or "quantile3".}
}
\details{
By default,  method = "quantile1", which will separately quantile normalize Methylated or
 Unmethylated intensities for Infinium I or II probes. The "quantile2" will 
quantile normalize combined Methylated or Unmethylated intensities for
 Infinium I or II probes. The "quantile3" will quantile normalize combined
 Methylated or Unmethylated intensities for Infinium I and II probes together.
}
\value{
The output is an an object of class \code{methDataSet} or \code{MethylSet}.
}
\references{
Pidsley, R., CC, Y.W., Volta, M., Lunnon, K., Mill, J. and Schalkwyk, L.C.
 (2013) A data-driven approach to preprocessing Illumina 450K methylation
 array data. BMC genomics, 14, 293.
}

\author{
Zongli Xu
}

\examples{
\donttest{

#for methDataSet
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat<-preprocessENmix(rgSet, bgParaEst="oob",nCores=6)
mdatq<-norm.quantile(mdat, method="quantile1")

#for MethySet
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
mdatq=norm.quantile(mdat,method="quantile1")
}
}
}
