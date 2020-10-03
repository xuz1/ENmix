\name{B2M}
\alias{B2M}
\title{
Converting methylation beta value to M value.
}
\description{
Convert methylation beta value to M value.
}
\usage{
B2M(x)
}
\arguments{
  \item{x}{An numeric matrix with values between 0 and 1}
}
\details{
Methylation beta value is calculated as beta=M/(M+U+a). M is methylated intensity,
U is unmethylated intensity, and a is a constant offset (by default , a=100).
M value is calculated as M=log2((M+a)/(U+a)). M or U is usually greater than 1000,
so a is negligible for most probes. if a=0, then M=log2(beta)/(1-beta).
}
\value{
A matrix of M values
}
\author{
Zongli Xu
}

\examples{
\dontrun{
if (require(minfiData)){
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat=getmeth(rgSet)
beta=getB(mdat,"Illumina")
m=B2M(beta)
}}
}

