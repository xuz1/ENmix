\name{ComBat.mc}
\alias{ComBat.mc}
\title{
A multi-processor wrapper for ComBat method.
}
\description{
A multi-processor wrapper for ComBat method. ComBat is a method to adjust batch effect where the batch covariate is known. 
}
\usage{
ComBat.mc(dat,batch,nCores = 1,...)
}
\arguments{
  \item{dat}{A data matrix with column for samples and row for probe.}
  \item{batch}{Batch covariate (multiple batches allowed)}
  \item{nCores}{Number of cores will be used for computation}
  \item{...}{See \code{ComBat} in \code{sva} package for extra options}
}
\value{
A data matrix with the same dimension as input data, adjusted for batch effects. 
Warning: Values for multimodal distributed CpGs could be over-adjusted. 
}
\references{
Johnson, WE, Rabinovic, A, and Li, C (2007). 
\emph{Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127.}
}
\author{
Zongli Xu
}

\seealso{
See \code{ComBat} in \code{sva} package for details.
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
mdat=norm.quantile(mdat,method="quantile1")
beta=bmiq.mc(mdat,nCores=10)
batch=factor(pData(mdat)$Slide)
betaC=ComBat.mc(beta,batch,nCores=6,mod=NULL)
}}
}




