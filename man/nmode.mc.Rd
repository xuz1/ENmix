\name{nmode.mc}
\alias{nmode.mc}
\title{
Estimating number of mode in methylaion data for each probe.
}
\description{
Due to SNPs in CpG probe region or other unknow factors, methylation beta values for some CpGs have multimodal distribution. This function is to identify this type of probes with obovious multimoal distribution. 
}
\usage{
nmode.mc(x, minN = 3, modedist=0.2, nCores = 1)
}
\arguments{
  \item{x}{A methylation beta value matrix with row for probes and column for samples.}
  \item{minN}{Minimum number of data points at each cluster}
  \item{modedist}{Minimum mode distance}
  \item{nCores}{Number of cores used for computation}
}
\details{
This function used an empirical approach to estimate number of mode in methylation beta value for each CpG probe. By default, the function requires the distance between modes have to be greater than 0.2 in methylation beta value, and each mode clusters should has at least 3 data points or 5\% of data points whichever is greater. 
}
\value{
A vector of integers
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Nucleic Acids Research 2015 
}
\author{
Zongli Xu
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat <- preprocessRaw(RGsetEx)
beta=getBeta(mdat, "Illumina")
nmode=nmode.mc(beta, minN = 3,modedist=0.2, nCores = 5)
}}
}

