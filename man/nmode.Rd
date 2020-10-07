\name{nmode}
\alias{nmode}
\title{
Estimating number of mode for each row of data
}
\description{
Due to SNPs in CpG probe region or other unknow factors, methylation beta 
values for some CpGs have multimodal distribution. This function is to
 identify this type of probes (so called gap probes) with obovious 
multimoal distribution. 
}
\usage{
nmode(x, minN = 3, modedist=0.2, nCores = 1)
}
\arguments{
  \item{x}{A methylation beta value matrix with row for probes and column 
for samples.}
  \item{minN}{Minimum number of data points at each cluster}
  \item{modedist}{Minimum distance between adjacent modes}
  \item{nCores}{Number of cores used for computation}
}
\details{
This function uses an empirical approach to estimate number of modes in 
methylation beta value for each CpG probe. By default, the function 
requires the distance between modes have to be greater than 0.2 in 
methylation beta value, and each mode clusters should has at least 3 
data points or 5\% of data points whichever is greater. 
}
\value{
A vector of integers indicating number of modes. Gap probes are probes 
with number of mode greater than 1.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina 
HumanMethylation450 BeadChip}. Nucleic Acids Research 2015 
}
\author{
Zongli Xu
}
\examples{
\donttest{

if (require(minfiData)) {
mdat <- preprocessRaw(RGsetEx)
beta=getBeta(mdat, "Illumina")
nmode=nmode(beta, minN = 3,modedist=0.2, nCores = 5)

path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat <- getmeth(rgSet)
beta=getB(mdat)
nmode=nmode(beta, minN = 3,modedist=0.2, nCores = 5)
}}}
