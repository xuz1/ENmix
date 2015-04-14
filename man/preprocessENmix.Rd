\name{preprocessENmix}
\alias{preprocessENmix}
\title{
The ENmix background correction for HumanMethylation450k BeadChip
}
\description{
ENmix models methylation signal intensities with a flexible exponential-normal mixture distribution, and models background noise with a truncated normal distribution. ENmix will split 450k BeadChip intensity data into 6 parts and separately model methylated and unmethylated intensities, 2 different color channels and 2 different probe designs.
}
\usage{
preprocessENmix(rgSet, bgParaEst = "oob", dyeCorr=FALSE, nCores = 2)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSetExtended}, \code{RGChannelSet} or \code{MethylSet}.}
  \item{bgParaEst}{Optional method to estimate background normal distribution parameters. This must be one of the strings: "oob","est", or "neg".}
  \item{dyeCorr}{Dye bias correction: "TRUE" or "FALSE"}
  \item{nCores}{Number of cores will be used for computation}
}
\details{
By default, ENmix will use out-of-band Infinium I intensities ("oob") to estimate normal distribution parameters to model background noise. Option "est" will use combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type. Option "neg" will use 600 chip internal controls probes to estimate background distribution parameters. If rgSet if a \code{MethylSet}, then only option "est" can be selected. 
}
\value{
An object of class \code{MethylSet}
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor, ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip. In review.
}
\author{
Zongli Xu and Liang Niu
}
\seealso{
  Package \code{minfi} for classes \code{\linkS4class{RGChannelSet}} and \code{\linkS4class{MethylSet}}
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
}}
}

