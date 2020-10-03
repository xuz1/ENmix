\name{preprocessENmix}
\alias{preprocessENmix}
\title{
The ENmix background correction
}
\description{
The ENmix background correction for HumanMethylation 450 and MethylationEPIC  BeadChip.
ENmix models methylation signal intensities with a flexible exponential-normal mixture
 distribution, and models background noise with a truncated normal distribution. ENmix
 will split BeadChip intensity data into 6 parts and separately model methylated and
 unmethylated intensities, 2 different color channels and 2 different probe types.
}
\usage{
preprocessENmix(rgSet, bgParaEst = "oob", dyeCorr="RELIC", QCinfo=NULL, exQCsample=TRUE,
                    exQCcpg=TRUE, exSample=NULL, exCpG=NULL, nCores = 2)
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}, \code{methDataSet}, 
     \code{RGChannelSetExtended}, \code{RGChannelSet} or \code{MethylSet}.}
  \item{bgParaEst}{Method to estimate background normal distribution parameters.
 Options are: "oob","est", or "neg".}
  \item{dyeCorr}{Dye bias correction method, "mean": correction based on averaged red/green
 ratio, or "RELIC": correction with RELIC method (default method), or "none": no dye bias correction.}
  \item{QCinfo}{If QCinfo object from function QCinfo() is provided, low quality 
samples (if exQCsample=TRUE) and CpGs (if exQCcpg=TRUE) will be excluded before
 background correction.}
  \item{exQCsample}{If TRUE, low quality samples listed in QCinfo will be excluded.}
  \item{exQCcpg}{If TRUE, low quality CpGs listed in QCinfo will be excluded.}
  \item{exSample}{User specified samples to be excluded before background correction}
  \item{exCpG}{User specified probes to be excluded before background correction}
  \item{nCores}{Number of cores will be used for computation}
}
\details{
By default, ENmix will use out-of-band Infinium I intensities ("oob") to estimate
 normal distribution parameters for modeling background noise. Option "est" will use
 combined methylated and unmethylated intensities to estimate background distribution
 parameters separately for each color channel and each probe type. Option "neg" will
 use 600 chip internal controls probes to estimate background distribution parameters.
 If rgSet if a \code{MethylSet}, then only option "est" can be selected. 
}
\value{
An object of class same with input data}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor, ENmix: a novel background
 correction method for Illumina HumanMethylation450 BeadChip. Nucleic Acids
 Research 2015.

Zongli Xu, Sabine A. S. Langie, Patrick De Boever, Jack A. Taylor1 and Liang 
Niu, RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip, 
BMC Genomics. 2017 
}
\author{
Zongli Xu and Liang Niu
}
\examples{
\dontrun{
if (require(minfiData)) {
#rgDataSet as input
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
#quality control information
qc<-QCinfo(rgSet)
#further excluding samples which are not in the qc$badsample list
ex_id=c("5723646053_R04C02")
#further excluding cpgs which are not in the qc$badCpG list
ex_cg=c("cg00000622", "cg00001245", "cg00001261")
mdat=preprocessENmix(rgSet,QCinfo=qc,exSample=ex_id,exCpG=ex_cg,nCores=6)


#RGChannelSet as input
mdat=preprocessENmix(RGsetEx,nCores=6)
}}
}

