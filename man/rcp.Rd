\name{rcp}
\alias{rcp}
\title{
Illumina methylation array probe type bias correction
}
\description{
Probe design type bias correction using Regression on Correlated Probes (RCP) method
}
\usage{
        rcp(mdat, dist=25, quantile.grid=seq(0.001,0.999,by=0.001), qcscore = NULL,
 nbthre=3, detPthre=0.000001)
}
\arguments{
  \item{mdat}{An object of class \code{methDataSet} or \code{MethylSet}.}
  \item{dist}{Maximum distance in base pair between type I and type II probe pairs 
for regression calibration}
  \item{quantile.grid}{Quantile grid used in linear regression}
  \item{qcscore}{Data quality infomation object, the output from function QCinfo.  
       If the object is provied, low quality data points as defined by detection 
       p value threshold (detPthre) or number of bead threshold (nbthre) will 
       be set as missing values.}
  \item{detPthre}{Detection P value threshold to define low qualitye data points}
  \item{nbthre}{Number of beads threshold to define low qualitye data points, 
       nbthre=3 in default.}
}
\details{
The function will first identify type I and type II probe pairs within a specified
 distance, and then perform linear regression calibration between the probe types.
 With the estimates the function will then adjust type II data using type I data 
 as references. 
}
\value{
A beta value matrix 
}
\references{
Liang Niu, Zongli Xu and Jack A. Taylor
\emph{RCP: a novel probe design bias correction method for Illumina Methylation 
BeadChip, Bioinformatics 2016}
}
\author{
Liang Niu, Zongli Xu
}
\examples{
\donttest{

if (require(minfiData)) {
#methDataSet as input
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
qc=QCinfo(rgSet)
mdat=preprocessENmix(rgSet,QCinfo=qc,nCores=6)
mdat=norm.quantile(mdat,method="quantile1")
beta=rcp(mdat)

#methylset as input
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),"extdata"),
 pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
qc=QCinfo(rgSet)
mdat=preprocessENmix(rgSet,QCinfo=qc,nCores=6)
mdat=norm.quantile(mdat,method="quantile1")
beta=rcp(mdat)
}}
}
