\name{rcp}
\alias{rcp}
\title{
Regression on Correlated Probes(RCP)
}
\description{
Probe design type bias correction using Regression on Correlated Probes (RCP) method
}
\usage{
        rcp(mdat, dist=25, quantile.grid=seq(0.001,0.999,by=0.001), qcscore = NULL, nbthre=3, detPthre=0.000001)
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}.}
  \item{dist}{Maximum distance in base pair between type I and type II probe pairs for regression calibration}
  \item{quantile.grid}{Quantile grid used in linear regression}
  \item{qcscore}{If the data quality infomation (the output from function QCinfo) is provied, low quality data points as defined by detection p value threshold (detPthre=0.000001) or number of bead threshold (nbthre=3) will be set to missing.}
  \item{detPthre}{Detection P value threshold to define low qualitye data points, detPthre=0.000001 in default.}
  \item{nbthre}{Number of beads threshold to define low qualitye data points, nbthre=3 in default.}

}
\details{
The function will first identify type I and type II probe pairs within specified distance, and then perform linear regression between the probe types to estimate regression coefficients. With the estimates the function will then calibrates type II data using type I data as references. 
}
\value{
A beta value matrix 
}
\references{
Liang Niu, Zongli Xu and Jack A. Taylor
\emph{RCP: a novel probe design bias correction method for Illumina Methylation BeadChip, Bioinformatics 2016}
}
\author{
Liang Niu, Zongli Xu
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat=preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6)
mdatq1=norm.quantile(mdat,method="quantile1")
beta=rcp(mdatq1)
}}
}

