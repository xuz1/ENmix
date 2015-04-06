\name{plotCtrl}
\alias{plotCtrl}
\title{
Plot internal controls of 450K BeadChip.
}
\description{
Intensity data are ploted for all internal control probe types on the Illumina Infinium HumanMethylation450 BeadChip. These figures can be used to check data quality and experimental procedures. 
}
\usage{
plotCtrl(rgSet)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
}
\value{
A set of jpeg figures.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Under review.
}
\author{
Zongli Xu
}

\examples{
if(FALSE){
if (require(minfiData)) {
plotCtrl(RGsetEx)
}}
}

