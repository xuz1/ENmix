\name{QCinfo}
\alias{QCinfo}
\title{
QC information.
}
\description{
Extract informations for data quanlity controls: detection P values and number of beads for each call of methylation beta value.
}
\usage{
QCinfo(rgSet)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSetExtended}.}
}
\value{
A list with 2 data matrices (for detection P value and number of beads) and a numeric array (average intensities for bisulfite conversion controls probes). 

Figure "qc_sample_1.jpg": scatter plot for Percent of low quality data per sample and Average bisulfite conversion intensity

Figure "qc_sample_2.jpg": histogram for Percent of low quality data per sample and Average bisulfite conversion intensity

Figure "qc_CpG.jpg": histogram for Percent of low quality data per CpG.
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
sheet <- read.450k.sheet(file.path(find.package("minfiData"),"extdata"), pattern = "csv$")
rgSet <- read.450k.exp(targets = sheet,extended = TRUE)
qcscore<-QCinfo(rgSet)
}}
}

