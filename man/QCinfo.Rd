\name{QCinfo}
\alias{QCinfo}
\title{
QC information.
}
\description{
Extract informations for data quanlity controls: detection P values and number of beads for each call of methylation beta value.
}
\usage{
QCinfo(rgSet, detPthre=0.05, nbthre=3, samplethre=0.01, CpGthre=0.05,
          bisulthre=NULL, outlier=TRUE, distplot=TRUE)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSetExtended}.}
  \item{detPthre}{Detection P value threshold to identify low quality data point}
  \item{nbthre}{Number of bead threshold to identify low quality data point}
  \item{samplethre}{Threshold to identify low quality samples, the percentage of low quality methylation data points across probes for each sample}
  \item{CpGthre}{Threshold to identify low quality probes, percentage of low quality methylation data points across samples for each probe}
  \item{bisulthre}{Threshold of bisulfite intensity for identification of low quality samples. By default, Mean - 3 x SD of sample bisufite control intensities will be used as the threshold.}
  \item{outlier}{If TRUE, outlier samples in total intensity or beta value distribution will be idenfied and classified as bad samples.}
  \item{distplot}{TRUE or FALSE, whether to produce beta value distribution plots before and after QC.}
}
\value{
detP:  a matrix of detection P values
nbead: a matrix for number of beads
bisul:  a vector of averaged intensities for bisulfite conversion controls
badsample: a list of low quality or outlier samples
badCpG:  a list of low quality CpGs
Figure "qc_sample.jpg": scatter plot for Percent of low quality data per sample and Average bisulfite conversion intensity
Figure "qc_CpG.jpg": histogram for Percent of low quality data per CpG.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Nucleic Acids Research 2015.
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

