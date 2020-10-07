\name{QCinfo}
\alias{QCinfo}
\title{
Extract QC information
}
\description{
Extract information for data quanlity control: detection P values, number of
 beads and averaged bisulfite conversion intensity. 
The function can also identify low quality samples and probes, as well as outlier
 samples in total intensity or beta value distribution.
}
\usage{
QCinfo(rgSet, detPthre=0.000001, detPtype="negative", nbthre=3, samplethre=0.05,
       CpGthre=0.05, bisulthre=NULL, outlier=TRUE, distplot=TRUE)
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}, or \code{RGChannelSetExtended}} 
  \item{detPthre}{Detection P value threshold to identify low quality data point}
  \item{detPtype}{Calculate detection P values based on negtive internal control
 ("negative") probes or out of the band ("oob") probes}
  \item{nbthre}{Number of bead threshold to identify data point of low quality}
  \item{samplethre}{Threshold to identify samples with low data quality, the percentage of
 low quality methylation data points across probes for each sample}
  \item{CpGthre}{Threshold to identify probes with low data quality, percentage of low 
quality methylation data points across samples for each probe}
  \item{bisulthre}{Threshold of bisulfite intensity for identification of low 
quality samples. By default, Mean - 3 x SD of sample bisufite control 
intensities will be used as a threshold.}
  \item{outlier}{If TRUE, outlier samples in total intensity or beta value 
distribution will be idenfied and classified as bad samples.}
  \item{distplot}{TRUE or FALSE, whether to produce beta value distribution 
plots before and after QC.}
}
\value{
detP:  a matrix of detection P values

nbead: a matrix for number of beads

bisul:  a vector of averaged intensities for bisulfite conversion controls per sample

badsample: a list of low quality or outlier samples

badCpG:  a list of low quality CpGs

outlier_sample:  a list of outlier samples in methylation beta value or totol 
                intensity distribution.

Figure "qc_sample.jpg": scatter plot of Percent of low quality data per sample 
vs. Average bisulfite conversion intensity

Figure "qc_CpG.jpg": histogram for Percent of low quality data per CpG.

Figure "freqpolygon_beta_beforeQC.jpg": distribution plot of input data, 
   samples colored in red are "bad" samples, list in badsample, including samples 
   with low data quality and outlier in methylaiton beta value or total intensity.

Figure "freqpolygon_beta_afterQC.jpg": distribution plot input data after filtering
   "bad" samples.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450
 BeadChip}. Nucleic Acids Research 2015.
}
\author{
Zongli Xu
}

\examples{
\donttest{

if (require(minfiData)) {
#rgDataSet as input
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
qc=QCinfo(rgSet)

#RGChannelSetExtended as input
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),"extdata"),
 pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
qc<-QCinfo(rgSet)
}}
}
