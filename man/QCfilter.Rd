\name{QCfilter}
\alias{QCfilter}
\title{
Sample, CpG, low quality data filter.
}
\description{
Filter low quality samples or CpGs, outlier samples or user specified samples or CpGs. 
}
\usage{
QCfilter(mdat,qcinfo=NULL,detPthre=0.000001,nbthre=3,samplethre=0.05,CpGthre=0.05,
         bisulthre=NULL,outlier=FALSE,outid=NULL, outCpG=NULL,plot=FALSE,rmlowq=FALSE)
}
\arguments{
  \item{mdat}{An object of class \code{methDataSet} or beta value matrix.}
  \item{qcinfo}{An object outputed from function \code{QCinfo}}
  \item{detPthre}{Detection P value threshold to identify low quality data point}
  \item{nbthre}{Number of bead threshold to identify low quality data point}
  \item{samplethre}{Threshold to identify low quality samples, the percentage of low quality methylation data points across probes for each sample}
  \item{CpGthre}{Threshold to identify low quality probes, percentage of low quality methylation data points across samples for each probe}
  \item{bisulthre}{Threshold of bisulfite intensity for identification of low quality samples. By default, Mean - 3 x SD of sample bisufite control intensities will be used as the threshold.}
  \item{outlier}{If TRUE, outlier samples will be excluded.}
  \item{outid}{A list of user specified samples to be excluded.}
  \item{outCpG}{A list of user specified CpGs to be excluded.}
  \item{plot}{TRUE or FALSE, whether to produce quality checking plots.}
  \item{rmlowq}{If TRUE and mdat is a beta value matrix, low quality data points as defined by detpthre and nbthre will be removed and set as missing data.}
}
\value{
An same type object as input object after excluding low quality samples and CpGs, or low quality data points.

Figure "qc_sample.jpg": scatter plot for Percent of low quality data per sample and Average bisulfite conversion intensity

Figure "qc_CpG.jpg": histogram for Percent of low quality data per CpG.

Figure "freqpolygon_beta_beforeQC.jpg": distribution plot before filtering. 

Figure "freqpolygon_beta_afterQC.jpg": distribution plot after filtering. 
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Nucleic Acids Research 2015.
}
\author{
Zongli Xu
}

\examples{
\dontrun{
if (require(minfiData)) {
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
qc=QCinfo(rgSet)
#for methDataSet
mdat=getmeth(rgSet)
mdat1=QCfilter(mdat,qcinfo=qc,outlier=TRUE)

#for methylation beta matrix
beta=getB(mdat)
beta=QCfilter(beta,qcinfo=qc,outlier=TRUE,rmlowq=TRUE)
}}
}

