\name{QCfilter}
\alias{QCfilter}
\title{
Filter samples and probes with low data quality.
}
\description{
Filter samples and probes with low data quality measured based on detection P values and number of beads.
}
\usage{
QCfilter(mdat, qcinfo, detPthre = 0.05, nbthre = 3, 
	samplethre = 0.01, CpGthre = 0.05,bisulthre=1,plot=FALSE,
        outid=NULL, outCpG=NULL)
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}}
  \item{qcinfo}{An object outputted from \code{QCinfo}}
  \item{detPthre}{Detection P value threshold}
  \item{nbthre}{Number of bead threshold}
  \item{samplethre}{Threshold for filtering samples, the percentage of low quality methylation data points across probes for each sample}
  \item{CpGthre}{Threshold for filtering probes, for percentage of low quality methylation data points across samples for each probe}
  \item{bisulthre}{Threshold for bisulfite to exclude samples.}
  \item{plot}{Whether to plot total intensity and beta value density plot before and after QC.}
  \item{outid}{User specified sample list to be filtered out. Valid ids are selected from colnames(mdat).}
  \item{outCpG}{User specified CpG list to be filtered out. Valid ids are selected from rownames(mdat).}
}
\value{
An object of class \code{MethylSet} with low data quality samples and probes filtered. 

If option plot=TRUE, four figures will be produced:  density_total_intensity_beforeQC.jpg, density_total_beta_beforeQC.jpg,
density_total_intensity_afterQC.jpg, density_total_beta_afterQC.jpg.
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
mdat=bgcorrect(rgSet,bgParaEst="oob",nCores=6)
#filter low quality samples and probes
mdat1=QCfilter(mdat,qcinfo=qcscore)
}}
}
