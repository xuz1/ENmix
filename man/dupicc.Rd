\name{dupicc}
\alias{dupicc}
\title{
Methylation duplicates evaluation.
}
\description{
Calcuate ICC (consistency and twoway) and Pearson's correlation coefficient for each CpG probes. Calculate Pearson's correlation coefficient and methylation differences between duplicate samples.  
}
\usage{
dupicc(dat,dupid,mvalue=FALSE,center=FALSE,ncores=2,qcflag=FALSE,qc=NULL,detPthre=0.05,nbthre=3,skipicc=FALSE)
}
\arguments{
  \item{dat}{Methylation beta value matrix}
  \item{dupid}{A data frame with two variables, id1 and id2, for corresponding duplicates ids in methylation matrix}
  \item{mvalue}{Whether to convert to M value for calculation of ICC}
  \item{center}{Whether methylation beta values will be centered for calculation of correlation between duplicate samples}
  \item{ncores}{Number of core will be used for calculation of ICC}
  \item{qcflag}{Whether to perform QC before calculation of ICC}
  \item{qc}{QC object from ENmix package}
  \item{detPthre}{Detection P value threshold to identify low quality data point}
  \item{nbthre}{Number of bead threshold to identify low quality data point}
  \item{skipicc}{If TRUE, ICC calculation will be skipped}
}
\value{
icc:  a data frame contain ICC and P values for each probes

dupcor: a data frame contain Pearson's correlation and averaged absolute difference between duplicates. 

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
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),"extdata"), pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
mdat=preprocessRaw(rgSet)
beta=getBeta(mdat,"Illumina")
dupidx=data.frame(id1=c("5723646052_R02C02","5723646052_R04C01","5723646052_R05C02"),
id2=c("5723646053_R04C02","5723646053_R05C02","5723646053_R06C02"))
iccresu<-dupicc(dat=beta,dupid=dupidx)
}}
}

