\name{repicc}
\alias{repicc}
\title{
Calculating intraclass correlation coefficient using replicate samples
}
\description{
The function can be used to calculate ICC for each CpG probe using 
balanced or unbalanced replicate samples. 
}
\usage{
repicc(dat,repid,mvalue=FALSE,nCores=2,qcflag=FALSE,qc=NULL,
       detPthre=0.05,nbthre=3)
}
\arguments{
  \item{dat}{Methylation beta value matrix}
  \item{repid}{A data frame with two variables, id and idx. id should be 
              the same with column name of "dat", idx is a variable
              to show the relationship between samples with same value for   
              samples from same individual.}
  \item{mvalue}{If TRUE, the beta value will be converted to M value for 
             calculation of ICC}
  \item{nCores}{Number of cores will be used for calculation of ICC}
  \item{qcflag}{Whether to perform QC before calculation of ICC}
  \item{qc}{QC object from function QCinfo}
  \item{detPthre}{If qcflag=TRUE, the methylation values with detection P value
       higher than the threshold will be removed before calculation}
  \item{nbthre}{If qcflag=TRUE, the methylation values with number of bead smaller}
}
\value{
   A data frame containing ICC for each probe
}
\references{
Zongli Xu, Jack A Taylor. \emph{Reliability of DNA methylation measures using 
Illumina methylation BeadChip}. Epigenetics 2020
}
\author{
Zongli Xu
}

\examples{
\donttest{
if (require(minfiData)){
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat=getmeth(rgSet)
beta=getB(mdat,"Illumina")
repid=data.frame(id=c("5723646052_R02C02","5723646052_R04C01","5723646052_R05C02","5723646053_R04C02","5723646053_R05C02","5723646053_R06C02"),idx=c(1,2,3,1,2,3))
iccresu<-repicc(dat=beta,repid=repid)
}}
}

