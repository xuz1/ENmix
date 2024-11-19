\name{dupicc}
\alias{dupicc}
\title{
Evaluation of measurement reliability using duplicate samples
}
\description{
The function can be used to evaluate duplicate samples by calculating: 
1) centered/un-centered Pearson's correlation coefficient between duplicates; 
2) abosolute difference between duplicates;
3) ICC for each CpG probes using oneway or twoway model.  
}
\usage{
dupicc(dat,dupid,mvalue=FALSE,center=TRUE,nCores=2,qcflag=FALSE,qc=NULL,
       detPthre=0.05,nbthre=3,skipicc=FALSE,corfig=FALSE,model="oneway")
}
\arguments{
  \item{dat}{Methylation beta value matrix}
  \item{dupid}{A data frame with two variables, id1 and id2. The two ids 
              in each row indicate a duplicate pair. These ids should be
              the same with column names of the input methylation matrix}
  \item{mvalue}{If TRUE, the beta value will be converted to M value for 
             calculation of ICC}
  \item{center}{If TRUE, the methylation beta values will be centered 
             for each CpG before calculation of ICC or correlation}
  \item{nCores}{Number of cores will be used for calculation of ICC}
  \item{qcflag}{Whether to perform QC before calculation of ICC}
  \item{qc}{QC object from function QCinfo}
  \item{detPthre}{If qcflag=TRUE, the methylation values with detection P value
       higher than the threshold will be removed before calculation}
  \item{nbthre}{If qcflag=TRUE, the methylation values with number of bead smaller 
       than the threshold will be removed}
  \item{skipicc}{If TRUE, ICC calculation will be skipped}
  \item{corfig}{If TRUE, a figure will be generated to demonstrate correlations 
         within duplicates or within non-duplicates}
  \item{model}{Using "oneway" or "twoway" model to calculate ICC} 
}
\value{
icc:  a data frame containing ICC and P values for each probe

dupcor: a data frame containing Pearson's correlation and 
        averaged absolute difference between duplicates. 

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
dupidx=data.frame(id1=c("5723646052_R02C02","5723646052_R04C01","5723646052_R05C02"),
id2=c("5723646053_R04C02","5723646053_R05C02","5723646053_R06C02"))
iccresu<-dupicc(dat=beta,dupid=dupidx)
}}
}

