\name{ipdmr}
\alias{ipdmr}
\title{
   Differentially methylated regions
}
\description{
   To identify differentially methylated regions using an interval P value method
}
\usage{
       ipdmr(data,include.all.sig.sites=TRUE,dist.cutoff=1000,bin.size=50,
             seed=0.05,region_plot=TRUE,mht_plot=TRUE,verbose=TRUE)
}
\arguments{
  \item{data}{A data frame with colname name "chr","start",
   "end","p" and "probe", indicating chromosome (1,2,3,...,X,Y), chromosome 
   start and end position, P value and probe names}
  \item{include.all.sig.sites}{Whether to use CpG singletons in calculation of FDR}
  \item{dist.cutoff}{Maximum distance in base pair to combine adjacent DMRs, and 
         the maximum distance between CpGs where auto-correlation will be calculated}
  \item{bin.size}{bin size for autocorrelation calculation}
  \item{seed}{FDR threshold for initial selection of DMR regions}
  \item{region_plot}{If TRUE, regional plots will be produced for each DMR}
  \item{mht_plot}{If TRUE, a p-value mahattan plot with marked DMRs will be produced}
  \item{verbose}{Whether to output detailed information}
}
\details{
  The input should be a data frame with column names "chr","start", "end","p" and "probe",
indicating chromosome, start and end position, P value and probe name. The function will use a 
novel interval p value method to identify differentially methylated regions. DMR results
 will be stored in a file with name resu_ipdmr.csv. If plot options were selected,
 two figure files will be generated: mht.jpg and region_plot.pdf.
}
\references{
Zongli Xu, Changchun Xie, Jack A. Taylor, Liang Niu, ipDMR: Identification
of differentially methyl-ated regions with interval p-values, Bioinfomatics 2020
}
\author{
Liang Niu, Zongli Xu
}

\examples{
\donttest{
dat=simubed()
names(dat)
#seed=0.1 is only for demonstration purpose, it should be smaller than 0.05 or 0.01 in actual study.
ipdmr(data=dat,seed=0.1) #seed=0.1
}
}

