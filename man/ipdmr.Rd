\name{ipdmr}
\alias{ipdmr}
\title{
   Differentially methylated region finder using interval P values
}
\description{
   To identify differentially methylated regions using an interval P value method
}
\usage{
       ipdmr(data,include.all.sig.sites=TRUE,dist.cutoff=1000,bin.size=50,
             seed=0.01,region_plot=TRUE,mht_plot=TRUE,verbose=TRUE)
}
\arguments{
  \item{data}{A data frame from bed format file with colname name "chr","start",
   "end","p","probe", "chr" indicate chromosome (1,2,3,...,X,Y), "start" and "end" are chromosome 
   position, "p" is for P value and "probe" for name of CpGs}
  \item{include.all.sig.sites}{Whether to use CpG singletons in calculate of FDR}
  \item{dist.cutoff}{Maximum distance in base pair to combine adjacent DMRs, and 
         the maximum distance between CpGs where auto-correlation will be calculated}
  \item{bin.size}{bin size for autocorrelation calculation}
  \item{seed}{FDR significance threshold for initial selection of DMR regions}
  \item{region_plot}{Whether to draw regional plot}
  \item{mht_plot}{Whether to draw mahattan plot}
  \item{verbose}{Whether to output detailed information}
}
\details{
  The input should be a data frame with column names "chr","start", "end","p","probe", "chr"
, indicating chromosome, 
start position, end position, P value and probe name. The function will use a 
novel interval p value method to identify differentially methylated regions. DMR results
 will be stored in a file with name "resu_ipdmr.csv". If plot options were selected,
 two figure files will be generated: mht.jpg and region_plot.pdf.
}
\references{
Zongli Xu, Changchun Xie, Jack A. Taylor, Liang Niu, ipDMR: Identification
of differentially methyl-ated regions with interval p-values, in review
}
\author{
Liang Niu, Zongli Xu
}

\examples{
\dontrun{
dat=simubed()
names(dat)
ipdmr(data=dat,seed=0.05)
}}


