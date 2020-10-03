\name{combp}
\alias{combp}
\title{
Identification of differentially methylated regions
}
\description{
To identify differentially methylated regions using a 
modified comb-p method
}
\usage{
       combp(data,dist.cutoff=1000,bin.size=310,seed=0.01,
             region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)
}
\arguments{
  \item{data}{A data frame with colname name "chr","start",
   "end","p" and "probe", indicating chromosome (1,2,3,...,X,Y),  
   chromosome start and end position, P value and probe names}
  \item{dist.cutoff}{Maximum distance in base pair to combine adjacent DMRs}
  \item{bin.size}{bin size for autocorrelation calculation}
  \item{seed}{FDR significance threshold for initial selection of DMR region}
  \item{region_plot}{If TRUE, regional plots will be generated}
  \item{mht_plot}{If TRUE, mahattan plot will be generated}
  \item{nCores}{Number of computer cores will be used in calculation}
  \item{verbose}{If TRUE, detailed running information will be printed}
}
\details{
The input should be a data frame with column names "chr","start", "end","p",
and "probe", indicating chromosome number, start position, end position, 
P value and probe name. The function use 
a modified comb-p method to identify differentially methylated regions. DMR 
results will be stored in a file with name resu_combp.csv. If plot options 
were selected, two figure files will be generated: mht.jpg and region_plot.pdf.
}
\references{
Pedersen BS1, Schwartz DA, Yang IV, Kechris KJ. Comb-p: software for combining, 
analyzing, grouping and correcting spatially correlated P-values. Bioinfomatics 2012

Zongli Xu, Changchun Xie, Jack A. Taylor, Liang Niu, ipDMR: Identification 
of differentially methyl-ated regions with interval p-values, Bioinformatics 2020

}
\author{
Liang Niu, Zongli Xu
}

\examples{
\dontrun{
dat=simubed()
names(dat)
combp(data=dat,seed=0.1)
}}



