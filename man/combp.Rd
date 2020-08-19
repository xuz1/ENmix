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
  \item{data}{A data frame from bed format file with colname name "V1","V2",
   "V3","V4","V5", V1 indicate chromosome (1,2,3,...,X,Y), V2 and V3 are chromosome 
   position, V4 is for P value and V5 for name of CpGs}
  \item{dist.cutoff}{Maximum distance in base pair to combine adjacent DMRs}
  \item{bin.size}{bin size for autocorrelation calculation}
  \item{seed}{FDR significance threshold for initial selection of DMR region}
  \item{region_plot}{Whether to draw regional plot}
  \item{mht_plot}{Whether to draw mahattan plot}
  \item{nCores}{Number of computer cores used in calculation}
  \item{verbose}{Whether to output detailed information}
}
\details{
The input should be a data frame with column name V1-V5, indicating chromosome,
 start position, end position, P values and probe names. The function will use 
a modified comb-p method to identify differentially methylated regions. DMR 
results will be stored in a file with name "resu_combp.csv". If plot options 
were selected, two figure files will be generated: mht.jpg and region_plot.pdf.
}
\references{
Pedersen BS1, Schwartz DA, Yang IV, Kechris KJ. Comb-p: software for combining, 
analyzing, grouping and correcting spatially correlated P-values. Bioinfomatics 2012

Zongli Xu, Changchun Xie, Jack A. Taylor, Liang Niu, ipDMR: Identification 
of differentially methyl-ated regions with interval p-values, in review

}
\author{
Liang Niu, Zongli Xu
}

