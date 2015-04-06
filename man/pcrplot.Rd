\name{pcrplot}
\alias{pcrplot}
\title{
Principal component regression plot
}
\description{
Principal component regression to explore methylation data variance structure and identifying possible confounding covariates. Principal components are derived using singular value decomposition method in beta value matrix. 
}
\usage{
pcrplot(beta, cov,npc=50)
}
\arguments{
  \item{beta}{A methylation beta value matrix with row for probes and column for samples.}
  \item{cov}{A data frame of covariates. Categorical variables should be converted to factors.}
  \item{npc}{The number of top principal components to plot}
}
\value{
A jpeg figure "svdscreeplot.jpg" to show the variations explained by each principal component.

A jpeg figure "pcr_diag.jpg" to show association strength between principal components and covariates with cell colors indicating different levels of association P values.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Under review
}
\author{
Zongli Xu
}
\examples{
if(FALSE){
if (require(minfiData)) {
mdat <- preprocessRaw(RGsetEx)
beta=getBeta(mdat, "Illumina")
group=pData(mdat)$Sample_Group
slide=factor(pData(mdat)$Slide)
cov=data.frame(group,slide)
pcrplot(beta,cov,npc=6)
}}
}

