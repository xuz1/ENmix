\name{pcrplot}
\alias{pcrplot}
\title{
Principal component regression plot
}
\description{
First, principal component analysis will be performed in the standadized input data
matrix (standadized for each row/CpG), and then the specified number of top
principal components (that explain most data variation) will be used to perform
linear regression with each specified variables. Regression P values will be
plotted for exploration of methylation data variance structure or identification 
of possible confounding variables for association analysis.
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
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Nucleic Acids Research 2015
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

