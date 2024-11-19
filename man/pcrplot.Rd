\name{pcrplot}
\alias{pcrplot}
\title{
Principal component regression plot
}
\description{
First, principal component analysis will be performed in the standadized input data
matrix (standadized for each row/CpG), and then the specified number of top
principal components (which explain most data variation) will be used to perform
linear regression with each specified variable. Regression P values will be
plotted for exploration of methylation data variance structure or identification 
of possible confounding variables in association analysis.
}
\usage{
pcrplot(beta, cov,npc=50)
}
\arguments{
  \item{beta}{A methylation beta value matrix with rows for probes and columns for
 samples. The input matrix should not contain any missing value.}
  \item{cov}{A data frame of covariates. Categorical variables should be converted
 to factors. The number of rows should equal to the number of columns in beta matrix}
  \item{npc}{The number of top ranked principal components to be plotted}
}
\value{
A jpeg figure "svdscreeplot.jpg" to show the variations explained by each principal
 component.

A jpeg figure "pcr_diag.jpg" to show association strength between principal
 components and covariates with cell colors indicating different levels of
 association P values.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 
BeadChip}. Nucleic Acids Research 2015
}
\author{
Zongli Xu
}
\examples{
\donttest{

if (require(minfiData)) {
mdat <- preprocessRaw(RGsetEx)
beta=getBeta(mdat, "Illumina")
group=pData(mdat)$Sample_Group
slide=factor(pData(mdat)$Slide)
cov=data.frame(group,slide)
pcrplot(beta,cov,npc=6)
}}
}
