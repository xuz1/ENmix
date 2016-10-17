\name{oxBS.MLE}
\alias{oxBS.MLE}
\title{
oxBS-MLE.
}
\description{
Find the Maximum Likelihood Estimate (MLE) of 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC)
using sequencing/array data from paired bisulfite and oxidative bisulfite treated DNA experiments.
}
\usage{
oxBS.MLE(beta.BS,beta.oxBS,N.BS,N.oxBS)
}
\arguments{
  \item{beta.BS}{A matrix of measurements (proportions of methylated signals in total signals) obtained from bisulfite (BS) experiments}
  \item{beta.oxBS}{A matrix of measurements obtained from oxidative bisulfite (oxBS) experiments}
  \item{N.BS}{A matrix of total signals from BS experiments}
  \item{N.oxBS}{A matrix of total signals from oxBS experiments}
}
\value{
A list with two elements:

\code{5mC}:  a matrix of estimated 5mC levels.

\code{5hmC}: a matrix for estimated 5hmC levels.
}
\details{
For all the inputs (\code{beta.BS}, \code{beta.oxBS}, \code{N.BS} and \code{N.oxBS}), the rows should be corresponding to CpG loci and the columns should be corresponding to samples.
The order of rows/columns in all four matrices should be consistent (otherwise \code{oxBS.MLE}  will stop with error messages).
Using a binomial model at each CpG locus in each sample, \code{oxBS.MLE} outputs a list with two matrices: a matrix of MLEs of 5mC levels and a matrix of MLEs of 5hmC levels.
The rows and columns of both matrices are consistent with the rows and columns of the input matrices.
For any CpG locus in any sample, if any of the four corresponding values (\code{beta.BS}, \code{beta.oxBS}, \code{N.BS} and \code{N.oxBS}) is \code{NA},
or \code{N.BS} is zero, or \code{N.oxBS} is zero, the MLE of both 5mC and 5hmC levels will be set as \code{NA}.
}
\references{
Zongli Xu, Jack A. Taylor, Yuet-Kin Leung, Shuk-Mei Ho and Liang Niu,
\emph{oxBS-MLE: An efficient method to estimate 5-methylcytosine and 5-hydroxymethylcytosine in paired bisulfite and oxidative bisulfite treated DNA}, under review.
}
\author{
Liang Niu and Zongli Xu
}
\examples{
# load example data
load(system.file("oxBS.MLE.RData",package="ENmix"))
# run oxBS.MLE
temp<-oxBS.MLE(beta.BS,beta.oxBS,N.BS,N.oxBS)
}

