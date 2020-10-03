\name{oxBS.MLE}
\alias{oxBS.MLE}
\title{
oxBS-MLE.
}
\description{
Maximum Likelihood Estimate (MLE) of 5-methylcytosine (5mC) and
 5-hydroxymethylcytosine (5hmC) using sequencing/array data from
 paired bisulfite and oxidative bisulfite
 treated DNA experiments.
}
\usage{
oxBS.MLE(beta.BS,beta.oxBS,N.BS,N.oxBS)
}
\arguments{
  \item{beta.BS}{A matrix of methylation beta values (proportion of
  methylated sites estimated as methylated intensity over total 
  intensity) obtained from bisulfite (BS) experiments}
  \item{beta.oxBS}{A matrix of methylation beta values obtained from 
  oxidative bisulfite (oxBS) experiments}
  \item{N.BS}{A matrix of total signals (sum of methylated and unmethylated
              intensity values) from BS experiments}
  \item{N.oxBS}{A matrix of total signals from oxBS experiments}
}
\value{
The output is a list with two elements:
\code{5mC}:  a matrix of estimated 5mC levels.
\code{5hmC}: a matrix for estimated 5hmC levels.
}
\details{
For all the inputs (\code{beta.BS}, \code{beta.oxBS}, \code{N.BS} and
 \code{N.oxBS}), the rows should be corresponding to CpG loci and the
 columns should be corresponding to samples.
The row/column names in all four matrices should be the same. 
For a specific CpG of a sample, if any one of the four 
 values (\code{beta.BS}, \code{beta.oxBS}, \code{N.BS} and \code{N.oxBS}) 
is \code{NA}, or \code{N.BS} is zero, or \code{N.oxBS} is zero, the MLE of both 5mC 
and 5hmC levels will be set as \code{NA}.
}
\references{
Zongli Xu, Jack A. Taylor, Yuet-Kin Leung, Shuk-Mei Ho and Liang Niu,
\emph{oxBS-MLE: an efficient method to estimate 5-methylcytosine and 
5-hydroxymethylcytosine in paired bisulfite and oxidative bisulfite 
treated DNA}, Bioinformatics. 2016
}
\author{
Liang Niu and Zongli Xu
}
\examples{
# load example data
load(system.file("oxBS.MLE.RData",package="ENmix"))
# run oxBS.MLE
resu<-oxBS.MLE(beta.BS,beta.oxBS,N.BS,N.oxBS)
dim(resu[["5mC"]])
dim(resu[["5hmC"]])
}

