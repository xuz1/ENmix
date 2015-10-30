\name{freqpoly}
\alias{freqpoly}
\title{
Frequency polygon plot 
}
\description{
Similar to histogram, frequency polygon plot can be used to display data distribution.
}
\usage{
freqpoly(mat, nbreaks=15, col="black", xlab="", ylab="Frequency", 
         type="l",append=FALSE,...)
}
\arguments{
  \item{mat}{A numeric vector}
  \item{nbreaks}{Number of bins for frequency counting}
  \item{col}{color code}
  \item{xlab}{x-axis lable}
  \item{ylab}{y-axis lable}
  \item{type}{character indicating the type of plotting; actually any of the 'type's as in 'plot.default'.}
  \item{append}{TRUE or FALSE, whether to create a new figure or append to the current figure.}
  \item{...}{Further arguments that get passed to the function "plot"}
}
\value{
Frequency polygon plot.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip}. Nucleic Acids Research 2015.
}
\author{
Zongli Xu
}

\examples{
if(FALSE){
if (require(minfiData)) {
mdat <- preprocessRaw(RGsetEx)
beta=getBeta(mdat, "Illumina")
freqpoly(beta[,1])
}}
}
