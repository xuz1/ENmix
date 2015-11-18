\name{multifreqpoly}
\alias{multifreqpoly}
\title{
Frequency polygon plot to display data distribution.
}
\description{
Produce Frequency polygon plot for each column of a numeric data matrix.
}
\usage{
multifreqpoly(mat, nbreaks=100, col=1:ncol(mat), xlab="", 
              ylab="Frequency",legend = list(x = "top", fill=col,
              legend = if(is.null(colnames(mat))) paste(1:ncol(mat)) 
              else colnames(mat)),...)
}
\arguments{
  \item{mat}{A numeric matrix}
  \item{nbreaks}{The number of bins for frequency counting}
  \item{col}{Line plot color code, the length should be equal to the number of columns in mat}
  \item{xlab}{x-axis lable}
  \item{ylab}{y-axis lable}
  \item{legend}{A list of arguments that get passed to the function "legend"}
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
multifreqpoly(beta,col=rep("black",ncol(beta)))
}}
}
