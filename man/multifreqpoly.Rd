\name{multifreqpoly}
\alias{multifreqpoly}
\title{
Multiple frequency polygon plot
}
\description{
Produce Frequency polygon plot for each column of a numeric data matrix. 
Similar to multidensity function, the plot can be used to inspect data distribution
but with much faster speed.
}
\usage{
multifreqpoly(mat, nbreaks=100, col=1:ncol(mat), xlab="", 
              ylab="Frequency",legend = list(x = "top", fill=col,
              legend = if(is.null(colnames(mat))) paste(1:ncol(mat)) 
              else colnames(mat)),append=FALSE,...)
}
\arguments{
  \item{mat}{A numeric matrix}
  \item{nbreaks}{The number of bins for frequency counting}
  \item{col}{Line plot color code, the length should be equal to the number of
 columns in mat}
  \item{xlab}{x-axis lable}
  \item{ylab}{y-axis lable}
  \item{legend}{A list of arguments that get passed to the function "legend"}
  \item{append}{TRUE or FALSE, whether to create a new figure or append to the current figure.}
  \item{...}{Further arguments that get passed to the function "plot"}
}
\value{
Frequency polygon plots. 
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor,
\emph{ENmix: a novel background correction method for Illumina HumanMethylation450
 BeadChip}. Nucleic Acids Research 2015.
}
\author{
Zongli Xu
}

\examples{
\donttest{

x=matrix(rnorm(10000),nrow=2000,ncol=5)
multifreqpoly(x,nbreaks=15,legend=list(x="topright",fill=1:ncol(x),legend=paste("V",1:5,sep="")))

if (require(minfiData)) {
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mraw <- getmeth(rgSet)
beta<-getB(mraw)
jpeg("dist_raw.jpg")
multifreqpoly(beta,col=1:ncol(beta))
dev.off()
}
}}
