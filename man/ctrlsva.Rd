\name{ctrlsva}
\alias{ctrlsva}
\title{
Non-negative control surrogate variables
}
\description{
Surrogate variables derived from intensity data for non-negative internal control probes. These variables can be modeled in association analysis to adjust for experimental batch effects. 
}
\usage{
ctrlsva(rgSet,percvar=0.9,npc=1,flag=1)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{percvar}{Minimum percentage of data variations can be explained by surrogate variables, range from 0 to 1,default is 0.9}
  \item{npc}{Number of surrogate variables, default is 1}
  \item{flag}{1: select number of surrogate variables based on argument \code{percvar}; 2: select number of surrogate variables based on argument \code{npc}}
}
\value{
ctrlsva:  a matrix of surrogate variables (columns) with row corresponding to samples
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
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),"extdata"), pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
sva<-ctrlsva(rgSet)
}}
}

