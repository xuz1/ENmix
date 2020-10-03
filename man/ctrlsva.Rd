\name{ctrlsva}
\alias{ctrlsva}
\title{
Non-negative internal control surrogate variables
}
\description{
Surrogate variables derived based on intensity data for non-negative internal control 
probes.
}
\usage{
ctrlsva(rgSet,percvar=0.95,npc=1,flag=1)
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet} or \code{RGChannelSet}.}
  \item{percvar}{Minimum percentage of data variations can be explained by 
surrogate variables, range from 0 to 1,default is 0.95}
  \item{npc}{Number of surrogate variables, default is 1}
  \item{flag}{1: select number of surrogate variables based on argument
 \code{percvar}; 2: select number of surrogate variables based on argument \code{npc}}
}

\value{
The function will return an numerical matrix with columns indicating surrogate variables and rows
 corresponding to samples.  These variables can be used in association analysis to adjust for
experimental batch effects.
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
\dontrun{
if (require(minfiData)) {
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
sva<-ctrlsva(rgSet)
}}
}

