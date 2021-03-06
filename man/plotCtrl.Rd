\name{plotCtrl}
\alias{plotCtrl}
\title{
Internal control plot
}
\description{
The function will generate a series of internal control plots that are similar to 
the plots from Illumina GenomeStudio software. Users should refer to GenomeStudio
online guide to interpret these figures. These figures can
 be used to check data quality and experimental procedures. 
}
\usage{
plotCtrl(rgSet,IDorder=NULL)
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet} or \code{RGChannelSet}.}
  \item{IDorder}{A list of sample ids in the order specified by user. The list can
 be a subset of sample ids in input dataset. If an id list is provided, all
 plots will be produced in the order of the list.}
}
\value{
A set of internal control figures in jpeg format.
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

if (require(minfiData)) {
#rgDataSet as input
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
plotCtrl(rgSet)

#RGChannelSet as input
pinfo=pData(RGsetEx)
IDorder=rownames(pinfo)[order(pinfo$Slide,pinfo$Array)]
plotCtrl(RGsetEx,IDorder)
}}
}
