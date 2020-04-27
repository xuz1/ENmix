\name{plotCtrl}
\alias{plotCtrl}
\title{
Internal control plot
}
\description{
The function will generate a series of internal control plots, which are similar to 
the plots from Illumina genomic studio software. Users should refer to the software's
online guide to interpret these figures. These figures can
 be used to check data quality and experimental procedures. 
}
\usage{
plotCtrl(rgSet,IDorder=NULL)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{IDorder}{A list of sample ids in the order user specified. The list can
 be a subset of the samples in input dataset. If an id list is provided, all
 plots will be produced in the order of the list. The default is NULL.}
}
\value{
A set of jpeg figures.
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
if(FALSE){
if (require(minfiData)) {
pinfo=pData(RGsetEx)
IDorder=rownames(pinfo)[order(pinfo$Slide,pinfo$Array)]
plotCtrl(RGsetEx,IDorder)
}}
}

