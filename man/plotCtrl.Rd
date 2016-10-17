\name{plotCtrl}
\alias{plotCtrl}
\title{
Plot internal controls of 450K  or MethylationEPIC BeadChip.
}
\description{
Intensity data are ploted for all internal control probe types on the Illumina Infinium HumanMethylation450 or MethylationEPIC BeadChip. These figures can be used to check data quality and experimental procedures. 
}
\usage{
plotCtrl(rgSet,IDorder=NULL)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{IDorder}{A list of sample ids in the order user specified. The list can be a subset of the samples in input dataset. If an id list is provided, all plots will be produced in the order of the list. The default is NULL.}
}
\value{
A set of jpeg figures.
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
pinfo=pData(RGsetEx)
IDorder=rownames(pinfo)[order(pinfo$Slide,pinfo$Array)]
plotCtrl(RGsetEx,IDorder)
}}
}

