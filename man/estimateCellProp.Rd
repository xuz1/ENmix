\name{estimateCellProp}
\alias{estimateCellProp}
\title{
Cell Proportion Estimation
}
\description{
Estimates relative proportion of underlying cell types in a sample based on reference
methylation data of pure cell types.  
}
\usage{
estimateCellProp(userdata,refdata="FlowSorted.Blood.450k",
                cellTypes=NULL,nonnegative = TRUE,nProbes=50,
                normalize=TRUE,refplot=FALSE)
}
\arguments{
  \item{userdata}{The input can be \code{rgDataSet},\code{methDataSet}\,code{MethylSet},\code{RGChannelSet} or methylation beta vaue matrix.}
  \item{refdata}{Reference data set will used: "FlowSorted.Blood.450k", "FlowSorted.DLPFC.450k", "FlowSorted.CordBlood.450k", "FlowSorted.CordBloodNorway.450k" or "FlowSorted.Blood.EPIC". }
  \item{cellTypes}{Specify cell type data in reference data will be used for deconvolution. if "NULL" all cell types will be used. }
  \item{normalize}{TRUE or FALSE, if TRUE, quantile normalization on methylated or unmethylated intensities will be performed.}
  \item{nonnegative}{TRUE or FALSE. If TRUE, the estimated proportions will be constrained to nonnegative values}
  \item{nProbes}{Number of best probes for each cell types will be used for the estimation.}
  \item{refplot}{TRUE or FALSE. IF TRUE, refdata distribution and heatmap will be plotted for inspection of reference dataset. }
}
\details{

This function use the method in Houseman et al (2012) to estimate cell type proportions using DNA methylation data. 

The following reference datasets can be used for the estimation. User should select a reference 
most resemble to user's data in tissue, age, and Illumina array type. 

FlowSorted.Blood.450k: consisting 450K methylation data for 60 blood samples from 6 male adults. 6 samples each for cell types: Bcell  CD4T  CD8T   Eos  Gran  Mono   Neu    NK  PBMC   WBC;  See Reinius et al. 2012 for details.

FlowSorted.CordBlood.450k: consisting 450k methylation data for 104 cord blood samples from 17 male and female individuals. Cell type (# samples) are: Bcell(15) CD4T(15)  CD8T(14) Gran(12) Mono(15) NK(14) nRBC(4) WholeBlood(15). See Bakulski et al. Epigenetics 2016 for details.

FlowSorted.CordBloodNorway.450k: consisting 450K methylation data for 77 cord blood samples from 11 individuals (6 girls and 5 boys). 11 samples for each of the cell types: Bcell  CD4T  CD8T  Gran  Mono    NK   WBC. See P Yousefi et al Environ. Mol. Mutagen 2015 for details.

FlowSorted.Blood.EPIC: consisting EPIC methylation data for 37 magnetic sorted blood cell references and 12 samples. See LA Salas et al. 2018 for details. 

FlowSorted.DLPFC.450k: consisting 450K methylation data for 58 brain tissue samples from 29 individuals. 15 females and 14 males, 6 Africans and 23 Caucasians, age between 13 to 79. 29 samples for each of the cell types: NeuN_neg and NeuN_pos. See Guintivano et al. 2013 for details.

}
\value{
  Matrix of composition estimates across all samples and cell types.
}

\references{
  EA Houseman, WP Accomando, DC Koestler, BC Christensen, CJ Marsit, HH
  Nelson, JK Wiencke and KT Kelsey.
  \emph{DNA methylation arrays as surrogate measures of cell mixture distribution.}
  BMC bioinformatics (2012) 13:86.
  doi:\href{http://www.dx.doi.org10.1186/1471-2105-13-86}{10.1186/1471-2105-13-86}.
}

\author{
Zongli Xu
}

\examples{
\dontrun{
require(minfidata)
path <- file.path(find.package("minfiData"),"extdata")
#based on rgDataset
rgSet <- readidat(path = path,recursive = TRUE)
celltype=estimateCellProp(userdata=rgSet,refdata="FlowSorted.Blood.450k",nonnegative = TRUE,normalize=TRUE)
#using methDataSet
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr="RELIC",
                      QCinfo=qc, exCpG=outCpG, nCores=6)
celltype=estimateCellProp(userdata=mdat,refdata="FlowSorted.Blood.450k",nonnegative = TRUE,normalize=TRUE)
mdat<-norm.quantile(mdat, method="quantile1")
#using beta value
beta<-rcp(mdat,qcscore=qc)
celltype=estimateCellProp(userdata=beta,refdata="FlowSorted.Blood.450k",nonnegative = TRUE)
}
}

