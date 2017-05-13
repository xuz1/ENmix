\name{mpreprocess}
\alias{mpreprocess}
\title{
A pipeline to perform background correction, dye bias correction, inter-array normalization and probe type bias correction for HumanMethylation 450 and MethylationEPIC BeadChip data.
}
\description{
Fuction mpreprocess is a pipeline to easy preprocess steps for Illumina DNA methylation BeadChip. It comprehensively removes background noise and correct bias due to array design, including background correction, dye bias correction, inter-array normalization, probe type bias correction. It will also identify and exclude low quality samples and probes, remove outlier values, and perform imputation.
}
\usage{
    mpreprocess(rgSet,nCores=1,bgParaEst="oob",dyeCorr="RELIC",
                qc=FALSE,qnorm=TRUE,qmethod="quantile1",
                foutlier=TRUE,rmcr=FALSE,impute=FALSE)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSetExtended}, \code{RGChannelSet} or \code{MethylSet}.}
  \item{nCores}{Number of cores will be used for computation}
  \item{bgParaEst}{Method to estimate background normal distribution parameters. This must be one of the strings: "oob","est", or "neg".}
  \item{dyeCorr}{Dye bias correction, "mean": correction based on averaged red/green ratio; or "RELIC": correction with RELIC method; or "none": no dye bias correction. The default is RELIC}
  \item{qc}{If TRUE, QC will be performed. Low quality samples and CpGs will be excluded before background correction.}
  \item{qnorm}{If TRUE, inter-array quantile normalization will be performed.}
  \item{qmethod}{Quantile normalization method. This should be one of the following strings: "quantile1", "quantile2", or "quantile3". See details in function norm.quantile.}
  \item{foutlier}{If TRUE, outlier and low quality values will be filtered out.}
  \item{rmcr}{TRUE: excluded rows and columns with more than 5\% of missing values. FALSE is in default}
  \item{impute}{Whether to impute missing values. If TRUE, k-nearest neighbor's methods will be used for imputation. FALSE is in default.}
}
\details{
Fuction mpreprocess is a pipeline that perform methylaiton data preprocessing and quality controls using functions: preprocessENmix, norm.quantile, rcp, QCinfo and rm.outlier. More customized preprocessing steps can be achieved using the individual functions.
}
\value{
a methylation beta value matrix.
}
\references{
Zongli Xu, Liang Niu, Leping Li and Jack A. Taylor, ENmix: a novel background correction method for Illumina HumanMethylation450 BeadChip. Nucleic Acids Research 2015.

Zongli Xu, Sabine A. S. Langie, Patrick De Boever, Jack A. Taylor1 and Liang Niu, RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip, in review 2016

Liang Niu, Zongli Xu and Jack A. Taylor: RCP: a novel probe design bias correction method for Illumina Methylation BeadChip, Bioinformatics 2016
}
\author{
Zongli Xu
}
\seealso{
  Package \code{minfi} for classes \code{\linkS4class{RGChannelSet}} and \code{\linkS4class{MethylSet}}
}
\examples{
if(FALSE){
if (require(minfiData)) {
beta=mpreprocess(RGsetEx,nCores=6)

sheet <- read.metharray.sheet(file.path(find.package("minfiData"),"extdata"), pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
beta=mpreprocess(rgSet,nCores=6,qc=TRUE,foutlier=TRUE,rmcr=TRUE,impute=TRUE)
}}
}

