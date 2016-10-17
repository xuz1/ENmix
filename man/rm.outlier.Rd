\name{rm.outlier}
\alias{rm.outlier}
\title{
Filtering out outlier and/or low quality values 
}
\description{
Setting outliers as missing value. Outlier was defined as value smaller than 3 times IQR from the lower quartile or larger than 3 times IQR from the upper quartile. If data quality information were provided, low quality data points will be set to missing first before looking for outliers. If specified, imputation will be performed using k-nearest neighbors method to impute all missing values. 
}
\usage{
rm.outlier(mat,byrow=TRUE,qcscore=NULL,detPthre=0.000001,nbthre=3,
           rmcr=FALSE,rthre=0.05,cthre=0.05,impute=FALSE,
           imputebyrow=TRUE,...)
}
\arguments{
  \item{mat}{An numeric matirx}
  \item{byrow}{TRUE: Looking for outliers row by row, or FALSE: column by column.}
  \item{qcscore}{If the data quality infomation (the output from function QCinfo) were provied, low quality data points as defined by detection p value threshold (detPthre) or number of bead threshold (nbthre) will be set to missing.} 
  \item{detPthre}{Detection P value threshold to define low qualitye data points, detPthre=0.000001 in default.}
  \item{nbthre}{Number of beads threshold define low qualitye data points, nbthre=3 in default.}
  \item{rmcr}{TRUE: excluded rows and columns with too many missing values as defined by rthre and cthre. FALSE is in default}
  \item{rthre}{Minimum of percentage of missing values for a row to be excluded}
  \item{cthre}{Minimum of percentage of missing values for a column to be excluded}
  \item{impute}{Whether to impute missing values. If TRUE, k-nearest neighbors methods will used for imputation. FALSE is in default. Warning: imputed values for multimodal distributed CpGs may not be correct.}
  \item{imputebyrow}{TRUE: impute missing values using similar values in row, or FALSE: in column}
  \item{...}{Arguments to be passed to the function impute.knn in R package "impute"}
}
\value{
An numeric matrix of same dimention as the input matrix.
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
qcscore<-QCinfo(rgSet)
mdat <- preprocessRaw(rgSet)
beta=getBeta(mdat, "Illumina")
#filter out outliers
b1=rm.outlier(beta)
#filter out low quality and outlier values
b2=rm.outlier(beta,qcscore=qcscore)
#filter out low quality and outlier values, remove rows and columns with too many missing values
b3=rm.outlier(beta,qcscore=qcscore,rmcr=TRUE)
#filter out low quality and outlier values, remove rows and columns with too many missing values, and then do imputation
b3=rm.outlier(beta,qcscore=qcscore,rmcr=TRUE,impute=TRUE)
}}
}
