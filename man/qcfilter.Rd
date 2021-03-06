\name{qcfilter}
\alias{qcfilter}
\title{
Filtering out low quality and outlier data 
}
\description{
Outlier was defined as values smaller than 3
 times IQR from the lower quartile or greater than 3 times IQR from the upper
 quartile. If data quality information were provided, low quality data points
 will be set as missing data first before looking for outliers. All outliers
and low quality data will be set as miss in output matrix. If set imput=TRUE, 
imputation will be performed using k-nearest neighbors method to impute all 
missing values. 
}
\usage{
qcfilter(mat,qcscore=NULL,rmoutlier=TRUE,byrow=TRUE,detPthre=0.000001,nbthre=3,
     rmcr=FALSE,rthre=0.05,cthre=0.05,impute=FALSE,imputebyrow=TRUE,fastimpute=FALSE,...)
}
\arguments{
  \item{mat}{An numeric matirx containing methylation beta values}
  \item{qcscore}{If the data quality infomation (the output from function QCinfo)
 were provied, low quality data points as defined by detection p value threshold 
(detPthre) and number of bead threshold (nbthre) will be set as missing value.} 
  \item{rmoutlier}{if TRUE, outliers data points will be set as missing data NA.}
  \item{byrow}{TRUE: Looking for outliers row by row, or FALSE: column by column.}
  \item{detPthre}{Detection P value threshold to define low qualitye data points, 
detPthre=0.000001 in default.}
  \item{nbthre}{Number of beads threshold to define low qualitye data points, nbthre=3
 in default.}
  \item{rmcr}{TRUE: exclude rows and columns with too many missing values as 
defined by rthre and cthre. FALSE is in default}
  \item{rthre}{Minimum of percentage of missing values for a row to be excluded}
  \item{cthre}{Minimum of percentage of missing values for a column to be excluded}
  \item{impute}{If TRUE, k-nearest neighbors methods will used for imputation.}
  \item{imputebyrow}{TRUE: impute missing values using similar values in row, 
or FALSE: in column}
  \item{fastimpute}{If TRUE, probe median will be used for fast imputation.}
  \item{...}{Arguments to be passed to the function impute.knn in R package "impute"}
}
\value{
The output is an numeric matrix.
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
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
qc=QCinfo(rgSet)
mdat=preprocessENmix(rgSet,QCinfo=qc,nCores=6)
mdat=norm.quantile(mdat,method="quantile1")
beta=rcp(mdat)
#filter out outliers data points only
b1=qcfilter(beta)
#filter out low quality and outlier data points
b2=qcfilter(beta,qcscore=qc)
#filter out low quality and outlier values, remove rows and columns with
# too many missing values
b3=qcfilter(beta,qcscore=qc,rmcr=TRUE)
#filter out low quality and outlier values, remove rows and columns with
# too many missing values, and then do imputation
b3=qcfilter(beta,qcscore=qc,rmcr=TRUE,impute=TRUE)
}}}
