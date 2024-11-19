\docType{class}
\name{methDataSet-class}
\alias{methDataSet-class}
\alias{methDataSet}
%methods
\alias{show, methDataSet-method}

\title{Class \code{"methDataSet"}}
\description{A class for storing Illumina methylation array 
methylated and unmethylated intensity data, and CpG annotation 
infomation.
}
\usage{
methDataSet(Meth = new("matrix"), Unmeth = new("matrix"),
    rowData=new("DataFrame"),...)
}
\arguments{
\item{Meth}{A matrix of methylated intensity values with row for CpGs
     and column for samples}
\item{Unmeth}{A matrix of unmethylated intensity values with row for CpGs
     and column for samples}
\item{rowData}{A dataframe contains CpG annotation information}
\item{...}{Other arguments for class SummarizedExperiment}
}

\value{
  An object of class \code{methDataSet}
}
\details{
CpG annotation information is incorporated in the object, and can 
be extracted using command \code{rowData}
}
\examples{
showClass("methDataSet")
}



