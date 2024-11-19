\docType{class}
\name{rgDataSet-class}
\alias{rgDataSet-class}
\alias{rgDataSet}
%methods
\alias{show, rgDataSet-method}

\title{Class \code{"rgDataSet"}}
\description{A class for storing Illumina methylation array 
raw intensity data of two color channels, and probe annotation 
infomation.}
\usage{
rgDataSet(Red = new("matrix"), Green = new("matrix"),
    NBeads = new("matrix"),rowData=new("DataFrame"),ictrl= new("DataFrame"),...)
}

\arguments{
\item{Red}{A matrix of Red channel intensity values with row for methylation probes
     and column for samples}
\item{Green}{A matrix of Green channel intensity values with row for methylation probes
     and column for samples}
\item{NBeads}{A matrix contains the number of beads used to generate intensity values
      on the Red and Green channels.}
\item{rowData}{A dataframe contains probe annotation information}
\item{ictrl}{A dataframe contais detailed information for Illumina internal control probes}
\item{...}{other arguments for class SummarizedExperiment}
}
\value{
  An object of class \code{rgDataSet}
}

\examples{
showClass("rgDataSet")
}

