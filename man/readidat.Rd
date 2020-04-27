\name{readidat}
\alias{readidat}
\title{
Parsing IDAT files from Illumina methylation arrays.
}
\description{
Read in IDAT files and create a rgDataSet with probe annotation
}
\usage{
readidat(path = NULL,manifestfile=NULL,recursive = FALSE, verbose = FALSE)
}
\arguments{
  \item{path}{Directory where idat files are located}
  \item{manifestfile}{The file name of array manifestfile, which can be downloaded 
    from Illumina website. Bioconductor manifest package will be used if not provided}
  \item{recursive}{if TRUE, idat files in the subdirectories will also be read in}
  \item{verbose}{if TRUE, detailed running info will be printed on screen}
}
\value{
An object of class \code{rgDataSet}, 
}
\details{
Probe annotation info can be extracted using command \code{rowData}
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
}}

