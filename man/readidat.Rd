\name{readidat}
\alias{readidat}
\title{
Parsing IDAT files for Illumina methylation arrays .
}
\description{
Read in IDAT files and create a rgDataSet with probe annotation
}
\usage{
readidat(path = NULL,manifestfile=NULL,recursive = FALSE, verbose = FALSE)
}
\arguments{
  \item{path}{Directory where idat files are located}
  \item{manifestfile}{Array manifestfile, which can be downloaded 
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
\donttest{

#Illumina methylation array manifestfile
#for MethylationEPIC v1.0 B5
system("wget http://webdata.illumina.com.s3-website-us-east-1.amazonaws.com/downloads
/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
system("unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
mf="infinium-methylationepic-v-1-0-b5-manifest-file.csv"

#for MethylationEPIC v1.0 B4
system("wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles
/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip")
system("unzip infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip")
mf="MethylationEPIC_v-1-0_B4.csv"

#for HumanMethylation450
system("wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles
/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv")
mf="HumanMethylation450_15017482_v1-2.csv"
if(require(minfiData)){
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,manifestfile=mf,recursive = TRUE)

#without providing manifest file, corresponding R manifestfile package will be used
rgSet <- readidat(path = path,recursive = TRUE)
}}
}

