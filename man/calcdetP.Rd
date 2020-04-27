\name{calcdetP}
\alias{calcdetP}
\title{
To calculate detection P values
}
\description{
Calculate detection P values based on negtive internal control probes or out of the
 band (oob) probes
}
\usage{
calcdetP(rgSet,detPtype = "negative")
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}}
  \item{detPtype}{Calculate detection P values based on negtive internal control
 ("negative") probes or out of the band ("oob") probes}
}
\value{
An numerical matrix of detection P values, with row for CpGs and column for samples
}
\references{
Wanding Zhou et al.
SeSAMe: reducing artifactual detection of DNA  methylation by Infinium BeadChips
 in genomic deletions, Nucleic Acids Research, 2018
}

\author{
Zongli Xu
}
\examples{
\dontrun{
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
detp=calcdetP(rgSet,detPtype = "negative")
detp2=calcdetP(rgSet,detPtype = "oob")
}
}

