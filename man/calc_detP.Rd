\name{calc_detP}
\alias{calc_detP}
\title{
Calculate detection P values
}
\description{
Calculate detection P values based negtive internal control probes or out of the band (oob) probes
}
\usage{
calc_detP(rgSet,detPtype = "negative")
}
\arguments{
  \item{rgSet}{An object of class \code{rgDataSet}}
  \item{detPtype}{Calculate detection P values based on negtive internal control ("negative")
    probes or out of the band ("oob") probes}
}
\value{
An numerical matrix of detection P values, with row for CpGs and column for samples
}
\references{
Wanding Zhou et al.
eSAMe: reducing artifactual detection of DNA  methylation by Infinium BeadChips in genomic deletions
Nucleic Acids Research, 2018
}

\author{
Zongli Xu
}
\examples{
\dontrun{
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
detp=calc_detP(rgSet,detPtype = "negative")
detp2=calc_detP(rgSet,detPtype = "oob")
}
}

