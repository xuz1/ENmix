\name{getB}
\alias{getB}
\title{
Extract methylation Beta values.
}
\description{
Extract Methylation Beta value, Beta = Meth / (Meth + Unmeth + offset)
}
\usage{
getB(mdat,type="Illumina",offset=100)
}
\arguments{
  \item{mdat}{An object of class \code{methDataSet} or  \code{MethylSet}.}
  \item{type}{type="Illumina" sets offset=100 as per Genome Studio software.}
  \item{offset}{Regularization factor in calculating beta ratio, 100 in default}
}
\value{
Methylation Beta value = Meth / (Meth + Unmeth + offset). Meth is methylated intensity matrix,
Unmeth is unmethylated intensity matrix.
}
\author{
Zongli Xu
}

\examples{
\dontrun{
if (require(minfiData)){
path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat=getmeth(rgSet)
beta=getB(mdat,"Illumina")
}}
}

