\name{getBeta}
\alias{getBeta}
\title{
Extract Beta value.
}
\description{
Extract Methylation Beta value, Beta = Meth / (Meth + Unmeth + offset)
}
\usage{
getBeta(mdat,type="Illumina",offset=100)
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}.}
  \item{type}{type="Illumina" sets offset=100 as per Genome Studio.}
  \item{offset}{offset in calculating beta ratio}
}
\value{
beta:  a matrix of beta values
}
\author{
Zongli Xu
}


