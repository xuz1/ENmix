\name{rm.cgsuffix}
\alias{rm.cgsuffix}
\title{
Remove suffix from CpG names and combine duplicated CpGs
}
\description{
Remove suffix from CpG names for EPIC v2 BeadChips and combined values for duplicated CpGs
}
\usage{
rm.cgsuffix(datMeth)
}
\arguments{
  \item{datMeth}{A methlation data matrix with row names for CpG id}
}
\details{
Remove suffix from CpG names for EPIC v2 BeadChips and combined values for duplicated CpGs
}
\value{
A matrix with number of rows equal or less than input data.
}
\author{
Zongli Xu
}

\examples{
\donttest{

#beta matrix with row for CpGs and column for samples
beta2=rm.cgsuffix(beta)
}}
