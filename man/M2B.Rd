\name{M2B}
\alias{M2B}
\title{
Convert methylation M value to methylation beta value.
}
\description{
Convert methylation M value to methylation beta value.
}
\usage{
M2B(x)
}
\arguments{
  \item{x}{An numeric matrix}
}
\details{
Methylation beta value is calculated as beta=M/(M+U+a). M is methylated intensity,
U is unmethylated intensity, and a is a constant offset (by default , a=100).
M value is calculated as M=log2((M+a)/(U+a)). M or U is usually greater than 1000,
so a is negligible for most probes. if a=0, then beta=2^M/(2^M+1).
}

\value{
A matrix of methylation Beta values.
}
\author{
Zongli Xu
}


