\name{readmanifest}
\alias{readmanifest}
\title{
Parsing Illumina methylation arrays manifest file.
}
\description{
Parsing Illumina methylation arrays manifest file.
}
\usage{
readmanifest(file)
}
\arguments{
  \item{file}{Illumina methylation array manifest file, downloaded from Illuminal website}
}
\value{
An object of dataframe caintains probe annotation information
}
\author{
Zongli Xu
}
\examples{
\dontrun{
manifestfile="manifest *.csv file path"
manifest=readmanifest(manifestfile)
}}

