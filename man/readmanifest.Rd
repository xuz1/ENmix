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
An object of dataframe caintaining probe annotation information
}
\author{
Zongli Xu
}
\examples{
\dontrun{

#for MethylationEPIC v1.0 B5
system("wget http://webdata.illumina.com.s3-website-us-east-1.amazonaws.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
system("unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
mf="infinium-methylationepic-v-1-0-b5-manifest-file.csv"
manifest=readmanifest(mf)

#for MethylationEPIC v1.0 B4
system("wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip")
system("unzip infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip")
mf="MethylationEPIC_v-1-0_B4.csv"
manifest=readmanifest(mf)

#for HumanMethylation450
system("wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv")
mf="HumanMethylation450_15017482_v1-2.csv"
manifest=readmanifest(mf)

}}

