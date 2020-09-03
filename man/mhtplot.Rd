\name{mhtplot}
\alias{mhtplot}
\title{
   P value manhattan plot
}
\description{
P value manhattan plot
}
\usage{
        mhtplot(probe=NULL,chr=NULL, pos=NULL, p=NULL,color="bg",sigthre=NULL,
     sigthre2=NULL,threlty=c(2,1),markprobe=NULL,markcolor="red", outf="mht",
     outfmt="jpg",reducesize=0)
}
\arguments{
  \item{probe}{probe name}
  \item{chr}{Chromosome, 1,2,...,22,X,Y}
  \item{pos}{Chromosome positions}
  \item{p}{P values}
  \item{color}{Color scheme of manhattan plot, "bg" indicate "black and gray"}
  \item{sigthre}{P value of significant threshold line}
  \item{sigthre2}{P value of second significant threshold line}
  \item{threlty}{Threshold line type, default is c(2,1)}
  \item{markprobe}{A list of CpGs to be marked out}
  \item{markcolor}{Color code for marked probe}
  \item{outf}{figure file name, default "mht.jpg"}
  \item{outfmt}{Output figure file format, "jpg" or "eps"}
  \item{reducesize}{Reduce eps file size by skipping some densely packed data points,
        larger the value, smaller the file size}
}
\details{Draw manhattan plot based on bed file format input file. Optional 
      and colored marks fo a subset of probes
}
\author{
Zongli Xu
}

