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
     outfmt="jpg",reducesize=0,...)
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
  \item{markcolor}{Color code for marked probe, "red" in default}
  \item{outf}{figure file name, default "mht"}
  \item{outfmt}{Output figure file format, "jpg" or "eps"}
  \item{reducesize}{A positive interger, larger the value, smaller the eps file size.
      Smaller file size is achieved by skipping some densely packed data points}
  \item{...}{Arguments pass to plot}
}
\author{
Zongli Xu
}

\examples{
\donttest{

dat=simubed()

thre1=1E-100
dat$fdr=p.adjust(mrgd$P, method="BH")
if(sum(dat$fdr<0.05)>0){thre1=max(dat$p[dat$fdr<0.05])}

thre2=1E-7

mprobe=dat$probe[dat$p<=thre1]
mhtplot(probe=dat$probe,chr=dat$chr,pos=dat$start,p=dat$p,sigthre=thre1,sigthre2=thre2,
        markprobe=mprobe,outf="mht_try",outfmt="jpg")
}
}
