\name{mhtplot}
\alias{mhtplot}
\title{
   P value manhattan plot
}
\description{
P value manhattan plot
}
\usage{
        mhtplot(dat,color="bg",sigthre=NULL,markprobe=NULL,
                 markcolor="red", outf="mht.jpg")
}
\arguments{
  \item{dat}{A data frame from bed format file with colname name "V1","V2",
   "V3","V4","V5", V1 indicate chromosome (1,2,3,...,X,Y), V2 and V3 are chromosome 
   position (V3 is an optional variable), V4 is for P value and V5 for name of CpGs}
  \item{color}{Color scheme of manhattan plot, "bg" indicate "black and gray"}
  \item{sigthre}{P value of significant threshold line}
  \item{markprobe}{A list of CpGs to be marked out}
  \item{markcolor}{Color code for marked probe}
  \item{outf}{figure file name, default "mht.jpg"}
}
\details{Draw manhattan plot based on bed file format input file. Optional 
      and colored marks fo a subset of probes
}
\author{
Zongli Xu
}

