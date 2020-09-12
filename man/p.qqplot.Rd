\name{p.qqplot}
\alias{p.qqplot}
\title{
P value Q-Q plot
}
\description{
P value Q-Q plot with optional confidence interval
}
\usage{
       p.qqplot(pvalues,outf="qq",outfmt="jpg",draw.conf=TRUE,
          conf.col="lightgray",conf.alpha=.95,pch=20,col="black",reducesize=0)
}
\arguments{
  \item{pvalues}{An numeric list of P values}
  \item{outf}{figure file name, default "qq.jpg"}
  \item{outfmt}{Output figure file format, "jpg" or "eps"}
  \item{draw.conf}{Whether to draw confidence interval of expected P values under NULL hypothesis}
  \item{conf.col}{Color code of confidence interval}
  \item{conf.alpha}{Confidence interval range, 0.95 in default}
  \item{pch}{Point type code}
  \item{col}{Point color code}
  \item{reducesize}{Reduce eps file size by skipping some densely packed data points,
        larger the value, smaller the file size}
}
\details{P value Q-Q plot with optional confidence interval
}
\author{
Zongli Xu
}

\examples{
\dontrun{
dat=simubed()
p.qqplot(pvalues=dat$p,draw.conf=TRUE,outf="qq_try",outfmt="jpg")
}}

