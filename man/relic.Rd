\name{relic}
\alias{relic}
\title{
RELIC dye bias correction method for Illumina HumanMethylation450 and MethylationEPIC  BeadChip
}
\description{REgression on Logarithm of Internal Control probes (RELIC) correct for dye bias on whole array by utilizing the intensity values of paired internal control probes that monitor the two color channels. 
}
\usage{
relic (mdat,at_red,cg_grn)
}
\arguments{
  \item{mdat}{An object of class \code{MethylSet}.}
  \item{at_red}{an intensity matrix for Illumina control probes "NORM_A" and "NORM_T"}
  \item{cg_grn}{an intensity matrix for Illumina control probes "NORM_C" and "NORM_G"}
}
\details{
The Illumina MethylationEPIC BeadChip contains 85 pairs of internal normalization control probes (name with prefix NORM_A, NORM_T, NORM_G or NORM_C), while its predecessor, Illumina HumanMethyl-ation450 BeadChip contains 93 pairs. RELIC first performs a regression on the logarithms of the intensity values of the normalization control probes to derive a quantitative relationship between red and green channels, and then uses the relation-ship to correct for dye-bias on intensity values for whole array. 
}
\value{
An object of class \code{MethylSet}
}
\references{
Zongli Xu, Sabine A. S. Langie, Patrick De Boever, Jack A. Taylor and Liang Niu, RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip, in review 2016
}
\author{
Zongli Xu and Liang Niu
}
\seealso{
  Package \code{preprocessENmix}
}
\examples{
if(FALSE){
if (require(minfiData)) {
    #background correction and dye bias correction
    mdat <- preprocessENmix(RGsetEx,bgParaEst="oob",nCores=6,dyeCorr =="RELIC")
    #dye bias correction only
    ctrls <- getProbeInfo(RGsetEx,type="Control")
    ctrls <- ctrls[ctrls$Address \%in\% featureNames(RGsetEx),]
    ctrl_r <- getRed(RGsetEx)[ctrls$Address,]
    ctrl_g <- getGreen(RGsetEx)[ctrls$Address,]
    CG.controls <- ctrls$Type \%in\% c("NORM_C","NORM_G")
    AT.controls <- ctrls$Type \%in\% c("NORM_A","NORM_T")
    cg_grn <- ctrl_g[CG.controls,]
    rownames(cg_grn) = ctrls$ExtendedType[CG.controls]
    at_red <- ctrl_r[AT.controls,]
    rownames(at_red) = ctrls$ExtendedType[AT.controls]
    mdat <- preprocessRaw(RGsetEx)
    mdat <- relic(mdat,at_red,cg_grn)
}}
}

