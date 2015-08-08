library(ENmix)
#read in data
sheet <- read.450k.sheet(file.path(find.package("minfiData"),
    "extdata"), pattern = "csv$")
rgSet <- read.450k.exp(targets = sheet, extended = TRUE)
#control plots
plotCtrl(rgSet)
#QC info
qc<-QCinfo(rgSet)
mraw <- preprocessRaw(rgSet)
beta<-getBeta(mraw, "Illumina")
#distribution plot
multifreqpoly(beta,main="Methylation Beta value distribution")
#Search for multimodal CpGs
#sample size in this example data is too small for this purpose!
bb=beta; bb[qc$detP>0.05 | qc$nbead<3]=NA #exclude low quality data first
nmode<-nmode.mc(bb, minN = 3, modedist=0.2, nCores = 6)
#background correction and dye bias correction
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE, nCores=6)
#exclude samples or CpGs with poor data quality
#user can also specify an extra list of CpG to be excluded, for example
outCpG = names(nmode)[nmode>1]
mdat<-QCfilter(mdat, qcinfo=qc, samplethre = 0.01, CpGthre = 0.05
    ,plot=TRUE, outid=NULL, outCpG=outCpG)
#inter-array normalization
mdat<-normalize.quantile.450k(mdat, method="quantile1")
#probe-type bias adjustment
beta<-bmiq.mc(mdat, nCores=6)
# Principal component regression analysis plot
cov<-data.frame(group=pData(mdat)$Sample_Group,
    slide=factor(pData(mdat)$Slide))
pcrplot(beta, cov, npc=6)
#filter out low quality and outlier values, remove rows and columns
#with too many missing value, and then do imputation
beta <- rm.outlier(beta,qcscore=qc,rmcr=TRUE,impute=TRUE)
#batch correction
#using M values instead of beta values maybe better at this step
batch<-factor(pData(mdat)[colnames(beta),]$Slide)
betaC<-ComBat.mc(beta, batch, nCores=6, mod=NULL)

                                                                                             
