library(ENmix)
#read in data
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),
    "extdata"), pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet, extended = TRUE)
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
outCpG = names(nmode)[nmode>1]
#background correction and dye bias correction
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE,
                      QCinfo=qc, exCpG=outCpG, nCores=6)
#inter-array normalization
mdat<-normalize.quantile.450k(mdat, method="quantile1")
#probe-type bias adjustment
beta<-rcp(mdat)
# Principal component regression analysis plot
cov<-data.frame(group=pData(mdat)$Sample_Group,
    slide=factor(pData(mdat)$Slide))
pcrplot(beta, cov, npc=6)
#filter out low quality and outlier data points for each probe;
#rows and columns with too many missing value can be removed if specify;
#Do imputation to fill missing data if specify.
beta <- rm.outlier(beta,qcscore=qc,rmcr=TRUE,impute=TRUE)



