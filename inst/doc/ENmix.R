library(ENmix)
#read in raw intensity data
sheet <- read.450k.sheet(file.path(find.package("minfiData"),
    "extdata"), pattern = "csv$")
rgSet <- read.450k.exp(targets = sheet, extended = TRUE)
#Control plots
pinfo=pData(rgSet)
IDorder=rownames(pinfo)[order(pinfo$Slide,pinfo$Array)]
plotCtrl(rgSet,IDorder)
#QC info
qc<-QCinfo(rgSet)
#Search for multimodal CpGs
#sample size in this example data is too small for this purpose!
#should not use beta matrix after ComBat analysis for this purpose!
mraw <- preprocessRaw(rgSet)
beta<-getBeta(mraw, "Illumina")
nmode<-nmode.mc(beta, minN = 3, modedist=0.2, nCores = 6)
#Frequency polygon plot to examining beta value distribution
anno=getAnnotation(rgSet)
beta1=beta[anno$Type=="I",]
beta2=beta[anno$Type=="II",]
library(geneplotter)
#comparisons between frequency polygon plot and density plot
#frequency polygon plot can more accurately reflect true distribution.
#and demostrate difference between Infinium I and II probes, and also easier to understand.
jpeg("dist.jpg",height=1200,width=800,quality=100)
par(mfrow=c(3,2))
multidensity(beta,main="Multidensity")
multifreqpoly(beta,main="Multifreqpoly")
multidensity(beta1,main="Multidensity: Infinium I")
multidensity(beta2,main="Multidensity: Infinium II")
multifreqpoly(beta1,main="Multifreqpoly: Infinium I")
multifreqpoly(beta2,main="Multifreqpoly: Infinium II")
dev.off()
#background correction
#user can also specify a list of bad CpGs to be excluded before background correction
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE, nCores=6)
#user specified CpG list to be excluded
outCpG = names(nmode)[nmode>1]
#Filter out Samples or CpGs with poor data quality
#users can also specify a list of samples or probes for exclusion, such as non-specific
#binding probes, probes affected by SNP, or multimodal distributed probes)
mdat<-QCfilter(mdat, qcinfo=qc, samplethre = 0.01, CpGthre = 0.05
    ,plot=TRUE, outid=NULL, outCpG=outCpG)
#between-array normalization
mdat<-normalize.quantile.450k(mdat, method="quantile1")
#probe type bias correction
beta<-bmiq.mc(mdat, nCores=6)
# Principal component regression analysis plot
cov<-data.frame(group=pData(mdat)$Sample_Group,
    slide=factor(pData(mdat)$Slide))
pcrplot(beta, cov, npc=6)
#filter out low quality and outlier values, remove rows and columns
#with too many missing value, and then do imputation
beta <- rm.outlier(beta,qcscore=qc,rmcr=TRUE,impute=TRUE)
#batch correction
batch<-factor(pData(mdat)[colnames(beta),]$Slide)
betaC<-ComBat.mc(beta, batch, nCores=6, mod=NULL)
                                                                                             
