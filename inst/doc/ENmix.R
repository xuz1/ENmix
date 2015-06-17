library(ENmix)
#read in raw intensity data
sheet <- read.450k.sheet(file.path(find.package("minfiData"),
    "extdata"), pattern = "csv$")
rgSet <- read.450k.exp(targets = sheet, extended = TRUE)
#Control plots
plotCtrl(rgSet)
#QC info
qc<-QCinfo(rgSet)
#Search for multimodal CpGs
#sample size in this example data is too small for this purpose!
#should not use beta matrix after ComBat analysis for this purpose!
mraw <- preprocessRaw(rgSet)
beta<-getBeta(mraw, "Illumina")
nmode<-nmode.mc(beta, minN = 3, modedist=0.2, nCores = 6)
#Frequency polygon plot to examining beta value distribution
multifreqpoly(beta)
#background correction
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE, nCores=6)
#user specified CpG list to be excluded
outCpG = names(nmode)[nmode>1]
#exclude non-specific binding probes or probes affected by SNP et al.
#outCpG = unique(c(outCpG,non-specific bind probes,snp probes,...))
#filter out low quality samples and probes
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
#batch correction
batch<-factor(pData(mdat)$Slide)
betaC<-ComBat.mc(beta, batch, nCores=6, mod=NULL)


