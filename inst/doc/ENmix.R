
library(ENmix)
#read in raw intensity data
sheet <- read.450k.sheet(file.path(find.package("minfiData"),
        "extdata"), pattern = "csv$")
rgSet <- read.450k.exp(targets = sheet,extended = TRUE)
#Control plots
plotCtrl(rgSet)
#QC info
qc<-QCinfo(rgSet)
#background correction
mdat=preprocessENmix(rgSet,bgParaEst="oob",nCores=6)
#filter out low quality samples and probes
mdat=QCfilter(mdat,qcinfo=qc, samplethre = 0.01, CpGthre = 0.05
    ,bisulthre=5000,plot=TRUE, outid=NULL, outCpG=NULL)
#Search for multimodal CpGs
#sample size in this example data is too small for this purpose!
#should not use beta matrix after ComBat analysis for this purpose!
beta=getBeta(mdat, "Illumina")
nmode=nmode.mc(beta, minN = 3,modedist=0.2, nCores = 6)
#between-array normalization
mdat=normalize.quantile.450k(mdat,method="quantile1")
#probe type bias correction
beta=bmiq.mc(mdat,nCores=6)
# Principal component regression analysis plot
cov=data.frame(group=pData(mdat)$Sample_Group,
        slide=factor(pData(mdat)$Slide))
pcrplot(beta,cov,npc=6)
#batch correction
batch=factor(pData(mdat)$Slide)
betaC=ComBat.mc(beta,batch,nCores=6,mod=NULL)

