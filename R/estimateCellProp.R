estimateCellProp <-function(userdata,refdata="FlowSorted.Blood.450k",cellTypes=NULL,nonnegative = TRUE,nProbes=50,normalize=TRUE,refplot=FALSE){
probeSelect="both"
  if(!is(userdata, "rgDataSet") & !is(userdata, "methDataSet") & !is.matrix(userdata) &
     !is(userdata, "MethylSet") & !is(userdata, "RGChannelSet")){
    stop("[predSex] The input must be a rgDataSet, a methDataSet or a matrix")}
    if(is(userdata, "rgDataSet")){
      userdata=getmeth(userdata)}else if(is(userdata, "RGChannelSet")){
      userdata=preprocessRaw(userdata)}else if(is.matrix(userdata)){
      if(min(userdata,na.rm=TRUE)<0 | max(userdata,na.rm=TRUE)>1){stop("[estimateCellProp] Methylation data should within [0,1]")}
      if(normalize){normalize=FALSE; message("rgDataSet or methDataSet required for normalization")}
    }

#edit
if(refdata=="FlowSorted.Blood.450k"){
  library(refdata, character.only = TRUE)
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
#CellType: Bcell  CD4T  CD8T   Eos  Gran  Mono   Neu    NK  PBMC   WBC, 6 each
  flag=refdata$CellType %in% c("PBMC","WBC","Eos","Neu")
  refdata=refdata[,!flag]
}else if(refdata=="FlowSorted.DLPFC.450k"){
  library(refdata, character.only = TRUE)
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
#CellType: NeuN_neg, NeuN_pos, 29 each
}else if(refdata=="FlowSorted.CordBlood.450k"){
  library(refdata, character.only = TRUE)
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
#CellType: Bcell CD4T  CD8T Gran Mono NK nRBC WholeBlood
#          15    15    14   12   15  14   4     15
 flag=refdata$CellType %in% c("WholeBlood")
 refdata=refdata[,!flag]
}else if(refdata=="FlowSorted.CordBloodNorway.450k"){
  library(refdata, character.only = TRUE)
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
#Bcell  CD4T  CD8T  Gran  Mono    NK   WBC, 11 each
 flag=refdata$CellType %in% c("WBC")
 refdata=refdata[,!flag]
}else if(refdata=="FlowSorted.Blood.EPIC"){
  library(refdata, character.only = TRUE)
#  library(ExperimentHub)
  hub <- ExperimentHub()
  query(hub, "FlowSorted.Blood.EPIC")
  FlowSorted.Blood.EPIC <- hub[["EH1136"]]
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
  flag=refdata$CellType %in% c("MIX")
  refdata=refdata[,!flag]
#need to check
}else if(refdata=="FlowSorted.CordBloodCombined.450k"){
  library(refdata, character.only = TRUE)
  hub <- ExperimentHub()
  query(hub, "FlowSorted.CordBloodCombined.450k")
  FlowSorted.CordBloodCombined.450k <- hub[["EH2256"]]
  refdata=get(refdata)
  refdata=preprocessRaw(refdata)
#table(refdata$CellType)
#Bcell  CD4T  CD8T  Gran  Mono    NK  nRBC   WBC
#   42    41    33    43    48    45    11    26
 flag=refdata$CellType %in% c("WBC")
 refdata=refdata[,!flag]
}

    if(!is.null(cellTypes)) {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(colData(refdata)$CellType %in% cellTypes)
        refdata <- refdata[,keep]
    }
    pd <- as.data.frame(colData(refdata))
    pd$CellType <- factor(pd$CellType)

if(sum(colnames(userdata) %in% colnames(refdata)) > 0){
   tmp=colnames(userdata)[colnames(userdata) %in% colnames(refdata)]
   stop("The following column names in user data are the same with refdata: ",tmp,"\nPlease change to different names")}

if(refplot){
#distribution of intensity data
    jpeg(filename="refdata_distribution.jpg",width=700, height=1000,quality = 100)
    par(mfrow=c(2,1))
    color=as.numeric(pd$CellType)
    tmp=rbind(assays(refdata)$Meth,assays(refdata)$Unmeth);rownames(tmp)=1:nrow(tmp)
    multifreqpoly(tmp,col=color,legend="",
    cex.main=1.5,main="Reference data intensity value distribution",xlab="Intensity value")
    ntype=length(levels(pd$CellType))
    legend("top",legend=levels(pd$CellType),lty=1,lwd=3,col=1:ntype,text.col=1:ntype,bty="n")
    rm(tmp)
#distribution of methylation data
    color=as.numeric(pd$CellType)
    multifreqpoly(assays(refdata)$Meth/(assays(refdata)$Meth+assays(refdata)$Unmeth+100),col=color,legend="",
    cex.main=1.5,main="Reference data Methylation value distribution",xlab="Methylation value")
    ntype=length(levels(pd$CellType))
    legend("top",legend=levels(pd$CellType),lty=1,lwd=3,col=1:ntype,text.col=1:ntype,bty="n")
    dev.off()
}

commonprobe=intersect(rownames(userdata),rownames(refdata))
userdata=userdata[commonprobe,]
refdata=refdata[commonprobe,]

#quantile normalization
if(normalize){
Meth=cbind(assays(refdata)$Meth,assays(userdata)$Meth)
Unmeth=cbind(assays(refdata)$Unmeth,assays(userdata)$Unmeth)
rname <- rownames(Meth)
cname <- colnames(Meth)
Meth  <-  preprocessCore::normalize.quantiles(Meth)
Unmeth  <-  preprocessCore::normalize.quantiles(Unmeth)
methy=Meth/(Meth+Unmeth+100)
rownames(methy)=rname
colnames(methy)=cname
ref=methy[,colnames(refdata)]
userdata=methy[,colnames(userdata)]
}else{
ref=getBeta(refdata)
if(is(userdata,"methDataSet")){userdata=getB(userdata)
}else if(is(userdata,"MethylSet")){userdata=getBeta(userdata)}
}

if(refplot){
#distribution of methylation data
    jpeg(filename="refdata_distribution_after_qc.jpg",width=700, height=500,quality = 100)
    color=as.numeric(pd$CellType)
    multifreqpoly(ref,col=color,legend="",
    cex.main=1.5,main="Reference data Methylation value distribution",xlab="Methylation value")
    ntype=length(levels(pd$CellType))
    legend("top",legend=levels(pd$CellType),lty=1,lwd=3,col=1:ntype,text.col=1:ntype,bty="n")
    dev.off()
}

#select a set of CpG probes best differentiating celltypes
#require(genefilter)
Ftest=genefilter::rowFtests(ref, pd$CellType)
ref=ref[!is.na(Ftest$p.value) & Ftest$p.value<1e-8,]

    tIndexes <- split(x=seq(length(pd$CellType)),f=factor(pd$CellType))
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(ref))
        x[i] <- 1
        genefilter::rowttests(ref, factor(x))
    })
    if (probeSelect == "any"){
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]
            c(rownames(yAny)[1:(nProbes*2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yUp <- y[order(y[,"dm"], decreasing=TRUE),]
            yDown <- y[order(y[,"dm"], decreasing=FALSE),]
            c(rownames(yUp)[1:nProbes], rownames(yDown)[1:nProbes])
        })
    }
    trainingProbes <- unique(unlist(probeList))
    trainingProbes=trainingProbes[!is.na(trainingProbes)]
    ref <- ref[trainingProbes,]
    userdata = userdata[trainingProbes,]
    
    refmean <- sapply(split(x=seq(length(pd$CellType)),f=factor(pd$CellType)), 
               function(i) rowMeans(ref[,i]))
    refmedian <- sapply(split(x=seq(length(pd$CellType)),f=factor(pd$CellType)), 
               function(i) Biobase::rowMedians(ref[,i]))
    rownames(refmedian)=rownames(ref)

if(refplot){
#if (!require("gplots")){stop("Please install gplots for methylation heatmap")}
jpeg("refdata_heatmap.jpg",width=700,height=700,quality=100)
gplots::heatmap.2(ref,labCol=pd$CellType,labRow="",col=colorRampPalette(c("blue", "white", "red"))(256),key=TRUE,trace="none")
dev.off()
}

#Estimate mean methylation level for each celltypes
if(FALSE)
{
    modelFix <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if(ncol(phenoDF) == 2) { # two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(ref))
    } else { # > 2 group solution

    modelBatch=NULL

    N <- dim(phenoDF)[1]
    phenoDF$y <- rep(0, N)
    xTest <- model.matrix(modelFix, phenoDF)
    sizeModel <- dim(xTest)[2]
    M <- dim(ref)[1]

    coefEsts <- matrix(NA, M, sizeModel)

    for(j in 1:M) {
        ii <- !is.na(ref[j,])
        phenoDF$y <- ref[j,]
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=phenoDF[ii,]))
                OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
            } else{OLS <- TRUE}

            if(OLS) {
                fit <- lm(modelFix, data=phenoDF[ii,])
                fitCoef <- fit$coef
            } else {
                fitCoef <- fit$coef$fixed
            }
            coefEsts[j,] <- fitCoef
        })
    }
    rownames(coefEsts) <- rownames(ref)
    colnames(coefEsts) <- names(fitCoef)
    }
}
coefEsts=refmean
#Estimate cell type proportions
Xmat = coefEsts
nCol = dim(Xmat)[2]
nSubj = dim(userdata)[2]

  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(userdata)
  colnames(mixCoef) = colnames(Xmat)

  if(nonnegative){
#   if(!require(quadprog))stop("Can not load package quadprog")

    Amat = diag(nCol)
    b0vec = rep(0,nCol)

    for(i in 1:nSubj){
      obs = which(!is.na(userdata[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = quadprog::solve.QP(Dmat, t(Xmat[obs,])%*%userdata[obs,i], Amat, b0vec)$sol
    }
  }else{
    for(i in 1:nSubj){
      obs = which(!is.na(userdata[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% userdata[obs,i])
    }
  }

  data.frame(Sample_Name=rownames(mixCoef),mixCoef)
}



