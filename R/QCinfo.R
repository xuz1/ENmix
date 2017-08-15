QCinfo <- function(rgSet,detPthre=0.000001,nbthre=3,samplethre=0.05,CpGthre=0.05,
      bisulthre=NULL,outlier=TRUE,distplot=TRUE)
{
    ##number of bead
    if(!is(rgSet, "RGChannelSetExtended"))
      stop("ERROR: rgSet is not an object of class 'RGChannelSetExtended'")

    nb <- assays(rgSet)$NBeads
    typeI <- getProbeInfo(rgSet, type = "I")
    bcA <- nb[typeI$AddressA, ]
    bcB <- nb[typeI$AddressB, ]
    bc_I <- bcA
    flag <- bcA>bcB
    bc_I[flag] <- bcB[flag];
    typeII <- getProbeInfo(rgSet, type = "II")
    locusNames <- getManifestInfo(rgSet, "locusNames")
    nbead <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
    dimnames = list(locusNames, sampleNames(rgSet)))
    nbead[typeI$Name,]<-bc_I
    nbead[typeII$Name,]<-nb[typeII$AddressA,]

    ##detection P value
    detP<-detectionP(rgSet)

    ## bisulfite conversion internal controls
    ctrls<-getProbeInfo(rgSet,type="Control")
    ctrls=ctrls[ctrls$Address %in% rownames(rgSet),]
    ctrl_r <- getRed(rgSet)[ctrls$Address,]
    ctrl_g <- getGreen(rgSet)[ctrls$Address,]
    cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & (ctrls$ExtendedType 
       %in% c("BS Conversion I C1","BS Conversion I-C2","BS Conversion I-C3")),]
    I_green=colMeans(ctrl_g[cc$Address,])
    cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & (ctrls$ExtendedType
       %in% c("BS Conversion I-C4","BS Conversion I-C5","BS Conversion I-C6")),]
    I_red=colMeans(ctrl_r[cc$Address,])
    cc=ctrls[ctrls$Type %in% c("BISULFITE CONVERSION II"),]
    II_red=colMeans(ctrl_r[cc$Address,])
    bisul=(I_green+I_red+II_red)/3

    #threshold of bisulfite conversion control intensity
    if(is.null(bisulthre)){bisulthre=mean(bisul,na.rm=TRUE)-3*sd(bisul,na.rm=TRUE)}

    ##low quality samples
    qcmat <- nbead<nbthre | detP>detPthre
    badValuePerSample <- apply(qcmat,2,sum)/nrow(qcmat)
    flag <- badValuePerSample > samplethre | bisul < bisulthre
    cat(sum(flag)," samples with percentage of low quality CpG value greater
     than ",samplethre, " or bisulfite intensity less than ", bisulthre, "\n")
    badsample=colnames(qcmat)[flag]

    ##low quality CpGs
    qcmat <- qcmat[,!flag]
    NbadValuePerCpG <- apply(qcmat,1,sum)
    badValuePerCpG <- NbadValuePerCpG/ncol(qcmat)
    flag2 <- badValuePerCpG>CpGthre & NbadValuePerCpG>1
    cat(sum(flag2)," CpGs with percentage of low quality value greater than ",
      CpGthre,"\n")
    badCpG <- rownames(qcmat)[flag2]
    qcmat=qcmat[!flag2,]

    #plotting quality scores
    cat("Ploting qc_sample.jpg ...")
    jpeg(filename="qc_sample.jpg",width=1000,height=1000,quality = 100)
    color=rep("black",length(badValuePerSample))
    color[flag]="red"
    plot(badValuePerSample,bisul,xlab="Percent of low quality data per sample",
       ylab="Average bisulfite conversion intensity",cex=1.5,col=color,
       main=paste(length(badsample)," samples were classified as low quality samples"))
    abline(h=bisulthre,lty=2,col="red")
    abline(v=samplethre,lty=2,col="red")
    dev.off()
    cat("Done\n")

    cat("Ploting qc_CpG.jpg ...")
    jpeg(filename="qc_CpG.jpg",width=1000,height=1000,quality = 100)
    par(mfrow=c(2,1))
    hist(badValuePerCpG,breaks=1000,xlab="Percent of low quality data per CpG",
     main=paste(length(badCpG)," CpGs were classified as low quality CpGs"))
    abline(v=CpGthre,lty=2,col="red")
    hist(badValuePerCpG,breaks=1000,xlim=c(0,0.1),
       xlab="Percent of low quality data per CpG",main="Zoom in view")
    abline(v=CpGthre,lty=2,col="red")
    dev.off()
    cat("Done\n")

    #Identifying outlier samples
    if(outlier)
    {
    cat("Identifying ourlier samples based on beta or total intensity values...\n")
    mdat=preprocessRaw(rgSet)
    mdat=mdat[rownames(qcmat),]
    mdat=mdat[,colnames(qcmat)]
    #outliers based on total intensity values
    mu <- assays(mdat)$Meth+assays(mdat)$Unmeth
    mumean=apply(mu,2,mean,na.rm=TRUE)
    q2575 <- quantile(mumean, probs=c(0.25,0.75), na.rm=TRUE)
    qr <- q2575["75%"]-q2575["25%"]
    low=q2575["25%"] - 3*qr
    flag1=  mumean < low
    cat("After excluding low quality samples and CpGs\n")
    cat(sum(flag1)," samples are outliers based on averaged total intensity value","\n")

    #outliers in beta value distribution
    beta=getB(mdat, type="Illumina")
    qq=apply(beta,2,function(x) quantile(x, probs=c(0.25,0.5,0.75), na.rm=TRUE))
    q2575 <- apply(qq,1,function(x) quantile(x, probs=c(0.25,0.75), na.rm=TRUE))
    qr <- q2575["75%",]-q2575["25%",]
    up=q2575["75%",] + 3*qr
    low=q2575["25%",] - 3*qr
    flag=qq > up | qq < low
    flag=apply(flag,2,sum)>0
    cat(sum(flag)," samples are outliers in beta value distribution","\n")
    flag=flag | flag1
    badsample=c(badsample,colnames(beta)[flag])
    outlier_sample=colnames(beta)[flag]
    cat(sum(flag)," outlier samples were added into badsample list\n")
    if(ncol(qcmat)<15){
    cat("WARNING: Sample size may be too small to correctly identify outlier samples!\n")
    cat("RECOMMAND: set outlier=FALSE or double check total intensity and beta value 
     distribution plots to confirm\n")}
    }

    if(distplot)
    {
    mdat=preprocessRaw(rgSet)
    beta=getB(mdat, type="Illumina")

    cat("Ploting freqpolygon_beta_beforeQC.jpg ...")
    jpeg(filename="freqpolygon_beta_beforeQC.jpg",width=1000,
     height=500,quality = 100)
    color=rep("black",ncol(beta))
    color[colnames(beta) %in% badsample]="red"
    multifreqpoly(beta,cex.lab=1.4,cex.axis=1.5, col=color, legend="",
    cex.main=1.5,main="Beta value distribution",
    xlab="Methylation beta value")
    dev.off()
    cat("Done\n")

    beta=beta[!(rownames(beta) %in% badCpG),]
    beta=beta[,!(colnames(beta) %in% badsample)]
    cat("Ploting freqpolygon_beta_afterQC.jpg ...")
    jpeg(filename="freqpolygon_beta_afterQC.jpg",width=1000,
     height=500,quality = 100)
    multifreqpoly(beta,cex.lab=1.4,cex.axis=1.5, col="black",legend="",
    cex.main=1.5,main="Beta value distribution",
    xlab="Methylation beta value")
    dev.off()
    cat("Done\n")
    }

    if(outlier)
    {list(detP=detP,nbead=nbead,bisul=bisul,badsample=badsample,badCpG=badCpG,
      outlier_sample=outlier_sample)
    }else{list(detP=detP,nbead=nbead,bisul=bisul,badsample=badsample,badCpG=badCpG)}
}


