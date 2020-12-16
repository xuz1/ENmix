QCfilter <-function(mdat,qcinfo=NULL,detPthre=0.000001,nbthre=3,
    samplethre=0.05,CpGthre=0.05,bisulthre=NULL,outlier=FALSE,
    outid=NULL, outCpG=NULL,plot=FALSE,rmlowq=FALSE)
{    
    if(!is(mdat,"methDataSet") & !is(mdat,"matrix")){
     stop("The input needs to be a methDataSet or a beta value matrix\n")}

    if(is.null(qcinfo)){stop("ERROR: Please provide QC information object
     from function QCinfo'")}
    if(outlier & is.null(qcinfo$outlier_sample)){
       stop("ERROR: No outlier sample information, please set
       outlier=FALSE\n")}
    qcmat<- qcinfo$nbead<nbthre | qcinfo$detP>detPthre
    bisul=qcinfo$bisul
    outlier_sample=qcinfo$outlier_sample

    rm(qcinfo)

    if(sum(!(colnames(mdat) %in% colnames(qcmat)))>0){stop("ERROR:
      some samples do not have QC information")}
    if(sum(!(rownames(mdat) %in% rownames(qcmat)))>0){stop("ERROR:
      some CpGs do not have QC information")}

    if(rmlowq & !is(mdat,"matrix"))
    {cat("Warning, input needs to be a data matrix to exclude low quality data points\n")
     cat("Reset rmlowq=FALSE\n"); rmlowq=FALSE}

    if(!identical(colnames(mdat),colnames(qcmat))){
    qcmat=qcmat[,colnames(mdat)]
    bisul=bisul[colnames(mdat)]
    outlier_sample=outlier_sample[outlier_sample %in% colnames(mdat)]
     }
    if(!identical(rownames(mdat),rownames(qcmat))){
    qcmat=qcmat[rownames(mdat),]
     }

    #threshold of bisulfite conversion control intensity
    if(is.null(bisulthre)){bisulthre=mean(bisul,na.rm=TRUE)-3*sd(bisul,
    na.rm=TRUE)}

    ##low quality samples
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

    ##outlier sample
    if(outlier & length(outlier_sample)>0){
    cat(length(outlier_sample)," outlier samples were excluded\n")
    badsample=unique(badsample,outlier_sample)
    }

    ##excluded user specified samples
    if(!is.null(outid)){
    cat(length(outid)," user specified samples were excluded\n")
    badsample=unique(badsample,outid)
    }

    ##excluded user specified CpG
    if(!is.null(outCpG)){
    cat(length(outCpG)," user specified CpGs were excluded\n")
    badCpG=unique(badCpG,outCpG)
    }

    ##summary
    cat("After excluding overlapped counts:\n")
    cat(length(badsample)," unique samples were excluded\n")
    cat(length(badCpG)," unique CpGs were excluded\n")


    if(plot)
    {
    #plotting quality scores
    cat("Ploting qc_sample.jpg ...")
    jpeg(filename="qc_sample.jpg",width=1000,height=1000,quality = 100)
    color=rep("black",length(badValuePerSample))
    color[flag]="red"
    plot(badValuePerSample,bisul,xlab="Percent of low quality data per sample",
       ylab="Average bisulfite conversion intensity",cex=1.5,col=color,
       main=paste(length(badsample)," samples were classified as low quality
       samples"))
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

    ##distribution plot before and after filtering
    if(is(mdat, "methDataSet")){
       mdat=getB(mdat, type="Illumina")
    }

    cat("Ploting freqpolygon_beta_beforeQC.jpg ...")
    jpeg(filename="freqpolygon_beta_beforeQC.jpg",width=1000,
     height=500,quality = 100)
    color=rep("black",ncol(mdat))
    color[colnames(mdat) %in% badsample]="red"
    multifreqpoly(mdat,cex.lab=1.4,cex.axis=1.5, col=color, legend="",
    cex.main=1.5,main="Beta value distribution",
    xlab="Methylation beta value")
    dev.off()
    cat("Done\n")

    mdat=mdat[!(rownames(mdat) %in% badCpG),]
    mdat=mdat[,!(colnames(mdat) %in% badsample)]
    cat("Ploting freqpolygon_beta_afterQC.jpg ...")
    jpeg(filename="freqpolygon_beta_afterQC.jpg",width=1000,
     height=500,quality = 100)
    multifreqpoly(mdat,cex.lab=1.4,cex.axis=1.5, col="black",legend="",
    cex.main=1.5,main="Beta value distribution",
    xlab="Methylation beta value")
    dev.off()
    cat("Done\n")
    }

    mdat=mdat[!(rownames(mdat) %in% badCpG),]
    mdat=mdat[,!(colnames(mdat) %in% badsample)]

    ## remove low quality data points
    if(rmlowq)
    {
        qcmat=qcmat[,colnames(mdat)]
        qcmat=qcmat[rownames(mdat),]

        mdat[qcmat]=NA
        cat(sum(qcmat), "low quality data points were removed\n")
    }
    mdat
}



