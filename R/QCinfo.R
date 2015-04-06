QCinfo <- function(rgSet)
{
    ##number of bead
    if(!is(rgSet, "RGChannelSetExtended"))
    stop("object needs to be of class 'RGChannelSetExtended'")

    nb <- assayDataElement(rgSet, "NBeads")
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

    #plot quanlity measures
    detPthre=0.05
    nbthre=3
    samplethre=0.05
    ##percent of low quality data per sample
    temp <- nbead<nbthre | detP>detPthre
    badValuePerSample <- apply(temp,2,sum)/nrow(temp)
    ##percent of low quality data per CpG
    ##remove low quality samples
    flag <- badValuePerSample>0.05
    temp <- temp[,!flag]
    badValuePerCpG <- apply(temp,1,sum)/ncol(temp)

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

    cat("Ploting qc_sample_1.jpg ...")
    jpeg(filename="qc_sample_1.jpg",width=1000,height=1000,quality = 100)
    summary(names(badValuePerSample)==names(bisul))
    plot(badValuePerSample,bisul,xlab="Percent of low quality data per sample",
       ylab="Average bisulfite conversion intensity",cex=1.5,main="")
    dev.off()
    cat("Done\n")
    cat("Ploting qc_sample_2.jpg ...")
    jpeg(filename="qc_sample_2.jpg",width=1000,height=1000,quality = 100)
    par(mfrow=c(2,1))
    hist(badValuePerSample,xlab="Percent of low quality data per sample",main="")
    hist(bisul,xlab="Average bisulfite conversion intensity",main="")
    dev.off()
    cat("Done\n")

    cat("Ploting qc_CpG.jpg ...")
    jpeg(filename="qc_CpG.jpg",width=1000,height=1000,quality = 100)
    par(mfrow=c(2,1))
    hist(badValuePerCpG,breaks=1000,xlab="Percent of low quality data per CpG",
         main="All CpG")
    hist(badValuePerCpG,breaks=1000,xlim=c(0,0.1),
       xlab="Percent of low quality data per CpG",main="Zoom in view")
    dev.off()
    cat("Done\n")

    list(detP=detP,nbead=nbead,bisul=bisul)
}
