QCfilter <-
function(mdat,qcinfo,detPthre=0.05,nbthre=3,samplethre=0.01,CpGthre=0.05,
         bisulthre=1,plot=FALSE, outid=NULL, outCpG=NULL)
{
    if(!is(mdat, "MethylSet")){stop("object needs to be of class 'MethylSet'")}
    if(sum(!rownames(qcinfo$nbead)==rownames(mdat))>0){stop("Error: different
	 matrix dimension, check qcinfo!")}
    if(sum(!colnames(qcinfo$nbead)==colnames(mdat))>0){stop("Error: different
	 matrix dimension, check qcinfo!")}

    ##density plot before QC
    if(plot)
    {
        cat("Ploting density_total_intensity_beforeQC.jpg ...")
        mu <- assayData(mdat)$Meth+assayData(mdat)$Unmeth
        jpeg(filename="density_total_intensity_beforeQC.jpg",width=1000,
             height=500,quality = 100)
        multidensity(mu,col="black",legend="",cex.lab=1.4,cex.axis=1.5,
                cex.main=1.5,main="Density plot of total intensity before QC",
                xlab="Total intensity")
        dev.off()
        cat("Done\n")
        cat("Ploting density_total_beta_beforeQC.jpg ...")
        beta <- getBeta(mdat, "Illumina")
        jpeg(filename="density_total_beta_beforeQC.jpg",width=1000,
             height=500,quality = 100)
        multidensity(beta,col="black",legend="",cex.lab=1.4,cex.axis=1.5,
                cex.main=1.5,main="Density plot of Beta value before QC",
                xlab="Methylation beta value")
        dev.off()
        cat("Done\n")
     }

    ##filter out low quanlity samples
    temp <- qcinfo$nbead<nbthre | qcinfo$detP>detPthre
    badValuePerSample <- apply(temp,2,sum)/nrow(temp)
    flag <- badValuePerSample > samplethre | qcinfo$bisul < bisulthre
    cat(sum(flag)," samples with percentage of low quanlity CpG value greater
	 than ",samplethre, " or bisulfite intensity less than ", bisulthre, "\n")
    cat("Low quality samples: ",colnames(mdat)[flag]," was excluded","\n")
    ##excluded user specified samples
    outid <- outid[outid %in% colnames(mdat)]
    outid <- outid[!(outid %in% colnames(mdat)[flag])]
    if(length(outid) > 0){
    cat("Samples excluded by user: ", outid, "\n")
    flag <- flag | (colnames(mdat) %in% outid)
    }
    mdat <- mdat[,!flag]
    ##filter out low quanlity CpGs
    temp <- temp[,!flag]
    NbadValuePerCpG <- apply(temp,1,sum)
    badValuePerCpG <- NbadValuePerCpG/ncol(temp)
    flag2 <- badValuePerCpG>CpGthre & NbadValuePerCpG>1
    cat(sum(flag2)," CpGs with percentage of low quanlity CpG value greater than ",
          CpGthre," were excluded","\n")
    ##exclude user specified CpGs
    outCpG <- outCpG[outCpG %in% rownames(mdat)]
    outCpG <- outCpG[!(outCpG %in% rownames(mdat)[flag2])]
    if(length(outCpG) > 0){
    cat("CpGs excluded by user:", outCpG, "\n")
    flag2 <- flag2 | (rownames(mdat) %in% outCpG)
    }
    mdat=mdat[!flag2,]

    ##density plot after QC
    if(plot)
    {
        cat("density_total_intensity_afterQC.jpg ...")
        mu <- assayData(mdat)$Meth+assayData(mdat)$Unmeth
        jpeg(filename="density_total_intensity_afterQC.jpg",width=1000,
             height=500,quality = 100)
        multidensity(mu,col="black",legend="",cex.lab=1.4,cex.axis=1.5,
                cex.main=1.5,main="Density plot of total intensity after QC",
                xlab="Total intensity")
        dev.off()
        cat("Done\n")
        
        cat("density_total_beta_afterQC.jpg ...")
        beta <- getBeta(mdat, "Illumina")
        jpeg(filename="density_total_beta_afterQC.jpg",width=1000,
             height=500,quality = 100)
        multidensity(beta,col="black",legend="",cex.lab=1.4,cex.axis=1.5,
                cex.main=1.5,main="Density plot of Beta value after QC",
                xlab="Methylation beta value")
        dev.off()
        cat("Done\n")
    }
    mdat
}
