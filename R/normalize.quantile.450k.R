normalize.q <- function(x)
{
    if(sum(is.na(x))>0)
    {stop("Multi-array quantile normalization does not allow missing value")}
    rname <- rownames(x)
    cname <- colnames(x)
    x <- normalize.quantiles(x)
    rownames(x) <- rname
    colnames(x) <- cname
    x
}


normalize.quantile.450k <- function(mdat,method="quantile1")
{
    if(!is(mdat, "MethylSet")){stop("object needs to be of class 'MethylSet'")}
    cat("Analysis is running, please wait...!","\n")

    if (method == "quantile1")
    {
    anno <- getAnnotation(mdat)
    mdat_I <- mdat[anno$Type == "I",]
    mdat_II <- mdat[anno$Type == "II",]
    assays(mdat_I)$Meth<-normalize.q(assays(mdat_I)$Meth)
    assays(mdat_I)$Unmeth<-normalize.q(assays(mdat_I)$Unmeth)
    assays(mdat_II)$Meth<-normalize.q(assays(mdat_II)$Meth)
    assays(mdat_II)$Unmeth<-normalize.q(assays(mdat_II)$Unmeth)
    methData <- rbind(assays(mdat_I)$Meth,assays(mdat_II)$Meth)
    unmethData <- rbind(assays(mdat_I)$Unmeth,assays(mdat_II)$Unmeth)
    CpGID <- rownames(mdat)
    SampleID <- colnames(mdat)
    methData <- methData[CpGID,]
    unmethData <- unmethData[CpGID,]
    assays(mdat)$Meth<-methData
    assays(mdat)$Unmeth<-unmethData
    }
    else if(method == "quantile2")
    {
    anno <- getAnnotation(mdat)
    mdat_I <- mdat[anno$Type == "I",]
    mdat_II <- mdat[anno$Type == "II",]
    mat <- rbind(assays(mdat_I)$Meth,assays(mdat_I)$Unmeth)
    mat  <-  normalize.q(mat)
    assays(mdat_I)$Meth <- mat[1:(nrow(mat)/2),]
    assays(mdat_I)$Unmeth <- mat[((nrow(mat)/2)+1):nrow(mat),]
    mat <- rbind(assays(mdat_II)$Meth,assays(mdat_II)$Unmeth)
    mat  <-  normalize.q(mat)
    assays(mdat_II)$Meth <- mat[1:(nrow(mat)/2),]
    assays(mdat_II)$Unmeth <- mat[((nrow(mat)/2)+1):nrow(mat),]
    methData <- rbind(assays(mdat_I)$Meth,assays(mdat_II)$Meth)
    unmethData <- rbind(assays(mdat_I)$Unmeth,assays(mdat_II)$Unmeth)
    CpGID <- rownames(mdat)
    SampleID <- colnames(mdat)
    methData <- methData[CpGID,]
    unmethData <- unmethData[CpGID,]
    assays(mdat)$Meth<-methData
    assays(mdat)$Unmeth<-unmethData

    }
    else if(method == "quantile3")
    {
    ##quantile normalization of combined intensity values
    ##similar with lumi lumimethyN
    mat <- rbind(assays(mdat)$Meth,assays(mdat)$Unmeth)
    mat  <-  normalize.q(mat)
    methData <- mat[1:(nrow(mat)/2),]
    unmethData <- mat[((nrow(mat)/2)+1):nrow(mat),]
    assays(mdat)$Meth<-methData
    assays(mdat)$Unmeth<-unmethData
    }
    mdat@preprocessMethod <- c(mu.norm="normalize.quantile", 
    preprocessMethod(mdat))
    mdat
}
